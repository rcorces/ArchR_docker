library(ArchR)
library(Seurat)
library(Cairo) 
set.seed(1)

addArchRThreads(threads = 1) 
inputFiles <- getTutorialData("Hematopoiesis")
addArchRGenome("hg19")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

# paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")

# getAvailableMatrices(projHeme1)

bioNames <- gsub("_R2|_R1|scATAC_","",projHeme1$Sample)
projHeme1$bioNames <- bioNames
bioNames <- bioNames[1:10]
cellNames <- projHeme1$cellNames[1:10]
projHeme1 <- addCellColData(ArchRProj = projHeme1, data = paste0(bioNames),
                            cells = cellNames, name = "bioNames2")

# plot quality control metrics
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

# plot sample statistics
# make ridge plot 
p1 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p2 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p4 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)

# Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

#added from short tutorial
projHeme1 <- addIterativeLSI(ArchRProj = projHeme1, useMatrix = "TileMatrix", name = "IterativeLSI")
projHeme1 <- addClusters(input = projHeme1, reducedDims = "IterativeLSI")
projHeme1 <- addImputeWeights(projHeme1) 

# saving and loading an ArchR project
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)

projHeme2 <- filterDoublets(projHeme1)

projHeme2 <- addIterativeLSI(
  ArchRProj = projHeme2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

# batch effect correction with harmony 
projHeme2 <- addHarmony(
  ArchRProj = projHeme2,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# clustering using Seurat
projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)
table(projHeme2$Clusters)
cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))

# clustering using scran
projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI",
  method = "scran",
  name = "ScranClusters",
  k = 15
)

###### UMAP ######
projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  #name = "UMAP", 
  #nNeighbors = 30, 
  #minDist = 0.5, 
  #metric = "cosine"
)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
#plot results using scran
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-ScranClusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

### t-SNE
projHeme2 <- addTSNE(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  name = "TSNE", 
  perplexity = 30
)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

### Dim reduction after Harmony ####
projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

projHeme2 <- addTSNE(
  ArchRProj = projHeme2, 
  reducedDims = "Harmony", 
  name = "TSNEHarmony", 
  perplexity = 30
)
p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-TSNE2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

### Identifying Marker Genes
markersGS <- getMarkerFeatures(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
  "CD14", "CEBPB", "MPO", #Monocytes
  "IRF8", 
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

### Visualizing marker genes on an Embedding ### 
p <- plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)
### Marker Genes Imputation with MAGIC
projHeme2 <- addImputeWeights(projHeme2)
markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)

p <- plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme2)
)
#Rearrange for grid plotting
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)

### Track Plotting with ArchRBrowser
markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", #B-Cell Trajectory
  "CD14", #Monocytes
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

p <- plotBrowserTrack(
  ArchRProj = projHeme2, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)

### Defining Cluster Identity with scRNA-seq
if(!file.exists("scRNA-Hematopoiesis-Granja-2019.rds")){
  download.file(
    url = "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/scRNA-Hematopoiesis-Granja-2019.rds",
    destfile = "scRNA-Hematopoiesis-Granja-2019.rds"
  )
}

seRNA <- readRDS("scRNA-Hematopoiesis-Granja-2019.rds")

### unconstrained integration
projHeme2 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "BioClassification",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
### constrained integration
cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
#From scRNA
cTNK <- paste0(paste0(19:25), collapse="|")
cNonTNK <- paste0(c(paste0("0", 1:9), 10:13, 15:18), collapse="|")
#Assign scATAC to these categories
clustTNK <- rownames(cM)[grep(cTNK, preClust)]
clustNonTNK <- rownames(cM)[grep(cNonTNK, preClust)]
#RNA get cells in these categories
rnaTNK <- colnames(seRNA)[grep(cTNK, colData(seRNA)$BioClassification)]
rnaNonTNK <- colnames(seRNA)[grep(cNonTNK, colData(seRNA)$BioClassification)]
groupList <- SimpleList(
  TNK = SimpleList(
    ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustTNK],
    RNA = rnaTNK
  ),
  NonTNK = SimpleList(
    ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustNonTNK],
    RNA = rnaNonTNK
  )    
)
projHeme2 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE, 
  groupList = groupList,
  groupRNA = "BioClassification",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co"
)
### Comparing Unconstrained and Constrained Integrations
pal <- paletteDiscrete(values = colData(seRNA)$BioClassification)
p1 <- plotEmbedding(
  projHeme2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal
)
p2 <- plotEmbedding(
  projHeme2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Co", 
  pal = pal
)
plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2", load = FALSE)

### Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
projHeme3 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "BioClassification",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)
projHeme3 <- addImputeWeights(projHeme3)
markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", #B-Cell Trajectory
  "CD14", #Monocytes
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

p1 <- plotEmbedding(
  ArchRProj = projHeme3, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme3)
)

p2 <- plotEmbedding(
  ArchRProj = projHeme3, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme3)
)

p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

plotPDF(plotList = p1, 
        name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = projHeme3, 
        addDOC = FALSE, width = 5, height = 5)

cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]

remapClust <- c(
  "01_HSC" = "Progenitor",
  "02_Early.Eryth" = "Erythroid",
  "03_Late.Eryth" = "Erythroid",
  "04_Early.Baso" = "Basophil",
  "05_CMP.LMPP" = "Progenitor",
  "06_CLP.1" = "CLP",
  "07_GMP" = "GMP",
  "08_GMP.Neut" = "GMP",
  "09_pDC" = "pDC",
  "10_cDC" = "cDC",
  "11_CD14.Mono.1" = "Mono",
  "12_CD14.Mono.2" = "Mono",
  "13_CD16.Mono" = "Mono",
  "15_CLP.2" = "CLP",
  "16_Pre.B" = "PreB",
  "17_B" = "B",
  "18_Plasma" = "Plasma",
  "19_CD8.N" = "CD8.N",
  "20_CD4.N1" = "CD4.N",
  "21_CD4.N2" = "CD4.N",
  "22_CD4.M" = "CD4.M",
  "23_CD8.EM" = "CD8.EM",
  "24_CD8.CM" = "CD8.CM",
  "25_NK" = "NK"
)
remapClust <- remapClust[names(remapClust) %in% labelNew]
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)
p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2")
plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

### Pseudo-bulk Replicates in ArchR
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2")

### calling peak in ArchR
pathToMacs2 <- findMacs2()
projHeme4 <- addReproduciblePeakSet(
  ArchRProj = projHeme4, 
  groupBy = "Clusters2", 
  pathToMacs2 = pathToMacs2
)
projHemeTmp <- addReproduciblePeakSet(
  ArchRProj = projHeme4, 
  groupBy = "Clusters2",
  peakMethod = "Tiles",
  method = "p"
)
saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4", load = FALSE)
projHeme5 <- addPeakMatrix(projHeme4)

### Identifying Marker Peaks with ArchR
#Our scRNA labels
table(projHeme5$Clusters2)
markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
### Plotting Marker Peaks in ArchR
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
pma <- markerPlot(seMarker = markersPeaks, name = "Erythroid", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markersPeaks, name = "Erythroid", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "Erythroid-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
p <- plotBrowserTrack(
  ArchRProj = projHeme5, 
  groupBy = "Clusters2", 
  geneSymbol = c("GATA1"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
  upstream = 50000,
  downstream = 50000
)
plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

### Pairwise Testing Between Groups
markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Erythroid",
  bgdGroups = "Progenitor"
)
pma <- markerPlot(seMarker = markerTest, name = "Erythroid", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Erythroid", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "Erythroid-vs-Progenitor-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

### Motif and Feature Enrichment
projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
plotPDF(ggUp, ggDo, name = "Erythroid-vs-Progenitor-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

### Motif Enrichment in Marker Peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme5,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)

### Bulk ATAC-seq
projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "ATAC")
enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme5,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "Codex")
enrichCodex <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme5,
  peakAnnotation = "Codex",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
plotPDF(heatmapCodex, name = "Codex-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)

### ChromVAR Deviatons Enrichment with ArchR
if("Motif" %ni% names(projHeme5@peakAnnotation)){
  projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
}
projHeme5 <- addBgdPeaks(projHeme5)
projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
p <- plotGroups(ArchRProj = projHeme5, 
                groupBy = "Clusters2", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(projHeme5)
)
p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

### footprinting with archr
motifPositions <- getPositions(projHeme5)
motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2")
seFoot <- getFootprints(
  ArchRProj = projHeme5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2")
seTSS <- getFootprints(
  ArchRProj = projHeme5, 
  positions = GRangesList(TSS = getTSS(projHeme5)), 
  groupBy = "Clusters2",
  flank = 2000
)

### Co-accessibility with ArchR
projHeme5 <- addCoAccessibility(
  ArchRProj = projHeme5,
  reducedDims = "IterativeLSI"
)
cA <- getCoAccessibility(
  ArchRProj = projHeme5,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA <- getCoAccessibility(
  ArchRProj = projHeme5,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = TRUE
)
cA <- getCoAccessibility(
  ArchRProj = projHeme5,
  corCutOff = 0.5,
  resolution = 1000,
  returnLoops = TRUE
)
cA <- getCoAccessibility(
  ArchRProj = projHeme5,
  corCutOff = 0.5,
  resolution = 10000,
  returnLoops = TRUE
)

markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", #B-Cell Trajectory
  "CD14", #Monocytes
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

p <- plotBrowserTrack(
  ArchRProj = projHeme5, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(projHeme5)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
        ArchRProj = projHeme5, 
        addDOC = FALSE, width = 5, height = 5)

projHeme5 <- addPeak2GeneLinks(
  ArchRProj = projHeme5,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = projHeme5,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g <- getPeak2GeneLinks(
  ArchRProj = projHeme5,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE
)

p2g <- getPeak2GeneLinks(
  ArchRProj = projHeme5,
  corCutOff = 0.45,
  resolution = 1000,
  returnLoops = TRUE
)

sap2g <- getPeak2GeneLinks(
  ArchRProj = projHeme5,
  corCutOff = 0.45,
  resolution = 10000,
  returnLoops = TRUE
)

p <- plotPeak2GeneHeatmap(ArchRProj = projHeme5, groupBy = "Clusters2")

seGroupMotif <- getGroupSE(ArchRProj = projHeme5, useMatrix = "MotifMatrix", groupBy = "Clusters2")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corGSM_MM <- correlateMatrices(
  ArchRProj = projHeme5,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)
corGIM_MM <- correlateMatrices(
  ArchRProj = projHeme5,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

### Trajectory Analysis with ArchR
trajectory <- c("Progenitor", "GMP", "Mono")
projHeme5 <- addTrajectory(
  ArchRProj = projHeme5, 
  name = "MyeloidU", 
  groupBy = "Clusters2",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)
p <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "cellColData", name = "MyeloidU")
plotPDF(p, name = "Plot-MyeloidU-Traj-UMAP.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 5, height = 5)

projHeme5 <- addImputeWeights(projHeme5)
p1 <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "GeneScoreMatrix", name = "CEBPB", continuousSet = "horizonExtra")
p2 <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "GeneIntegrationMatrix", name = "CEBPB", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
trajMM  <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "MotifMatrix", log2Norm = FALSE)

p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "Plot-MyeloidU-Traj-Heatmaps.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 6, height = 8)
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
