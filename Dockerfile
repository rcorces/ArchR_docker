FROM rocker/rstudio:4.1.1
MAINTAINER paudocker02 <paulina.paiz@gladstone.ucsf.edu>

# Install system dependencies for R
RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    apt-transport-https \
    build-essential \
    gfortran \
    libatlas-base-dev \
    libbz2-dev \
    libcairo2-dev \
    libxml2-dev \
    # libcurl4-openssl-dev \
    libicu-dev \
    liblzma-dev \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libpcre3-dev \
    libtcl8.6 \
    libtiff5 \
    libtk8.6 \
    libx11-6 \
    libxt6 \
    libxt-dev \
    locales \
    tzdata \
    libglib2.0-dev \
    zlib1g-dev \
    meson \
    pkg-config \
    gtk-doc-tools \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libbz2-dev \
    libgsl0-dev 

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-dev \
    python3-pip

# Install system dependencies for the tidyverse R packages
# RUN apt-get install -y \
 #   make
    #libcurl4-openssl-dev
    #libssl-dev
    #pandoc
    # libxml2-dev
    
# Install Python
# FROM python:3.9.2-buster

# Install Python packages
RUN pip3 install MACS2

# Copy the current directory contents into the container at /app, instead of bind-mount or volume. 
# COPY . /app

# executed on build 
# RUN git clone https://github.com/paupaiz/ArchR
# WORKDIR /ArchR

ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
# RUN R -e Sys.setenv(RENV_DOWNLOAD_FILE_METHOD = getOption("download.file.method"))

WORKDIR /home/ArchR
COPY renv.lock renv.lock
# RENV_PATHS_CACHE_HOST=/opt/local/renv/cache
# RENV_PATHS_CACHE_CONTAINER=/renv/cache
# docker run --rm \
#    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
#    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
#    R --vanilla --slave -e 'renv::activate(); renv::restore()'
# will add source("renv/init.R") to the project .Rprofile, thereby instructing newly-launched R sessions to automatically load the current project.
# RUN R -e 'renv::activate()' 
# RUN R -e 'renv::isolate()'
# RUN R -e 'install.packages("devtools")'
# RUN R -e 'devtools::install_version("RcppAnnoy", "0.0.16", repos="http://cran.us.r-project.org")'
# specify BiocNeiighbors version
#  RUN R -e 'renv::install("bioc::BiocNeighbors")' 
RUN R -e 'renv::restore()'

# Get and install system dependencies
# RUN R -e "install.packages('remotes')" \
 # && R -e "remotes::install_github('r-hub/sysreqs')"

# FAILED RUN chown -R rstudio . \
# RUN R -e 'renv::restore()'
# RUN R -e "devtools::install_github("paupaiz/ArchR", ref="master", repos = BiocManager::repositories())"

# Default to starting R
# CMD ["R"]
