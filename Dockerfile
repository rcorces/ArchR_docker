FROM ubuntu:xenial
# Install system dependencies for R
RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    apt-transport-https \
    build-essential \
    curl \
    gfortran \
    libatlas-base-dev \
    libbz2-dev \
    libcairo2 \
    libcurl4-openssl-dev \
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
    locales \
    tzdata \
    zlib1g-dev
    
# Install system dependencies for the tidyverse R packages
RUN apt-get install -y \
    make
    libcurl4-openssl-dev
    libssl-dev
    pandoc
    libxml2-dev
    
ARG R_VERSION=4.0.3
RUN wget https://cdn.rstudio.com/r/ubuntu-1604/pkgs/r-${R_VERSION}_1_amd64.deb
RUN apt-get install -y gdebi-core
RUN gdebi r-${R_VERSION}_1_amd64.deb

# pull in a manifest file and restore it
COPY renv.lock ./
RUN R -e 'renv::restore()'

    
