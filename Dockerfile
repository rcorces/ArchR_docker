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
    libcairo2-dev \
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
    libgsl0-dev \
    gsl-bin
        
# Install system dependencies for the tidyverse R packages
RUN apt-get install -y \
    make
    libcurl4-openssl-dev
    libssl-dev
    pandoc
    libxml2-dev
    
# Install Python
FROM python:3.9.2-buster

# Install Python packages
RUN pip install MACS2

# Install GSL
#RUN wget https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
#RUN tar -zxvf gsl-2.6.tar.gz; cd gsl-2.6; ./configure && make && make install

# Install R
ARG R_VERSION=4.0.3
FROM r-base:${R_VERSION}

# Install R packages
# pull in an renv manifest file and restore it to load pacakges
COPY renv.lock ./
RUN R -e 'renv::restore()'

    
