# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/g2c_pipeline:nf_graf_plink_alpha

FROM ubuntu:latest

MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"


# Install make, zlib, file, and git
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      default-libmysqlclient-dev \
      curl \
      vim \
      wget

# Install bcftools dependencies
RUN apt-get update && \
    apt-get install -y gcc make zlib1g-dev libbz2-dev liblzma-dev \
                       libcurl4-openssl-dev libssl-dev libncurses5-dev && \
    apt-get clean

# Install bcftools
ARG BCFTOOLS_VERISON=1.19
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xvf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    ./configure --without-curses && \
    make && \
    make install && \
    cd -

# Launch bash at runtime
CMD ["/bin/bash"]
