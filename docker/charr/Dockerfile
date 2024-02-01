# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/g2c_pipeline:nf_graf_plink_alpha

#FROM ubuntu:22.04
#FROM ubuntu:focal-20220531
FROM continuumio/miniconda3
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Create a conda environment and activate it
#RUN conda create --name myenv python=3.11 && \
#    conda activate myenv


# Install make, zlib, file, and git
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      default-libmysqlclient-dev \
      curl \
      vim \
      wget \
      rm

#RUN apt install -y python3


# Clone NCBI Grafpop
RUN wget https://github.com/HTGenomeAnalysisUnit/SCE-VCF/releases/download/v0.1.2/sceVCF

RUN chmod +x sceVCF
RUN mv sceVCF /bin/sceVCF

RUN wget -O /bin/echtvar https://github.com/brentp/echtvar/releases/download/v0.1.9/echtvar \
    && chmod +x /bin/echtvar

# Install bcftools dependencies
RUN apt-get update && \
    apt-get install -y gcc make zlib1g-dev libbz2-dev liblzma-dev \
                       libcurl4-openssl-dev libssl-dev libncurses5-dev && \
    apt-get clean

# Install bcftools
#ARG BCFTOOLS_VERISON=1.19
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xvf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    ./configure --prefix=/ && \
    make && \
    make install

RUN rm -r bcftools-1.19
RUN rm bcftools-1.19.tar.bz2

# Launch bash at runtime
CMD ["/bin/bash"]
