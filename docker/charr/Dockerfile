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
      wget

#RUN apt install -y python3


# Clone NCBI Grafpop
RUN wget https://github.com/HTGenomeAnalysisUnit/SCE-VCF/releases/download/v0.1.2/sceVCF

#RUN chmod a+x sceVCF

# Launch bash at runtime
CMD ["/bin/bash"]
