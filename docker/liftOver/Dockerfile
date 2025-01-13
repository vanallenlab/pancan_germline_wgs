# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/liftover

FROM ubuntu:22.04
#FROM ubuntu:focal-20220531
#FROM continuumio/miniconda3
LABEL org.opencontainers.image.authors="Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Update and install necessary tools to fix repository issues and install dependencies
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get update && \
    apt-get install -y \
        gnupg2 \
        curl \
        ca-certificates \
        libcurl4-openssl-dev \
        gcc make zlib1g-dev libbz2-dev liblzma-dev \
        libssl-dev libncurses5-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the liftOver binary to a directory in the PATH
COPY liftOver /usr/local/bin/liftOver

# Launch bash at runtime
CMD ["/bin/bash"]
