# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/liftover

#FROM ubuntu:24.04
#FROM ubuntu:focal-20220531
FROM continuumio/miniconda3
LABEL org.opencontainers.image.authors="Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Install bcftools dependencies
RUN apt-get update -y
RUN apt-get install -y gcc make zlib1g-dev libbz2-dev liblzma-dev \
                       libcurl4-openssl-dev libssl-dev libncurses5-dev
RUN apt-get clean

# Copy the liftOver binary to a directory in the PATH
COPY liftOver /usr/local/bin/liftOver

# Launch bash at runtime
CMD ["/bin/bash"]
