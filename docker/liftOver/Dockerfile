# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/liftover

#FROM ubuntu:22.04
#FROM ubuntu:focal-20220531
FROM continuumio/miniconda3
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"


# Install liftover dependencies
RUN apt-get update && \
    apt-get install -y gcc make zlib1g-dev libbz2-dev liblzma-dev \
                       libcurl4-openssl-dev libssl-dev libncurses5-dev && \
    apt-get clean

#Copy Liftover in
COPY liftOver /usr/local/bin/

# Launch bash at runtime
CMD ["/bin/bash"]
