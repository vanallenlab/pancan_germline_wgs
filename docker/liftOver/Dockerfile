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

RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver -O /usr/local/bin/liftOver

# Make liftOver executable
RUN chmod +x /usr/local/bin/liftOver

# Clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Launch bash at runtime
CMD ["/bin/bash"]
