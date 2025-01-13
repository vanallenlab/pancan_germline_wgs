# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/bcftools

#FROM ubuntu:22.04
#FROM ubuntu:focal-20220531
FROM continuumio/miniconda3
LABEL org.opencontainers.image.authors="Noah Fields <Noah_Fields@dfci.harvard.edu>"



COPY liftOver /opt/

# Launch bash at runtime
CMD ["/bin/bash"]
