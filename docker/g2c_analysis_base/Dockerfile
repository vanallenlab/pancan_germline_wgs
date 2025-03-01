# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This Dockerfile builds the base of the image used for all analysis after initial
# variant calling (e.g., starting with quality control/filtering and proceeding
# through final formal analyses)

# Base image: conda-forge's miniconda build with mamba pre-installed (for building python environment prior to conda-pack)
FROM condaforge/mambaforge:24.9.2-0 AS build
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install git
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends git

# Clone pancan_germline_wgs repo and check out specified hash
ARG REPO_HASH=main
RUN git clone \
    https://github.com/vanallenlab/pancan_germline_wgs.git \
    opt/pancan_germline_wgs \
    --branch=$REPO_HASH

# Create default python environment (environment must be present in build context)
ENV CONDA_SUBDIR=linux-64
RUN conda config --env --set subdir linux-64 && \
    mamba update --all && \
    mamba env create --file=/opt/ped_germline_SV/config/environment.pedsv.yml

# Launch bash as entrypoint
CMD ["/bin/bash"]
