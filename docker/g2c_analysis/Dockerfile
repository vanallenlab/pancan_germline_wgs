# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This Dockerfile builds the main image used for all analysis after initial
# variant calling (e.g., starting with quality control/filtering and proceeding
# through final formal analyses)

# Bare bones ubuntu is used as base
ARG BASE_HASH=latest
FROM vanallenlab/g2c_analysis_base:$BASE_HASH AS build
FROM ubuntu:focal-20230301
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Essentially all dependencies, libraries, and environments are installed and
# configured in g2c_analysis_base
COPY --from=build /venv /venv

# Make bash default shell & load conda environment on container startup
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN echo "source /venv/bin/activate" >> ~/.bashrc
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install linux libraries
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      less \
      curl \
      wget \
      git

# Clone GATK-SV repo, check out specified hash, and only keep svtk
ARG GATKSV_TAG=v1.0.1
RUN git clone \
    https://github.com/broadinstitute/gatk-sv.git \
    opt/gatk-sv \
    --branch=$GATKSV_TAG && \
    mv opt/gatk-sv/src/svtk opt/ && \
    rm -rf opt/gatk-sv

# Clone RLCtools repo and check out specified hash
ARG RLCTOOLS_TAG=main
RUN git clone \
    https://github.com/RCollins13/RLCtools.git \
    opt/RLCtools \
    --branch=$RLCTOOLS_TAG

# Clone pancan_germline_wgs repo and check out specified hash
ARG G2C_TAG=main
RUN git clone \
    https://github.com/vanallenlab/pancan_germline_wgs.git \
    opt/pancan_germline_wgs \
    --branch=$G2C_TAG

# Activate virtual environment before installing helper libraries
RUN source activate g2c

# Install svtk
RUN pip install -e opt/svtk

# Install RLCtools
ARG RLCTOOLS_VERSION=0.1
RUN Rscript -e "install.packages('opt/RLCtools/RLCtools_$RLCTOOLS_VERSION.tar.gz', repos=NULL, type='source')"

# Install G2C helper libraries
ARG G2CR_VERSION=0.2.0
RUN Rscript -e "install.packages('opt/pancan_germline_wgs/src/G2CR_$G2CR_VERSION.tar.gz', repos=NULL, type='source')"
RUN pip install opt/pancan_germline_wgs/src/g2cpy

# Aliases for convenience
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Launch bash as entrypoint
CMD ["/bin/bash"]
