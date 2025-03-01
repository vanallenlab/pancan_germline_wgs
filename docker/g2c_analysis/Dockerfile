# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This Dockerfile builds the main image used for all analysis after initial
# variant calling (e.g., starting with quality control/filtering and proceeding
# through final formal analyses)

# Essentially all dependencies, libraries, and environments are installed and
# configured in g2c_analysis_base, which is used as the base image for this build
FROM condaforge/mambaforge:24.9.2-0 AS build
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Aliases for convenience
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Install git
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends git

# Clone pancan_germline_wgs repo and check out specified hash
ARG REPO_HASH=main
RUN git clone \
    https://github.com/vanallenlab/pancan_germline_wgs.git \
    opt/pancan_germline_wgs \
    --branch=$REPO_HASH

# Install missing R libraries as necessary
RUN Rscript -e 'install.packages("argparse", repos="https://cran.rstudio.com")'

# Launch bash as entrypoint
CMD ["/bin/bash"]
