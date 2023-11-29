# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This Dockerfile builds the primary image used for most core genomic data 
# processing, filtering, and quality control operations that rely on code from
# this git repo

# It is based on the primary sv-pipeline image from GATK-SV, which is available here:
# us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline
ARG BASE_IMAGE=us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2023-09-13-v0.28.3-beta-af8362e3
FROM $BASE_IMAGE as base

# Aliases for convenience
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Clone pancan_germline_wgs repo into image and check out specified hash
ARG REPO_HASH=main
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends git
RUN cd /opt && \
    git clone https://github.com/vanallenlab/pancan_germline_wgs.git && \
    cd pancan_germline_wgs && \
    git checkout $REPO_HASH

# Install missing R libraries as necessary
RUN Rscript -e 'install.packages("argparse", repos="https://cran.rstudio.com")'

# Launch bash as entrypoint
CMD ["/bin/bash"]
