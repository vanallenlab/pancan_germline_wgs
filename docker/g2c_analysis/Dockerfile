# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This Dockerfile builds the main image used for all analysis after initial
# variant calling (e.g., starting with quality control/filtering and proceeding
# through final formal analyses)

# Bare bones ubuntu is used as base
FROM ubuntu:focal-20230301 AS runtime
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Essentially all dependencies, libraries, and environments are installed and
# configured in g2c_analysis_base, which is used as the base image for this build
COPY --from=g2c_analysis_base /venv /venv

# Make bash default shell & load conda environment on container startup
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN echo "source /venv/bin/activate" >> ~/.bashrc
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install git
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends git

# Clone GATK-SV repo, check out specified hash, and only keep svtk
ARG SVTK_TAG=v1.0.1
RUN git clone \
    https://github.com/broadinstitute/gatk-sv.git \
    opt/gatk-sv \
    --branch=$SVTK_TAG && \
    mv opt/gatk-sv/src/svtk opt/ && \
    rm -rf opt/gatk-sv

# Install svtk
RUN pip install -e opt/svtk

# Clone RLCtools repo and check out specified hash
ARG RLCTOOLS_TAG=main
RUN git clone \
    https://github.com/RCollins13/RLCtools.git \
    opt/RLCtools \
    --branch=$RLCTOOLS_TAG

# Install RLCtools
ARG RLCTOOLS_VERSION=0.1
RUN Rscript -e 'install.packages("opt/RLCtools/RLCtools_$RLCTOOLS_VERSION.tar.gz", repos=NULL, type="source")'

# Clone pancan_germline_wgs repo and check out specified hash
ARG G2C_TAG=main
RUN git clone \
    https://github.com/vanallenlab/pancan_germline_wgs.git \
    opt/pancan_germline_wgs \
    --branch=$G2C_TAG

# Install G2C helper libraries
ARG G2CR_VERSION=0.2.0
RUN Rscript -e 'install.packages("opt/pancan_germline_wgs/src/G2CR/G2CR_$G2CR_version.tar.gz", repos=NULL, type="source")'
pip install code/src/g2cpy

# Aliases for convenience
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Launch bash as entrypoint
CMD ["/bin/bash"]
