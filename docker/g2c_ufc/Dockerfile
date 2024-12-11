# The Germline Somatic Convergence (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0
# Germline_Somatic Docker
# Use an Ubuntu base image

MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"

FROM condaforge/mambaforge:22.11.1-4 AS build
ENV CONDA_SUBDIR=linux-64
RUN conda config --env --set subdir linux-64 && \
    mamba update --all && \
    mamba install pandas seaborn matplotlib
RUN pip install qqman

# Copy your utility script to the image
COPY *.py /opt/

# Set the default command to run at runtime
CMD ["/bin/bash"]
