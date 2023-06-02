# The Genomic Architecture of Human Cancers
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Simple Dockerfile for basic image with fastqsplitter installed

FROM continuumio/miniconda3:23.3.1-0
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install fastqsplitter
RUN conda update --all
RUN conda install -c bioconda fastqsplitter

# Launch bash at runtime
CMD ["/bin/bash"]
