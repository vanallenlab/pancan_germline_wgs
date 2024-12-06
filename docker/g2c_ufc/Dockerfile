# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for python3 w/ pandas

FROM ubuntu:latest
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"


# Install python packages dependencies
RUN apt-get update && \
    apt-get install -y python3-numpy python3-pandas \
                       python3-scipy python3-seaborn python3-pip && \
    apt-get clean

# Copy your utility script to the image
COPY *.py /opt/

# Install qqman via pip
RUN python3 -m pip3 install qqman

# Launch bash at runtime
CMD ["/bin/bash"]
