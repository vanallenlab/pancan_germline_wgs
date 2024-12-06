# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for python3 w/ pandas

FROM ubuntu:latest
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Install python packages dependencies
RUN apt-get update && \
    apt-get install -y python3-pip && \
    apt-get clean

# Install Python packages using pip
RUN python3 -m pip install numpy pandas qqman

# Set up a virtual environment in /opt/myenv
#RUN python3 -m venv /opt/myenv

# Install qqman and all its dependencies inside the virtual environment
#RUN /opt/myenv/bin/pip install qqman

# Set environment variables so the virtual environment is used by default
#ENV PATH="/opt/myenv/bin:$PATH"

# Copy your utility script to the image
COPY *.py /opt/

# Launch bash at runtime
CMD ["/bin/bash"]
