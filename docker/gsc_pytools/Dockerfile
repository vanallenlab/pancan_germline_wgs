# The Germline Somatic Convergence (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Use an Ubuntu base image
FROM ubuntu:latest
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Install Python 3.9 and dependencies
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.9 python3.9-distutils && \
    apt-get install -y python3.9-venv && \
    apt-get install -y python3-pip && \
    apt-get install -y python3-numpy python3-pandas python3-firthlogist \
                       python3-scipy python3-seaborn python3-statsmodels && \
    apt-get clean

# Upgrade pip to the latest version
RUN python3.9 -m pip install --upgrade pip

# Copy your utility script to the image
COPY gsc_util.py /opt/

# Set the default command to run at runtime
CMD ["/bin/bash"]
