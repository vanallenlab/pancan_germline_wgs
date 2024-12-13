# The Germline Somatic Convergence (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0
# Germline_Somatic Docker
# Use an Ubuntu base image
FROM ubuntu:latest
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Install Python 3.9 and pip
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.9 python3-pip python3.9-distutils && \
    apt-get clean

# Install Python packages using pip
RUN python3.9 -m pip install numpy pandas firthlogist scipy seaborn statsmodels scikit-learn

# Copy your utility script to the image
COPY gsc_util.py /opt/
COPY *.py /opt/

# Set the default command to run at runtime
CMD ["/bin/bash"]
