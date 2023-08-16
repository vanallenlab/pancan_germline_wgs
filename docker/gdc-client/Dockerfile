# The Germline Genomics of Cancer (G2C)
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Simple Dockerfile for basic image with NCI gdc-client installed

FROM ubuntu:focal-20220531
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Update libraries
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      less \
      curl \
      wget \
      file \
      unzip

# Install gdc-client
RUN cd /opt/ && \
    wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip && \
    unzip gdc-client_v1.6.1_Ubuntu_x64.zip && \
    rm gdc-client_v1.6.1_Ubuntu_x64.zip && \
    ./gdc-client --help && \
    cd /

# Launch bash at runtime
CMD ["/bin/bash"]
