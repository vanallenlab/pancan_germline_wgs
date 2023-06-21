# The Germline Genomics of Cancer (G2C)
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Simple Dockerfile for basic image with fastq-pair installed

FROM ubuntu:focal-20220531
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install make, zlib, file, and git
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      less \
      default-libmysqlclient-dev \
      curl \
      wget \
      file \
      git

# Install cmake
RUN rm -rf /var/lib/apt/lists/* && \
    wget \
      https://github.com/Kitware/CMake/releases/download/v3.24.1/cmake-3.24.1-Linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh && \
    chmod u+x /tmp/cmake-install.sh && \
    mkdir /opt/cmake-3.24.1 && \
    /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake-3.24.1 && \
    rm /tmp/cmake-install.sh && \
    ln -s /opt/cmake-3.24.1/bin/* /usr/local/bin
RUN cmake --version

# Install fastq-pair
RUN cd /opt && \
    git clone https://github.com/linsalrob/fastq-pair.git && \
    cd fastq-pair && \
    mkdir build && cd build && \
    cmake .. && \
    make && make install && \
    cd /

# Launch bash at runtime
CMD ["/bin/bash"]
