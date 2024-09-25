# The Germline Somatic Convergence (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for python3 w/ pandas

FROM ubuntu:latest
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"


# Install python packages dependencies
RUN apt-get update && \
    apt-get install -y python3-numpy python3-pandas \
                       python3-scipy python3-seaborn && \
    apt-get clean

COPY gsc_util.py /opt/

#RUN pip3 install firthlogist

# Launch bash at runtime
CMD ["/bin/bash"]