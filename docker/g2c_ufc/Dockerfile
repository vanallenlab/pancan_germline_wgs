# The Germline Somatic Convergence (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0
# Germline_Somatic Docker
# Use an Ubuntu base image
FROM ubuntu:latest
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Install python packages dependencies
RUN apt-get update && \
    apt-get install -y python3-numpy python3-pandas \
                       python3-scipy python3-seaborn python3-qqman && \
    apt-get clean

# Install Python packages using pip
RUN python3 -m pip install numpy pandas qqman scipy seaborn statsmodels scikit-learn

# Copy your utility script to the image
COPY *.py /opt/

# Set the default command to run at runtime
CMD ["/bin/bash"]
