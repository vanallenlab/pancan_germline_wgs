# The Germline Genomics of Cancer (G2C)
# Contact: Ryan Collins <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/g2c_pipeline:nf_graf_plink_alpha

#FROM ubuntu:22.04
#FROM ubuntu:focal-20220531
FROM continuumio/miniconda3
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Create a conda environment and activate it
#RUN conda create --name myenv python=3.11 && \
#    conda activate myenv


# Install make, zlib, file, and git
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      default-libmysqlclient-dev \
      curl \
      vim \
      wget

#RUN apt install -y python3

# Install Plink on Mac
#RUN wget \
#      https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip && \
#      unzip plink_linux_x86_64_20231018.zip

# Print the current working directory
#RUN echo "Current working directory: $(pwd)" > log.txt

# Clone NCBI Grafpop
RUN wget -O GrafPop1.0.tar.gz \
      https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetZip.cgi?zip_name=GrafPop1.0.tar.gz && \
      tar -xzvf GrafPop1.0.tar.gz

# Move the grafpop executable
RUN mv grafpop /bin

# Copy the AncInferSNPs.txt to the home directory
RUN cp data/AncInferSNPs.txt .

# Remove all traces of hg19
RUN rm data/TG_10_1000_g37*

# Remove the tar.gz file for lighter Docker
RUN rm GrafPop1.0.tar.gz

COPY infer_ancestry.py /opt/
RUN chmod a+x /opt/infer_ancestry.py

# Launch bash at runtime
CMD ["/bin/bash"]
