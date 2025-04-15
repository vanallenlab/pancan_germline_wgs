FROM condaforge/mambaforge:22.11.1-4 AS build
ENV CONDA_SUBDIR=linux-64
RUN conda config --env --set subdir linux-64 && \
    mamba install pandas seaborn matplotlib

RUN pip install networkx
RUN pip install qqman

# Copy your utility script to the image
COPY *.py /opt/