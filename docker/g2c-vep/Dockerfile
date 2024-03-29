# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Dockerfile for annotating VCFs with VEP

FROM ubuntu:focal-20220531
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install basic linux libraries
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      gzip \
      less \
      default-libmysqlclient-dev \
      curl \
      wget \
      file \
      git \
      perl

# Make VEP cache and plugin directories
RUN mkdir /opt/vep_cache /opt/vep_plugins
ENV VEP_CACHE=/opt/vep_cache
ENV VEP_PLUGINS=/opt/vep_plugins

# Install other libraries & perl modules required for VEP
ENV DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_PRIORITY=critical
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      bioperl \
      libbz2-dev \
      liblzma-dev \
      libmodule-build-perl \
      libdbi-perl \
      libdbd-mysql-perl \
      libarchive-zip-perl \
      libarchive-extract-perl \
      libtry-tiny-perl \
      libwww-perl \
      libscalar-util-numeric-perl \
      libset-intervaltree-perl

# Install htslib
ARG HTSLIB_VERSION=1.18
RUN cd /opt/ && \
    wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd -

# Install Bio::DB::HTS using preinstalled htslib (from prior step)
RUN cd /opt/ && \
    git clone https://github.com/Ensembl/Bio-DB-HTS.git && \
    cd Bio-DB-HTS && \
    perl Build.PL && \
    ./Build && \
    export PERL5LIB=$PERL5LIB:$(pwd -P)/lib:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/:$(pwd -P)/blib/arch/auto/Bio/DB/HTS/Faidx && \
    cd t && \
    for f in $(ls *.t); do perl $f; done && \
    perl -e "use Bio::DB::HTS" && \
    cd -
ENV PERL5LIB=$PERL5LIB:/opt/Bio-DB-HTS/lib:/opt/Bio-DB-HTS/blib/arch/auto/Bio/DB/HTS/:/opt/Bio-DB-HTS/blib/arch/auto/Bio/DB/HTS/Faidx

# Install VEP
ARG VEP_VERSION=110
RUN cd /opt && \
    git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    perl INSTALL.pl \
      --NO_HTSLIB \
      --NO_TEST \
      --CACHE_VERSION $VEP_VERSION \
      --AUTO ap \
      --SPECIES homo_sapiens \
      --PLUGINS UTRAnnotator,dbNSFP,SpliceAI \
      --CONVERT \
      --CACHEDIR $VEP_CACHE/ \
      --PLUGINSDIR $VEP_PLUGINS/ && \
    cd -

# Alias VEP on path
ENV PATH=/opt/ensembl-vep/:$PATH

# Launch bash as entrypoint
CMD ["/bin/bash"]
