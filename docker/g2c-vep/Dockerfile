# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Dockerfile for annotating VCFs with VEP

FROM ubuntu:focal-20220531
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install essential system tools and Perl environment
ENV DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_PRIORITY=critical
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      libperl-dev \
      libcurl4-openssl-dev \
      libssl-dev \
      libdbi-perl \
      libdbd-mysql-perl \
      default-libmysqlclient-dev \
      bioperl \
      curl \
      wget \
      file \
      git \
      perl \
      gzip \
      less \
      libmodule-build-perl \
      libarchive-zip-perl \
      libarchive-extract-perl \
      libtry-tiny-perl \
      libwww-perl \
      libscalar-util-numeric-perl \
      libset-intervaltree-perl \
      libdevel-cycle-perl \
      libpadwalker-perl \
      libtest-warn-perl \
      libipc-run-perl \
      cpanminus

# Install additional Perl modules needed for VEP + plugins
RUN cpanm --quiet --notest \
    Bio::Root::Version \
    Test::Warnings \
    PadWalker \
    Devel::Cycle \
    Set::IntervalTree \
    Scalar::Util::Numeric \
    Bio::EnsEMBL::Utils::IO

# Make VEP cache and plugin directories
RUN mkdir /opt/vep_cache /opt/vep_plugins
ENV VEP_CACHE=/opt/vep_cache
ENV VEP_PLUGINS=/opt/vep_plugins

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

# Clone Ensembl repos
RUN cd /opt && \
    git clone https://github.com/Ensembl/ensembl.git && \
    git clone https://github.com/Ensembl/ensembl-test.git && \
    git clone https://github.com/Ensembl/Bio-DB-HTS.git && \
    git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd -

ENV PERL5LIB=/opt/ensembl/modules:/opt/ensembl-test/modules:/opt/Bio-DB-HTS/lib

# Install Bio::DB::HTS
RUN cd /opt/Bio-DB-HTS && \
    perl Build.PL && \
    ./Build && \
    rm -f t/10leak.t t/11housekeeping_apache2.t && \
    ./Build test && \
    ./Build install && \
    perl -e "use Bio::DB::HTS" && \
    cd -

# Install VEP
ARG VEP_VERSION=113
RUN cd /opt/ensembl-vep && \
    perl INSTALL.pl \
      --NO_HTSLIB \
      --NO_TEST \
      --CACHE_VERSION $VEP_VERSION \
      --AUTO ap \
      --SPECIES homo_sapiens \
      --PLUGINS UTRAnnotator,dbNSFP,SpliceAI,AlphaMissense,PrimateAI \
      --CONVERT \
      --CACHEDIR $VEP_CACHE/ \
      --PLUGINSDIR $VEP_PLUGINS/ && \
    cd -

# Install LOFTEE (not available directly through VEP install)
RUN cd /opt && \
    git clone -b grch38 https://github.com/konradjk/loftee.git && \
    cp -rf loftee/* $VEP_PLUGINS/ && \
    cd -

# Alias VEP on path
ENV PATH=/opt/ensembl-vep/:$PATH

# Launch bash as entrypoint
CMD ["/bin/bash"]
