FROM continuumio/miniconda:4.5.11

# File Author / Maintainer
MAINTAINER Toni Hermoso Pulido <toni.hermoso@crg.eu> 

# https://anaconda.org/bioconda/ucsc-fetchchromsizes
# https://anaconda.org/bioconda/ucsc-oligomatch
COPY kentUtils.yml /
RUN conda env create -f /kentUtils.yml && conda clean -a
ENV PATH /opt/conda/envs/kentUtils/bin:$PATH

# Volume where to place input and output files
VOLUME /share

