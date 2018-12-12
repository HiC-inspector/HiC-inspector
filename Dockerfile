FROM biocorecrg/debian-perlbrew-pyenv:stretch

# File Author / Maintainer
MAINTAINER Toni Hermoso Pulido <toni.hermoso@crg.eu>

# Install Perl Modules
RUN cpanm Getopt::Long Pod::Usage Data::Dumper POSIX Benchmark Parallel::ForkManager

RUN apt-get update
RUN apt-get install -y gawk libxml2-dev libcurl4-openssl-dev

# Install R 3.5
RUN apt-get install -y software-properties-common apt-transport-https
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' &> /dev/null
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/' -y

RUN apt-get update
RUN apt-get install -y --allow-unauthenticated r-base r-base-dev

COPY deps.R /usr/local

RUN Rscript /usr/local/deps.R

ARG BOWTIE_VERSION=1.2.1.1
ARG BEDTOOLS_VERSION=2.27.1

# Install bowtie 
RUN cd /usr/local; curl --fail --silent --show-error --location --remote-name https://github.com/BenLangmead/bowtie/releases/download/v$BOWTIE_VERSION/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip
RUN cd /usr/local; unzip -d /usr/local bowtie-${BOWTIE_VERSION}-linux-x86_64.zip

RUN cd /usr/local; rm bowtie-${BOWTIE_VERSION}-linux-x86_64.zip

# Let's put in PATH
RUN cd /usr/local/bin; ln -s ../bowtie-${BOWTIE_VERSION}/bowtie* .


# Installing bedtools 
RUN cd /usr/local; curl --fail --silent --show-error --location --remote-name https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz
RUN cd /usr/local; tar zxf bedtools-${BEDTOOLS_VERSION}.tar.gz; cd /usr/local/bedtools2; make

RUN cd /usr/local; rm bedtools-${BEDTOOLS_VERSION}.tar.gz

# Let's put in PATH
RUN cd /usr/local/bin; cp -prf ../bedtools2/bin/* . ; rm -rf /usr/local/bedtools2

# Install fetchChromSizes
RUN cd /usr/local/bin; curl --fail --silent --show-error --location --remote-name http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes; chmod a+rx fetchChromSizes

# VOLUMES
VOLUME /utils
VOLUME /input
VOLUME /output

# Program
RUN mkdir /soft
COPY hic-inspector.pl /soft
COPY conf.pl /soft
COPY scripts /soft/scripts
RUN ln -s /utils /soft/utils
RUN ln -s /soft/hic-inspector.pl /usr/local/bin/hic-inspector.pl

# Clean
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*

# Let's place workdir in /soft
WORKDIR /soft

