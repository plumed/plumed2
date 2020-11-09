# Start from latest ubuntu
FROM centos:7

RUN yum -y update \
 && yum -y group install "Development Tools" \
 && yum -y install gawk vim zlib-devel gsl-devel fftw-devel

RUN useradd -ms /bin/bash plumed
USER plumed
WORKDIR /home/plumed

# import plumed code.
# assumes a file plumed2.tgz is present in the Dockerfile directory
COPY plumed2.tgz /home/plumed

# build and test plumed (no install)
RUN tar xzf plumed2.tgz \
 && cd plumed2 \
 && ./configure \
 && make -j 4 \
 && make check

