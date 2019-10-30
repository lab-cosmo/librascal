FROM ubuntu:18.04

# Install generic dependencies
RUN apt-get update && \
    apt-get install -y \
        software-properties-common \
        git \
        cmake \
        vim \
        python3.6 python3-pip \
        libboost-test-dev \
        doxygen pandoc

# - "ppa:ubuntu-toolchain-r/test" is used for gcc-9
RUN add-apt-repository ppa:ubuntu-toolchain-r/test

# Install compilers
RUN apt-get update && apt-get install -y gcc-9 g++-9
ENV CC=gcc-9 CXX=g++-9
