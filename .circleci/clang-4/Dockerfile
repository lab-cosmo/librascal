FROM ubuntu:18.04

# Install generic dependencies
RUN apt-get update && \
    apt-get install -y \
        git \
        cmake \
        vim \
        python3.6 python3-pip \
        libboost-test-dev \
        doxygen pandoc

# Install compilers
RUN apt-get update && apt-get install -y clang-4.0
ENV CC=clang-4.0 CXX=clang++-4.0
