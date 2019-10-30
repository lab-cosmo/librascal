FROM ubuntu:18.04

# Install generic dependencies
RUN apt-get update && \
    apt-get install -y \
        software-properties-common \
        wget \
        git \
        cmake \
        vim \
        python3.6 python3-pip \
        libboost-test-dev \
        doxygen pandoc

# - "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-9 main" is the
#   official LLVM repository for clang-9
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add -  && \
    add-apt-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-9 main"

# Install compilers
RUN apt-get update && apt-get install -y clang-9
ENV CC=clang-9 CXX=clang++-9
