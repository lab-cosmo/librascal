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

# - "ppa:ubuntu-toolchain-r/test" is used for gcc-9
# - "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-9 main" is the
#   official LLVM repository for clang-9
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add -  && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    add-apt-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-9 main"

# Install compilers
RUN apt-get update && \
    apt-get install -y \
        gcc-5 g++-5 \
        gcc-9 g++-9 \
        clang-4.0 \
        clang-9
