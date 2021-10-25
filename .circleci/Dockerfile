FROM ubuntu:18.04

# Install generic dependencies
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
        software-properties-common \
        wget \
        git \
        cmake \
        libboost-test-dev \
        doxygen \
        pandoc \
        valgrind

# - "ppa:deadsnakes/ppa" provides other Python version
# - "ppa:ubuntu-toolchain-r/test" is used for gcc-10
# - "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-9 main" is the
#   official LLVM repository for clang-9
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add -  && \
    add-apt-repository ppa:deadsnakes/ppa && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    add-apt-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-9 main"
RUN apt-get update

# install Python 3.8
RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
        python3.8 \
        python3.8-dev \
        python3-pip

# set python 3.8 as the default so that it is picked up by cmake
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 2

# install pip in python 3.8
RUN python3.8 -m pip install --upgrade pip

# install compilers
RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
        gcc-5 g++-5 \
        gcc-10 g++-10 \
        clang-4.0 \
        clang-9 \
        clang-format-9
