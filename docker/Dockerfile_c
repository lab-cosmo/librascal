FROM ubuntu:16.04

RUN apt-get update \
    && apt-get -y install \
    build-essential \
    apt-utils \
    gcc \
    curl \
    cmake \
    wget \
    bzip2 \
    git \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean all


# RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \ 
#     && bash /tmp/miniconda.sh -bfp /usr/local \ 
#     && rm -rf /tmp/miniconda.sh \ 
#     && conda install -y python=3 \ 
#     && conda update conda \ 
#     && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \ 
#     && conda clean --all --yes 
    
# ENV PATH /opt/conda/bin:$PATH

# RUN conda config --add channels conda-forge
# RUN conda install -c conda-forge numpy git -y

# RUN conda install -c conda-forge sphinx breathe sphinx_rtd_theme doxygen -y
## libboost-test1.58.0
RUN apt-get update \
    && apt-get -y install \
    python3-pip \ 
    libboost-test-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean all

RUN pip3 install numpy