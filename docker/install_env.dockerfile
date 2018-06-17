FROM ubuntu:16.04

RUN apt-get update \
    && apt-get -y install \
    apt-utils \
    gcc \
    cmake \
    wget \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean all

# to build the image with the same uid and gid as the user
# docker build --build-arg UID=$(id -u) --build-arg GID=$(id -g) \
# -f Dockerfile -t test .

ARG UNAME=user
ARG UID=1000
ARG GID=1000
RUN groupadd -g $GID $UNAME
RUN useradd -m -u $UID -g $GID -s /bin/bash $UNAME
USER $UNAME
CMD /bin/bash

ENV HOME /home/$UNAME

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \ 
    && bash /tmp/miniconda.sh -bfp $HOME/ \ 
    && rm -rf /tmp/miniconda.sh \ 
    && conda install -y python=3 \ 
    && conda update conda \ 
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \ 
    && conda clean --all --yes 

ENV PATH $HOME/miniconda/bin:$PATH

RUN conda config --add channels conda-forge

# To get Eigen and pybind11
RUN conda install -c conda-forge -y git 

# To build the Tests
RUN conda install -c conda-forge -y boost
ENV Boost_NO_SYSTEM_PATHS=OFF 
ENV BOOST_ROOT=$HOME/miniconda

# To build the Doc
RUN conda install -c conda-forge -y sphinx breathe sphinx_rtd_theme doxygen 


# ## libboost-test1.58.0
# RUN apt-get update \
#     && apt-get -y install \
#     python3-pip \
#     libboost-test-dev \
#     && rm -rf /var/lib/apt/lists/* \
#     && apt-get clean all

# RUN pip3 install numpy
