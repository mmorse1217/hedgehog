FROM ubuntu:18.04 as hedgehog-deps

# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    build-essential \ 
    expat \
    extra-cmake-modules \
    git \
    gettext \
    libblas-dev \
    libexpat1-dev \
    liblapack-dev \
    libtool \
    libvtk7-dev \
    sudo \
    mpich \
    python \
    wget \
    zlib1g-dev &&\
    apt-get purge -y curl && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /libs 


# Install the latest version of cmake
RUN wget -O - https://raw.githubusercontent.com/mmorse1217/terraform/master/programs/cmake.sh | bash

#install Petsc 3.7.2
WORKDIR /libs
RUN mkdir petsc-3.7.2 && \
    wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.2.tar.gz &&\
    tar xvf petsc-3.7.2.tar.gz  && \
    rm petsc-3.7.2.tar.gz && \
    cd /libs/petsc-3.7.2 && \
    ./configure \
    --with-pic=true --with-shared-libraries=0 --prefix=/usr/local\
    --with-debugging=0 --with-64-bit-indices=1 --with-single-library=1 && \ 
    make && \
    make install 
    

# Install FFTW3: install single and double precision libraries for PVFMM
RUN wget http://www.fftw.org/fftw-3.3.4.tar.gz && \
    tar xvf fftw-3.3.4.tar.gz && \
    rm fftw-3.3.4.tar.gz && \
    cd fftw-3.3.4 && \
    ./configure --enable-type-prefix --with-pic --enable-openmp && \
    make && \
    make install && \
    make clean && \
    ./configure --enable-float --enable-type-prefix --with-pic  --enable-openmp && \
    make && \
    make install



# download and install PvFMM
RUN git clone https://github.com/dmalhotra/pvfmm && \
    mkdir -p /libs/pvfmm/build && \
    cd pvfmm/build && \
    cmake .. && \
    make && \
    make install
# TODO change precomputed directory to be the host machine not the container

# Install p4est 
RUN wget http://p4est.github.io/release/p4est-1.1.tar.gz && \
    tar xvf p4est-1.1.tar.gz && \
    rm p4est-1.1.tar.gz && \
    cd /libs/p4est-1.1 && \
    ./configure CC=mpicc CXX=mpic++ F77=mpif77 FC=mpif90 --enable-mpi --prefix=/usr/local && \
    make && \
    make install 

# Install blendsurf
RUN git clone https://github.com/mmorse1217/blendsurf.git &&\
    mkdir -p blendsurf/build/ && \
    cd blendsurf/build/ && \
    cmake ..  && \
    make  &&\
    make install


# Install Geogram
RUN cd /libs && \
    wget https://gforge.inria.fr/frs/download.php/file/37442/geogram_1.6.2.tar.gz &&\
    tar xvf geogram_1.6.2.tar.gz && \
    rm geogram_1.6.2.tar.gz && \
    cd geogram_1.6.2/ &&\
    mkdir build && \
    cd build &&\
    cmake -DGEOGRAM_WITH_VORPALINE=OFF -DVORPALINE_PLATFORM=Linux64-gcc\
    -DGEOGRAM_WITH_EXPLORAGRAM=OFF -DGEOGRAM_LIB_ONLY=ON\
    -DGEOGRAM_WITH_GRAPHICS=OFF .. &&\
    make &&\
    make install



RUN cd /libs && \
    git clone https://gitlab.com/libeigen/eigen.git && \
    cd eigen && \
    mkdir -p build/ && \
    cd build && \
    cmake .. && \
    make && \ 
    make install

# Install nanospline
RUN cd /libs && \
    git clone https://github.com/qnzhou/nanospline/ && \
    cd nanospline/ && \
    mkdir -p build/ && \
    cd build/ && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make 
    #make install

COPY patchwork/ /libs/patchwork

#Install patchwork via copy
RUN cd /libs/patchwork && \
    mkdir -p build/ && \
    cd build/ && \
    cmake ..  && \
    make && \
    make install

ENV CC=gcc CXX=g++ NANOSPLINE_DIR=/libs/nanospline


WORKDIR /hedgehog
CMD ["/bin/bash"]

FROM hedgehog-deps as hedgehog-dev

ENV VIM_DEV=1 DEBIAN_FRONTEND=noninteractive \
    PATH="~/miniconda3/bin:${PATH}"  \
    TERM=xterm-256color 
WORKDIR /terraform 

RUN git clone https://github.com/mmorse1217/terraform --recursive /terraform && \
    bash dotfiles/setup.sh && \
    bash programs/python.sh && \
    bash vim/build_from_source.sh && \  
    bash vim/lang-servers/setup.sh && \  
    bash vim/lang-servers/python-language-server.sh && \  
    bash vim/lang-servers/clangd.sh   && \
    bash vim/install_plugins.sh && \
    mkdir -p /hedgehog
#RUN apt-get update -y && apt install -y sudo git
#RUN bash dotfiles/setup.sh 
#
#RUN bash programs/python.sh 
#
## build vim
#RUN bash vim/build_from_source.sh  
#
## setup language servers for c++/python
#RUN bash vim/lang-servers/setup.sh  
#RUN bash vim/lang-servers/python-language-server.sh  
#RUN bash vim/lang-servers/clangd.sh  
#
## install plugins
#RUN bash vim/install_plugins.sh
#RUN mkdir /hedgehog
WORKDIR /hedgehog

CMD ["/bin/bash"]
