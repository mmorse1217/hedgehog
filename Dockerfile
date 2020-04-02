FROM ubuntu:18.04 as hedgehog-deps

# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    build-essential \ 
    cmake \ 
    expat \
    extra-cmake-modules \
    git \
    gettext \
    libblas-dev \
    libexpat1-dev \
    liblapack-dev \
    libtool \
    libvtk7-dev \
    mpich \
    python \
    wget \
    zlib1g-dev &&\
    apt-get purge -y curl && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /libs 



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





CMD ["/bin/bash"]
FROM hedgehog-deps as hedgehog-dev
RUN mkdir /hedgehog
WORKDIR /hedgehog
COPY patchwork/ /libs/patchwork

#Install patchwork
RUN cd /libs/patchwork && \
    mkdir -p build/ && \
    cd build/ && \
    cmake ..  && \
    make && \
    make install

#Install PVFMM wrapper
COPY utils/pvfmm-utils/ /libs/pvfmm-utils/
RUN cd /libs/pvfmm-utils/ && \
    mkdir build/ && \
    cd build/ && \
    cmake -DCMAKE_MODULE_PATH=cmake .. && \
    make VERBOSE=1 &&\
    make install

ENV CC=gcc CXX=g++
CMD ["bash"]


