FROM ubuntu:18.04 as hedgehog-deps

# Install dependencies
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
# TODO COMPILE PETSC!!!
    

# ADD FFTW from mobo directory to the container and install 
#RUN wget http://www.fftw.org/fftw-3.3.4.tar.gz && \
#    tar xvf fftw-3.3.4.tar.gz && \
#    rm fftw-3.3.4.tar.gz && \
#    cd fftw-3.3.4 && \
#    ./configure --prefix=`pwd` --enable-openmp && \
#    make && \
#    make install
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


# Install p4est 
RUN wget http://p4est.github.io/release/p4est-1.1.tar.gz && \
    tar xvf p4est-1.1.tar.gz && \
    rm p4est-1.1.tar.gz && \
    cd /libs/p4est-1.1 && \
    ./configure CC=mpicc F77=mpif77 FC=mpif90 --enable-mpi && \
    make && \
    make install 

# Install blendsurf
RUN git clone https://github.com/mmorse1217/blendsurf.git &&\
    mkdir -p blendsurf/build/ && \
    cd blendsurf/build/ && \
    cmake ..  && \
    make  &&\
    make install

# Install patchwork
#RUN git clone https://github.com/mmorse1217/patchwork.git &&\
#    mkdir -p patchwork/build/ && \
#    cd patchwork/build/ && \
#    cmake -DCMAKE_MODULE_PATH=/usr/share/cmake-3.10/Modules/ ..  && \
#    make 
# TODO replace with subrepo + COPY + compile

ENV HEDGEHOG_LIBS=/libs \
    BLENDSURF_DIR=${HEDGEHOG_LIBS}/blendsurf \
    P4EST_DIR=${HEDGEHOG_LIBS}/p4est-1.1 \
    PVFMM_DIR=${HEDGEHOG_LIBS}/pvfmm \
    FFTW_DIR=${HEDGEHOG_LIBS}/fftw-3.3.4

CMD ["/bin/bash"]
FROM hedgehog-deps as hedgehog-build

#TEST try compiling ebi 
RUN mkdir -p ${HEDGEHOG_LIBS}/utils/pvfmm-utils
#RUN mkdir -p ${MOBO_DIR}/src/
ADD utils/pvfmm-utils ${HEDGEHOG_LIBS}/utils/pvfmm-utils
#ADD src ${MOBO_DIR}/src
#RUN ls ${MOBO_DIR}
CMD ["bash"]

