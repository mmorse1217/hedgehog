FROM ubuntu:14.04

# Install dependecies
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    build-essential \
    expat \
    freeglut3-dev \ 
    g++-4.7 \ 
    gcc-4.7 \
    gettext \
    gfortran-4.7 \ 
    gfortran-4.4=4.4.7-8ubuntu1 \
    git \
    libblas-dev \
    libexpat1-dev \
    liblapack-dev \
    libtool \
    make=3.81-8.2ubuntu3 \
    mercurial \
    openssh-client \
    openssh-server \
    tar \
    vim

WORKDIR /
RUN mkdir hedgehog-libs 
RUN cd hedgehog-libs && pwd
ENV HEDGEHOG_LIBS /hedgehog-libs

#install MPICH 1.2.1

# Download and unzip
WORKDIR ${HEDGEHOG_LIBS}
RUN wget http://www.mpich.org/static/downloads/1.2.1/mpich2-1.2.1.tar.gz
RUN tar xvf mpich2-1.2.1.tar.gz 

WORKDIR ${HEDGEHOG_LIBS}/mpich2-1.2.1

# Configure and install

RUN ./configure CC=gcc-4.7 CXX=g++-4.7 F77=gfortran-4.7 F90=gfortran-4.7 --enable-sharedlibs=gcc --prefix=${HEDGEHOG_LIBS}/mpich2-1.2.1
RUN make
RUN make install -i
# ^^^ install with -i due to an error with make install that causes libmpich.a to not get copied :(
# add it to PATH so we can use mpicc/mpif90


ENV PATH ${HEDGEHOG_LIBS}/mpich2-1.2.1/bin:$PATH
ENV MPI_HOME ${HEDGEHOG_LIBS}/mpich2-1.2.1

# setup mpd, for whatever reason
RUN cd /etc && touch /etc/mpd.conf && echo "MPD_SECRETWORD=ooo_so_secret" > /etc/mpd.conf && chmod 600 /etc/mpd.conf


#install Petsc 3.7.2
WORKDIR ${HEDGEHOG_LIBS}
RUN mkdir petsc-3.7.2
RUN wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.2.tar.gz 
RUN tar xvf petsc-3.7.2.tar.gz 

WORKDIR ${HEDGEHOG_LIBS}/petsc-3.7.2
RUN ./configure PETSC_ARCH=linux-gnu PETSC_DIR=${HEDGEHOG_LIBS}/petsc-3.7.2 --with-debugging=0 --with-64-bit-indices=1 --with-single-library=1 --download-fblaslaplack --with-mpi-dir=$MPI_HOME

RUN (mpd &) && (make all test)
# ADD FFTW from mobo directory to the container and install 
WORKDIR ${HEDGEHOG_LIBS}
RUN wget http://www.fftw.org/fftw-3.3.4.tar.gz
RUN tar xvf fftw-3.3.4.tar.gz

WORKDIR ${HEDGEHOG_LIBS}/fftw-3.3.4
RUN ./configure --prefix=`pwd`  --enable-openmp
RUN make
RUN make install
ENV FFTW_DIR ${HEDGEHOG_LIBS}/fftw-3.3.4

# download and install PvFMM
WORKDIR ${HEDGEHOG_LIBS}
RUN git clone https://github.com/dmalhotra/pvfmm
WORKDIR ${HEDGEHOG_LIBS}/pvfmm
# last working version
RUN git checkout 7ee406618aa7e6d17480658ac462fb2c88fe427d
RUN autoreconf --install
RUN ./configure --with-fftw=${FFTW_DIR} CC=gcc-4.7 MPICC=mpicc F77=gfortran-4.7 --prefix=`pwd`
RUN make
RUN make install

ENV PVFMM_DIR ${HEDGEHOG_LIBS}/pvfmm

# Install p4est
WORKDIR ${HEDGEHOG_LIBS}
RUN wget http://p4est.github.io/release/p4est-1.1.tar.gz
RUN tar xvf p4est-1.1.tar.gz
WORKDIR ${HEDGEHOG_LIBS}/p4est-1.1
RUN ./configure CC=mpicc F77=mpif77 FC=mpif90 --prefix=`pwd` --enable-mpi
RUN make
RUN make install 

# Install blendsurf...
RUN mkdir ${HEDGEHOG_LIBS}/blendsurf3
COPY blendsurf3/ ${HEDGEHOG_LIBS}/blendsurf3
WORKDIR ${HEDGEHOG_LIBS}/blendsurf3
RUN make clean -C ${HEDGEHOG_LIBS}/blendsurf3
ENV HEDGEHOG_CC gcc
ENV HEDGEHOG_CXX gcc
RUN make -C ${HEDGEHOG_LIBS}/blendsurf3

ENV BLENDSURF_DIR ${HEDGEHOG_LIBS}/blendsurf3
ENV FACEMAP_DIR ${HEDGEHOG_LIBS}/face_map
RUN apt-get update
RUN apt-get install -y libpng-dev


# Build Face-map...
#ENV P4EST_DIR ${HEDGEHOG_LIBS}/p4est-1.1
ENV PETSC_ARCH linux-gnu
ENV PETSC_DIR ${HEDGEHOG_LIBS}/petsc-3.7.2
#ENV DOCKER 1
#ENV MACHINE_NAME docker 
#
#RUN mkdir -p ${HEDGEHOG_LIBS}/face_map
#COPY face_map/ ${HEDGEHOG_LIBS}/face_map 
#WORKDIR ${HEDGEHOG_LIBS}/face_map
#RUN make clean -C ${HEDGEHOG_LIBS}/face_map
#RUN make clean -C ${HEDGEHOG_LIBS}/face_map/tests
#ENV PATH /usr/lib:$PATH
RUN ls ${PETSC_DIR}/${PETSC_ARCH}/include
RUN ls ${PETSC_DIR}/${PETSC_ARCH}/lib
#RUN ls ${MPI_HOME}/lib
#RUN ls ${MPI_HOME}/include
#RUN make -j2 -C ${HEDGEHOG_LIBS}/face_map lib
#RUN make -j2 -C ${HEDGEHOG_LIBS}/face_map/tests 

## TODO fix linking order in Face-map
# set some Mobo-related environment variables
#ENV MOBO_DIR /hedgehog
#ENV MACHINE_NAME energon

#TEST try compiling ebi 
RUN mkdir -p ${HEDGEHOG_LIBS}/utils/pvfmm-utils
#RUN mkdir -p ${MOBO_DIR}/src/
ADD utils/pvfmm-utils ${HEDGEHOG_LIBS}/utils/pvfmm-utils
#ADD src ${MOBO_DIR}/src
#RUN ls ${MOBO_DIR}
#RUN ls ${MOBO_DIR}/pvfmm/pvfmm-utils
#RUN ls ${MOBO_DIR}/src/ebi
# install CMake
WORKDIR /
RUN wget https://github.com/Kitware/CMake/releases/download/v3.13.4/cmake-3.13.4-Linux-x86_64.sh 
RUN yes Y | sh cmake-3.13.4-Linux-x86_64.sh --prefix /usr/local --exclude-subdir
# Compile pvfmm wrapper
#ENV PETSC_ROOT ${PETSC_DIR} 
#WORKDIR ${HEDGEHOG_LIBS}/utils/pvfmm-utils
#RUN make clean
#RUN make -j2 lib 
#RUN make -j2 pvfmm_test 

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#WORKDIR ${MOBO_DIR}/src/ebi/
#
#ENV PVFMM_INC ${PVFMM_DIR}/include
#ENV PVFMM_LIB ${PVFMM_DIR}/lib
#ENV DOCKER 1
#WORKDIR ${MOBO_DIR}/src/ebi
#ENV PVFMM_INC ${PVFMM_DIR}/include
#ENV PVFMM_LIB ${PVFMM_DIR}/lib 
RUN rm -f ${HEDGEHOG_LIBS}/fftw-3.3.4.tar.gz ${HEDGEHOG_LIBS}/mpich2-1.2.1.tar.gz ${HEDGEHOG_LIBS}/petsc-3.7.2.tar.gz ${HEDGEHOG_LIBS}/p4est-1.1.tar.gz
#RUN apt-get install -y vim
##ENV FORTRAN_LIB /usr/lib/gcc/x86_64-linux-gnu/4.7/
#RUN echo 'alias mobo_makegrep="make -j 4 BOPT=O 2>&1 >/dev/null | grep -i"' >> ~/.bashrc
CMD ["bash"]

