-include $(PVFMM_DIR)/MakeVariables
#VENDOR   = intel

FFLAGS   = -lstfmm3d -lgfortran -mkl -static-intel -openmp -lmath
INCDIR   = -I$(PVFMM_INC) -I. -I${FFTW_DIR}/include

CXX      = ${CXX_PVFMM}
CXXFLAGS = ${CXXFLAGS_PVFMM} ${INCDIR} -O3 -mavx
LDLIBS  := -L${PVFMM_LIB} -L${FFTW_DIR}/lib  -L${FFTW_DIR}  ${LDFLAGS_PVFMM}

OBJS     = fmm_interface.o
MODS     = $(OBJS:.o=.mod)
COBJS    = pvfmm_interface.o
BIN      = fmm_interface_test.x

ifeq (${VENDOR},intel)
  FC       = mpif90
  FFLAGS   = -fopenmp 
  LDFLAGS  += -cxxlib
else 
  FC       = gfortran
  FFLAGS   = -fopenmp -Wall -Wcast-align
endif
LDLIBS   +=  -lstdc++ 

#all: ${BIN}


pvfmm_bis_test: pvfmm_bis_test.o pvfmm_bis_interface.o pvfmm_interface.o
	$(CXX)  $(CXXFLAGS) $^ $(LDLIBS) -o $@ 
		
pvfmm_test:  pvfmm_interface.o pvfmm_test.o
	$(CXX)  $(CXXFLAGS) $^ $(LDLIBS) -o $@ 
pvfmm_parallel_test:  pvfmm_interface.o pvfmm_parallel_test.o
	$(CXX)  $(CXXFLAGS) $^ $(LDLIBS) -o $@ 

dense_eval_test:  pvfmm_interface.o dense_eval_test.o
	$(CXX)  $(CXXFLAGS) $^ $(LDLIBS) -o $@ 
pvfmm_tree_test:  pvfmm_interface.o pvfmm_tree_test.o
	$(CXX)  $(CXXFLAGS) $^ $(LDLIBS) -o $@ 

lib: pvfmm_interface.o
	ar rvs libpvfmmwrap.a $^ 
pvfmm_interface.o: pvfmm_interface.h

#pvfmm_test: pvfmm_test.cpp pvfmm_interface.cc
pvfmm_parallel_test.o: pvfmm_parallel_test.cpp pvfmm_interface.cc

#libpvfmmint.a: pvfmm_interface.o
#ar rvs $@ $^
#$(CXX) $(CXXFLAGS) $^ -shared -o  $@

#pvfmm_interface.o: pvfmm_interface.cc pvfmm_interface.h
#$(CXX) $(CXXFLAGS) -Wall -fpic -c $^ -o $@

%.x: %.o ${OBJS} ${COBJS}
	${FC} ${LDFLAGS} $? ${LDLIBS} -o $@

%.o: %.f90 ${MODS}
	${FC} -c ${FFLAGS} ${CPPFLAGS} $<

%.mod: %.f90
	${FC} -c ${FFLAGS} ${CPPFLAGS} $<

clean:
	${RM} -r *.o *.mod ${BIN}

.SECONDARY:
