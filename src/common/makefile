all:ebilibs

EBI_CPPSRC = svdrep.cpp	vecmatop.cpp ebi.cpp ebi_uti.cpp memAlloc.cpp Dcuhre.cpp CmptTbls.cpp kernel3d.cpp exsol3d.cpp nummat.cpp  interpolate.cpp stats.cpp utils.cpp
#syms_alt.cpp

#EBI_CSRC = ffunc.c dsgrid_modified.c
#LIB_FSRC = f_dcuhre.f #qshep3d.f

#c_csauser.c dsgrid3d.c suprld.c splded.c splcwd_new.c cfaerr.c csa3xd.c ffunc.c

EBI_SRC = ${EBI_CPPSRC} ${EBI_CSRC}

EBI_HPP =

EBI_OBJ=${EBI_SRC:.cpp=.o} $(LIB_FSRC:.f=.o) ${EBI_CSRC:.c=.o}

#EBI_TST=tt.cpp ss.cpp
EBI_TST=

LIBBASE=libebi

include ${MOBO_DIR}/src/ebi/machine_makefiles/makefile.${MACHINE_NAME}.in

depend: ${EBI_SRC} ${LIB_FSRC} ${EBI_CSRC} ${EBI_HPP} ${EBI_TST}
	-${RM} makedepend.${PETSC_ARCH}
	-${DEPCC} ${DEPORTS} ${PETSC_INCLUDE} ${EBI_INCLUDE} ${EBI_SRC} ${EBI_TST} > makedepend.${PETSC_ARCH}


include makedepend.${PETSC_ARCH}

# DO NOT DELETE

