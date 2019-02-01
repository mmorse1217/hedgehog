all:ebilibs

EBI_SRC =	analytic_surf.cpp patch_surf.cpp patch_surf_analytic.cpp patch_surf_blended.cpp patch_surf_face_map.cpp patch_samples.cpp face_map_subpatch.cpp p4est_interface.cpp

EBI_HPP =	

EBI_OBJ =	${EBI_SRC:.cpp=.o}

EBI_TST =

LIBBASE =	libebi

include ${MOBO_DIR}/src/ebi/machine_makefiles/makefile.${MACHINE_NAME}.in

depend: ${EBI_SRC} ${EBI_HPP} ${EBI_TST}
	-${RM} makedepend.${PETSC_ARCH} 
	-${DEPCC} ${DEPORTS} ${PETSC_INCLUDE} ${EBI_INCLUDE} ${EBI_SRC} ${EBI_TST} > makedepend.${PETSC_ARCH}


include makedepend.${PETSC_ARCH}

# DO NOT DELETE
