all:ebilibs

EBI_SRC =	solver_interface.cpp evaluator_interface.cpp evaluator_far.cpp evaluator_on_surface.cpp evaluator_near.cpp evaluator_near_extrapolate.cpp solver_utils.cpp solver_gmres_double_layer.cpp markgrid.cpp collocation_patch_samples.cpp evaluator_near_extrapolate_average.cpp evaluator_near_interpolate.cpp spatial_grid.cpp

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
