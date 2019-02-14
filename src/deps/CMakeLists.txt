set(Hedgehog_dep_src
    petsc/vector.cpp
)
set(Hedgehog_dep_header
    petsc/petsc.hpp
    petsc/vector.hpp
)
PREPEND(Hedgehog_dep_src)
PREPEND(Hedgehog_dep_header)

set(Hedgehog_dep ${Hedgehog_dep_src} 
    ${Hedgehog_dep_header} PARENT_SCOPE)
