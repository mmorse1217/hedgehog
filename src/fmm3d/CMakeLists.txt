set(Hedgehog_fmm_src
    fmm_interface.cpp
    pvfmm_bis_interface.cpp
)
set(Hedgehog_fmm_header
    fmm_interface.hpp
    pvfmm_bis_interface.hpp
    pvfmm/pvfmm_laplace_kernel.hpp
    pvfmm/pvfmm_stks_kernel.hpp
    pvfmm/pvfmm_navier_kernel.hpp
    pvfmm/pvfmm_mod_helmholtz_kernel.hpp
)
PREPEND(Hedgehog_fmm_src)
PREPEND(Hedgehog_fmm_header)

set(Hedgehog_fmm ${Hedgehog_fmm_src} 
    ${Hedgehog_fmm_header} PARENT_SCOPE)
