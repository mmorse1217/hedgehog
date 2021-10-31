#include "nummat.hpp"
BEGIN_EBI_NAMESPACE
template<>
NumMat<double>::~NumMat<double>() { 
    if(_owndata) { 
        if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
    } else if(_v && __teardown){
        PetscScalar* tmp = const_cast<PetscScalar*>(_data);
        VecRestoreArray(_v, &tmp );
    }
}
NumMat<double> get_local_vector(int m, int n, Vec v){
    double* v_ptr;
    VecGetArray(v, &v_ptr);
    DblNumMat local_subvector(m, n, false, v_ptr);
    local_subvector._v = v;
    return local_subvector;
}

template<>
void NumMat<double>::restore_local_vector(){
    VecRestoreArray(_v, &_data);
}
END_EBI_NAMESPACE
