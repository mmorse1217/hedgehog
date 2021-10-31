#include "numvec.hpp"
BEGIN_EBI_NAMESPACE
template<>
NumVec<double>::~NumVec<double>() { 
    if(_owndata) { 
        if(_m>0) { delete[] _data; _data = NULL; } 
    } else if(_v){
        PetscScalar* tmp = const_cast<PetscScalar*>(_data);
        VecRestoreArray(_v, &tmp );
    }
}
NumVec<double> get_local_vector(int m, Vec v){
    double* v_ptr;
    VecGetArray(v, &v_ptr);
    DblNumVec local_subvector(m, false, v_ptr);
    local_subvector._v = v;
    return local_subvector;
}

template<>
void NumVec<double>::restore_local_vector(){
    VecRestoreArray(_v, &_data);
}

END_EBI_NAMESPACE
