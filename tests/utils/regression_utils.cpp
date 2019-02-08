#include "regression_utils.hpp"
#include "common/utils.hpp"
#include "../catch.hpp"

using Ebi::get_local_vector;

void Regression::dump(Vec v, string data_name){
    DblNumMat m = get_local_vector(1, Petsc::get_vec_size(v), v);
    dump(m, data_name);
    m.restore_local_vector();
}

void Regression::load(Vec& v, string data_name){
    DblNumMat m = get_local_vector(1, Petsc::get_vec_size(v), v);
    load(m, data_name);
    m.restore_local_vector();
}

void Regression::compare(Vec computed_vec, Vec true_vec){
    DblNumMat computed_mat = get_local_vector(1, Petsc::get_vec_size(computed_vec), computed_vec);
    DblNumMat true_mat = get_local_vector(1, Petsc::get_vec_size(true_vec), true_vec);
    
    // should be the same size
    assert(computed_mat.n() == true_mat.n());

    compare(computed_mat, true_mat);

    computed_mat.restore_local_vector();
    true_mat.restore_local_vector();
}

void Regression::dump(DblNumMat m, string data_name){
    cout.precision(16);
    // write to file tests/regression_data/data_name.reg
    Debug::save_mat(m, "tests/regression_data/"+data_name);
}

void Regression::load(DblNumMat& m, string data_name){
    cout.precision(16);
    // read from file tests/regression_data/data_name.reg
    Debug::load_mat(m, "tests/regression_data/"+data_name);
}

void Regression::compare(DblNumMat computed_mat, DblNumMat true_mat){
    assert(computed_mat.m() == true_mat.m());
    assert(computed_mat.n() == true_mat.n());
    for (int i = 0; i < computed_mat.m(); i++) {
        for (int j = 0; j < computed_mat.n(); j++) {
           CHECK(fabs(computed_mat(i,j) - true_mat(i,j)) <= 1e-16);
        }
    }
}
