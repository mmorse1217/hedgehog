#ifndef __RESULTS_WRITER_HPP__
#define __RESULTS_WRITER_HPP__
#include "common/nummat.hpp"
void write_petsc_vec_to_csv(Vec vec, string filename);
void write_petsc_vecs_to_csv(vector<Vec> vecs, string filename);


#endif
