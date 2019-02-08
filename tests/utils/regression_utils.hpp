#ifndef __REGRESSION_UTILS_HPP__
#define __REGRESSION_UTILS_HPP__
#include "common/nummat.hpp"
namespace Regression {
using Ebi::DblNumMat;
    void dump(Vec v, string data_name);
    void dump(DblNumMat m, string data_name);

    void load(Vec& v, string data_name);
    void load(DblNumMat& m, string data_name);

    void compare(Vec computed_vec, Vec true_vec);
    void compare(DblNumMat computed_mat, DblNumMat true_mat);
};
#endif
