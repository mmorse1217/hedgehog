#ifndef __WRITE_LATEX_HPP__
#define __WRITE_LATEX_HPP__
#include "nummat.hpp"
namespace hedgehog {
    string mat_to_latex_table(DblNumMat data, vector<string> column_names);
};
#endif
