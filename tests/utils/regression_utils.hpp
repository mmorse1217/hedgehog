#pragma once
#include "common/nummat.hpp"
#include <bdry3d/patch_samples.hpp>
namespace Regression {
using Ebi::DblNumMat;
    void dump(Vec v, string data_name);
    void dump(DblNumMat m, string data_name);

    void load(Vec& v, string data_name);
    void load(DblNumMat& m, string data_name);

    void compare(Vec computed_vec, Vec true_vec);
    void compare(DblNumMat computed_mat, DblNumMat true_mat);
    string build_prefix(string class_name);
    void dump_regression_data(vector<Vec> computed_data, vector<string> file_names, string class_name);
    void compare_to_regression_data(vector<Vec> computed_data, vector<string> file_names, string class_name);
    void test_vec(Vec vec, string file);
    void setup_face_map(unique_ptr<Ebi::PatchSurf> &surface,
                        unique_ptr<Ebi::PatchSamples> &samples);
    void setup_blended(unique_ptr<Ebi::PatchSurf> &surface,
                       unique_ptr<Ebi::PatchSamples> &samples);
    void setup_analytic(unique_ptr<Ebi::PatchSurf> &surface,
                       unique_ptr<Ebi::PatchSamples> &samples);
};
