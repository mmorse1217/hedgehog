#include "regression_utils.hpp"
#include "common/utils.hpp"
#include "../catch.hpp"
#include <bdry3d/patch_surf_analytic.hpp>
#include <bdry3d/patch_surf_blended.hpp>
#include <bdry3d/patch_surf_face_map.hpp>
#include <bdry3d/patch_samples.hpp>
using hedgehog::get_local_vector;

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
           CHECK(fabs(computed_mat(i,j) - true_mat(i,j)) <= 25*2.2204e-16);
        }
    }
}

void Regression::compare_to_regression_data(vector<Vec> computed_data,
                                            vector<string> file_names,
                                            string class_name) {
  assert(computed_data.size() == file_names.size());
  const int num_vecs = computed_data.size();

  string prefix = Regression::build_prefix(class_name);

  for (int i = 0; i < num_vecs; i++) {
    string file_name = file_names[i];
    Vec computed_vec_to_store = computed_data[i];
    Regression::test_vec(computed_vec_to_store, prefix + file_name);
  }
}
void Regression::dump_regression_data(vector<Vec> computed_data,
                                      vector<string> file_names,
                                      string class_name) {
  assert(computed_data.size() == file_names.size());
  const int num_vecs = computed_data.size();

  string prefix = Regression::build_prefix(class_name);

  for (int i = 0; i < num_vecs; i++) {
    string file_name = file_names[i];
    Vec computed_vec_to_store = computed_data[i];
    Regression::dump(computed_vec_to_store, prefix + file_name);
  }
}
string Regression::build_prefix(string class_name){
    string prefix = class_name + "/";
    prefix += Test::get_domain() + "/";
    // Recall: 0 = analytic
    //         1 = blended
    //         2 = face-map
    const int bdtype = Options::get_int_from_petsc_opts("-bdtype");

    if(bdtype == 2){
        prefix += "face_map/";
        if(Options::get_double_from_petsc_opts("-bd3d_facemap_adaptive")){
            assert(fabs(Options:: get_double_from_petsc_opts("-bd3d_facemap_fit_accuracy") - 1e-4) <=1e-16);
            prefix += "adaptive/";
        } else {
            prefix += "no_ref/";

        }
    } else if(bdtype == 1){
        prefix += "blended/";
    } else if(bdtype == 0){
        prefix += "analytic/";
    } else { 
        assert(0);//????
    }
    cout << "prefix: " << prefix << endl;
    cout << "bdtype: " << bdtype << endl;
    return prefix;
}
void Regression::test_vec(Vec vec, string file){
    Vec vec_computed;
    VecDuplicate(vec, &vec_computed);
    Regression::load(vec_computed, file);
    Regression::compare(vec_computed, vec);
}

void setup_samples(hedgehog::PatchSurf* surface,  unique_ptr<hedgehog::PatchSamples>& samples){
    vector<int> partition(surface->patches().size(), 0);
    samples = std::move(unique_ptr<hedgehog::PatchSamples>(new hedgehog::PatchSamples("", "")));
    samples->bdry() = surface;
    samples->patch_partition() = partition;
    samples->setup();
}
void Regression::setup_face_map(unique_ptr<hedgehog::PatchSurf>& surface, 
        unique_ptr<hedgehog::PatchSamples>& samples){
    unique_ptr<hedgehog::PatchSurf> face_map(new hedgehog::PatchSurfFaceMap("BD3D_", "bd3d_"));
    dynamic_cast<hedgehog::PatchSurfFaceMap*>(face_map.get())->_surface_type = 
        Options::get_string_from_petsc_opts("-bd3d_filename") == "wrl_meshes/wrl/newtorus.wrl" ?
        hedgehog::PatchSurfFaceMap::POLYNOMIAL :
        hedgehog::PatchSurfFaceMap::BLENDED;

    dynamic_cast<hedgehog::PatchSurfFaceMap*>(face_map.get())->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    dynamic_cast<hedgehog::PatchSurfFaceMap*>(face_map.get())->refine_test();
    surface = std::move(face_map);
    
    setup_samples(surface.get(), samples);
}
void Regression::setup_blended(unique_ptr<hedgehog::PatchSurf>& surface, 
        unique_ptr<hedgehog::PatchSamples>& samples){
    
    unique_ptr<hedgehog::PatchSurf> blended(new hedgehog::PatchSurfBlended("BD3D_", "bd3d_"));
    blended->setFromOptions();
    blended->setup();
    surface = std::move(blended);
    
    setup_samples(surface.get(), samples);
}
void Regression::setup_analytic(unique_ptr<hedgehog::PatchSurf>& surface, 
        unique_ptr<hedgehog::PatchSamples>& samples){

  unique_ptr<hedgehog::PatchSurf> analytic( new hedgehog::PatchSurfAnalytic("BD3D_", "bd3d_"));
  analytic->setFromOptions();
  analytic->setup();
  surface = std::move(analytic);

  setup_samples(surface.get(), samples);
}
