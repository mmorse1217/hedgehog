#include "../catch.hpp"
#include "common/stats.hpp"
#include "bie3d/markgrid.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/on_surface_point.hpp"
#include <omp.h>
using namespace hedgehog;
using Markgrid::NearFieldMap;
using Markgrid::compute_closest_points_on_patches;


TEST_CASE("Benchmarking ", "[bench]"){
    SECTION("Benchmark closest point optimization"){
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        // Load patch + polynomial from file
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/small_flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/small_flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_files/poly/small_flat_patch.poly");
        
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();
        
        int num_samples_per_dim = 5;
        int num_samples = num_samples_per_dim*num_samples_per_dim; 
        double step = 1./(num_samples_per_dim-1);
        DblNumMat targets(DIM, num_samples);
       /* 
                targets(0, 0) = 0;
                targets(1, 0) = 0;
               targets(2, 0) = 0;*/ 
        for(int i = 0; i < num_samples_per_dim; i++){
            for(int j = 0; j < num_samples_per_dim; j++){
                int index = i*num_samples_per_dim + j;
                // points in [-.2, .2]
                targets(0, index) = .4*i*step - .2;
                targets(1, index) = .4*j*step - .2 ;
                targets(2, index) = 0;
            }
        }

        int num_iter = 1e5;

        double eval_time = omp_get_wtime();

        vector<uint> patches_to_check(1,0); // only one patch
        for(int it = 0; it < num_iter; it++){
            for(int i = 0; i < num_samples_per_dim; i++){
                for(int j = 0; j < num_samples_per_dim; j++){
                    int index = i*num_samples_per_dim + j;
                    Point3 target_point(targets.clmdata(index));

                    NearFieldMap near_patches_to_point = 
                        compute_closest_points_on_patches(target_point, 
                                index, 
                                face_map,
                                patches_to_check);
                }
            }
        }
        eval_time = omp_get_wtime() - eval_time;
        cout.precision(8);
        cout << "number of iterations: " << num_iter*num_samples<< endl;
        cout << "total evaluation time: " << eval_time << endl;
        cout << "mean evaluation time: " << eval_time/double(num_iter*num_samples) << endl;

    }
}
