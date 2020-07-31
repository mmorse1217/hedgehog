#include "../catch.hpp"
#include "common/nummat.hpp"
#include "bie3d/solver_utils.hpp"
#include "bie3d/markgrid.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "bdry3d/patch_samples.hpp"

using namespace hedgehog;

TEST_CASE("Compare face-map vs. blendsurf","[blendsurf][face-map][geom][no-fmm]"){
    /*
    // 0: Initialize surfaces and patch samplings 
    PatchSurfBlended* blended_surface = new PatchSurfBlended("BD3D_", "bd3d_");
    //PatchSurfAnalytic* blended_surface = new PatchSurfAnalytic("BD3D_", "bd3d_");
    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;

    PatchSamples* blended_patch_samples = new PatchSamples("","");
    PatchSamples* face_map_patch_samples = new PatchSamples("","");

    // Setup 'em up
    //PetscOptionsSetValue(NULL, "-bd3d_filename", "../wrl_files/sphere.wrl");
    PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/cube.wrl");
    PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/cube.wrl");
    blended_surface->setFromOptions();
    blended_surface->setup();
    face_map->setFromOptions();
    face_map->setup();

    vector<int> patch_partition_face_map(face_map->patches().size(), 0 );
    vector<int> patch_partition_blended(blended_surface->patches().size(), 0);

    blended_patch_samples->bdry() = blended_surface;
    face_map_patch_samples->bdry() = face_map;

    face_map_patch_samples->patch_partition() = patch_partition_face_map;
    blended_patch_samples->patch_partition() = patch_partition_blended;

    blended_patch_samples->setup();
    face_map_patch_samples->setup(true);

    //-1: Unit test code written to do the test below
    SECTION("test closest_sample_points()"){
        // Test points on axes
        int num_points = 6;
        DblNumMat test_points(DIM, num_points);
        clear(test_points);

        for(int i = 0; i < 3; i++){
            test_points(i, i) = 1.;
            test_points(i, i+3) = -1.;
        }

        DblNumMat blended_samples_local = get_local_vector(DIM, 
                num_local_points(blended_patch_samples->sample_point_3d_position()),
                blended_patch_samples->sample_point_3d_position());
        DblNumMat face_map_samples_local = get_local_vector(DIM, 
                num_local_points(face_map_patch_samples->sample_point_3d_position()),
                face_map_patch_samples->sample_point_3d_position());

        vector<int> true_closest_points_blended(num_points, -1);
        vector<int> true_closest_points_face_map(num_points, -1);
        
        for(int i = 0; i < num_points; i++){
            Point3 test_point(test_points.clmdata(i));
            double distance = DBL_MAX;

            // Compute the closest blended sample point
            for(int j = 0; j < blended_samples_local.n(); j++){
                Point3 current_sample(blended_samples_local.clmdata(j));
                double distance_to_sample = (test_point - current_sample).length();
                true_closest_points_blended[i]  = distance_to_sample < distance ?
                                j : 
                                true_closest_points_blended[i];
                distance  = distance_to_sample < distance ?
                                distance_to_sample : 
                                distance;
            }
            
            // Compute the closest face map sample point
            distance = DBL_MAX;
            for(int j = 0; j < face_map_samples_local.n(); j++){
                Point3 current_sample(face_map_samples_local.clmdata(j));
                double distance_to_sample = (test_point - current_sample).length();
                true_closest_points_face_map[i]  = distance_to_sample < distance ?
                                j : 
                                true_closest_points_face_map[i];
                distance  = distance_to_sample < distance ?
                                distance_to_sample : 
                                distance;
            }
        }

        // Compute via Markgrid::closest_sample_points()
        vector<int> closest_blended_point_indices = 
            Markgrid::closest_sample_points(blended_patch_samples, test_points);
        vector<int> closest_face_map_point_indices = 
            Markgrid::closest_sample_points(face_map_patch_samples, test_points);
   
        // test
        CHECK( closest_blended_point_indices.size() == num_points);
        CHECK( closest_face_map_point_indices.size() == num_points);
        
        for(int i = 0; i < num_points; i++){
            CHECK(true_closest_points_blended[i] == closest_blended_point_indices[i]);
            CHECK(true_closest_points_face_map[i] == closest_face_map_point_indices[i]);
        }
        // be a good citizen
        face_map_samples_local.restore_local_vector();
        blended_samples_local.restore_local_vector();
    }
    //cout.precision(16);
    SECTION("test compute_closest_on_surface_point()"){
        // Test points on axes
        int num_points = 6;
        DblNumMat test_points(DIM, num_points);
        clear(test_points);

        for(int i = 0; i < 3; i++){
            test_points(i, i) = 1.;
            test_points(i, i+3) = -1.;
        }

        DblNumVec distance_to_closest_point(num_points);
        IntNumVec patch_containing_closest_point(num_points);
        NumVec<Point2> closest_point_parametric_value(num_points);
        
        Markgrid::compute_closest_on_surface_point(
                test_points,
                blended_patch_samples,
                distance_to_closest_point,
                patch_containing_closest_point,
                closest_point_parametric_value);
        for(int i = 0; i < num_points; i++){
            cout << "distance: " << distance_to_closest_point(i) << ", " 
                << "patch: " << patch_containing_closest_point(i) << ", " 
                << "param value: " << closest_point_parametric_value(i) << endl;

        }

    }
    SECTION("Geometry comparison"){
        // 1: Choose a few random sample points in parametric domain of 
        //    each face-map patch (i.e. random samples are evaluated on each patch)
        
        // test points per patch to check
        int num_point_per_dim = 4;
        int num_points_per_patch = num_point_per_dim*num_point_per_dim;
        //int num_points_per_patch = 1;
        NumVec<Point2> test_points_parametric(num_points_per_patch);
        
        vector<Patch*> patches = face_map_patch_samples->bdry()->patches();
        int num_patches = patches.size();

        // TODO generalize for num_points_per_patch != 1; just testing some
        // point for now
        test_points_parametric(0) = Point2(.5, .5); 
        for(int i = 0; i < num_point_per_dim; i++){
            for(int j = 0; j < num_point_per_dim; j++){
                double x = float(i)/num_point_per_dim + .1;
                double y = float(j)/num_point_per_dim+ .1;
                test_points_parametric(i*num_point_per_dim + j) = Point2(x,y);
            }
        }
        // 2: Compute the approximated on-surface geometric information for each 
        //    sample, (i.e. position, partial derivatives, etc.)
        
        NumVec<Point3> face_map_positions(num_points_per_patch*num_patches);
        DblNumMat face_map_position_reformatted(DIM, num_points_per_patch*num_patches);

        NumVec<Point3> face_map_du(num_points_per_patch*num_patches);
        NumVec<Point3> face_map_dv(num_points_per_patch*num_patches);
        
        // Evaluate geometry data at all test_points_parametric values via
        // face_map and store it
        for(int pi = 0; pi < num_patches; pi++){
            FaceMapPatch* current_patch = (FaceMapPatch*) patches[pi];

            for(int i = 0; i < num_points_per_patch; i++){
                Point3 position_and_derivatives[3];
                Point2 current_point_parametric_value = test_points_parametric(i);
                
                // Evaluate
                current_patch->xy_to_patch_coords(
                        current_point_parametric_value,
                        PatchSamples::EVAL_VL|PatchSamples::EVAL_FD,
                        (double*) position_and_derivatives);
                // Store
                face_map_positions(pi*num_points_per_patch + i) = position_and_derivatives[0];
                face_map_du(pi*num_points_per_patch +i) = position_and_derivatives[1];
                face_map_dv(pi*num_points_per_patch +i) = position_and_derivatives[2];
                
                for(int d = 0; d < DIM; d++){
                    face_map_position_reformatted(d, pi*num_points_per_patch+ i) = position_and_derivatives[0](d);
                }
            }
        }
        // 3: Compute the nearest sample point position on the blended surface
        // 4: Compute the nearest on-surface point position (and parametric value)
        //    on the blended surface by factoring out the closest point code from markgrid
        DblNumVec distance_to_closest_point(num_points_per_patch*num_patches);
        IntNumVec patch_containing_closest_point(num_points_per_patch*num_patches);
        NumVec<Point2> closest_point_parametric_value(num_points_per_patch*num_patches);
        
        Markgrid::compute_closest_on_surface_point(
                face_map_position_reformatted,
                blended_patch_samples,
                distance_to_closest_point,
                patch_containing_closest_point,
                closest_point_parametric_value);


        // 5: Evaluate the geomtric data at the closest point on the blended surface
        vector<Patch*> blended_patches = blended_patch_samples->bdry()->patches();
        
        NumVec<Point3> blended_positions(num_points_per_patch*num_patches);
        NumVec<Point3> blended_du(num_points_per_patch*num_patches);
        NumVec<Point3> blended_dv(num_points_per_patch*num_patches);
        for(int i = 0; i < num_patches*num_points_per_patch; i++){
            int pi = patch_containing_closest_point(i);
            BlendedPatch* pith_patch = (BlendedPatch*) blended_patches[pi];
            //AnalyticPatch* pith_patch = (AnalyticPatch*) blended_patches[pi];
            Point2 parameteric_value = closest_point_parametric_value(i);

            //eval geometry data
            Point3 blended_position_and_derivatives[3];
            pith_patch->xy_to_patch_coords(parameteric_value,
                    PatchSamples::EVAL_VL|PatchSamples::EVAL_FD,
                    (double*) blended_position_and_derivatives);
            
            blended_positions(i) = blended_position_and_derivatives[0];
            blended_du(i) = blended_position_and_derivatives[1];
            blended_dv(i) = blended_position_and_derivatives[2];

        }
            vector<int> pos_digit_count(16, 0);
            vector<int> nor_digit_count(16, 0);
        
        // 6: Compare blended and face-map geometric data
        for(int i = 0; i < num_patches*num_points_per_patch; i++){
            Point3 current_blended_pos = blended_positions(i);
            Point3 current_blended_du = blended_du(i);
            Point3 current_blended_dv = blended_dv(i);
            Point3 current_face_map_pos = face_map_positions(i);
            Point3 current_face_map_du = face_map_du(i);
            Point3 current_face_map_dv = face_map_dv(i);
            
            Point3 blended_normal = cross(current_blended_du, current_blended_dv);
            double blended_normal_mag = blended_normal.l2();
            blended_normal /= blended_normal_mag;
            Point3 face_map_normal = cross(current_face_map_du, current_face_map_dv);
            double face_map_normal_mag = face_map_normal.l2();
            face_map_normal /= face_map_normal_mag;
            
            double position_rel_error =  
                (current_blended_pos - current_face_map_pos).l2()/current_blended_pos.l2();
            double normal_rel_error = 
                (blended_normal - face_map_normal).l2()/blended_normal.l2();
            //cout << log10(position_rel_error) << endl;
            //cout << log10(normal_rel_error) << endl;
            int pos_index = int(-log10(position_rel_error));
            int nor_index = int(-log10(normal_rel_error));
            pos_digit_count[pos_index] += 1;
            nor_digit_count[nor_index] += 1;
        }
            cout << "position accuracy:" << endl;
            for(int j = 0; j < 16; j++){
                cout << j << " digits:" << pos_digit_count[j] << endl; 
            }
            cout << "normal accuracy:" << endl;
            for(int j = 0; j < 16; j++){
                cout << j << " digits:" << nor_digit_count[j] << endl; 
            }
    }*/
}
