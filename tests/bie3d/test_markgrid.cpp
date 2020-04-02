#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "bie3d/markgrid.hpp"
#include "bie3d/spatial_grid.hpp"
#include "common/stats.hpp"
#include "bdry3d/on_surface_point.hpp"
#include "common/utils.hpp"
#include "../utils/evaluation_utils.hpp"
#include "common/vtk_writer.hpp"
#include <time.h>
using namespace Ebi;
using Markgrid::SpatialGrid;
using Markgrid::NearFieldMap;

void test_closest_point_single_patch(DblNumMat targets, PatchSurfFaceMap* face_map,
        DblNumMat actual_closest_points){
    int num_targets = targets.n();
    assert(num_targets == actual_closest_points.n());
    NumVec<OnSurfacePoint> closest_points(num_targets);
    Markgrid::mark_near_field(targets, face_map, closest_points);

    Vec closest_points_vec = compute_3d_position_of_on_surface_points(
            closest_points, face_map);
    DblNumMat closest_points_3d_position_local = get_local_vector(DIM, num_targets, closest_points_vec);
    for(int i = 0; i < num_targets; i++){
        Point3 computed_position(closest_points_3d_position_local.clmdata(i));
        Point3 true_position(actual_closest_points.clmdata(i));
        Point3 target(targets.clmdata(i));

        CHECK( fabs((target - computed_position).length() - (target - true_position).length() ) <=1e-13);
        CHECK((computed_position - true_position).length() <= true_position.length()*1e-13 + 1e-13);

    }

    closest_points_3d_position_local.restore_local_vector();
    VecDestroy(&closest_points_vec);

}

void test_closest_point_two_patches(DblNumMat targets, PatchSurfFaceMap* face_map,
        DblNumMat actual_closest_points){
    int num_targets = targets.n();
    assert(num_targets == actual_closest_points.n());
    NearFieldMap closest_points;
    Markgrid::mark_near_field(targets, face_map, closest_points);


    

    for(int i = 0; i < num_targets; i++){
        vector<OnSurfacePoint> near_points = closest_points[i];
        Point3 true_position(actual_closest_points.clmdata(i));
        Point3 target(targets.clmdata(i));
        //assert(near_points.size() == 2);
        
        OnSurfacePoint first_closest_point = near_points[0]; 
        OnSurfacePoint second_closest_point = near_points[1]; 
        NumVec<OnSurfacePoint> temp(1);
        temp(0) = first_closest_point;
        vector<int> pids = {0,1};
        DblNumMat target_temp(3,1,false,target.array());
        dump_vtk_data_for_paraview(target_temp,
                temp, i,
                pids, face_map );

    // Check that the distance from the closest point on each patch is
    // the same
    CHECK( fabs(first_closest_point.distance_from_target  - 
                second_closest_point.distance_from_target) <= 1e-14);

    FaceMapSubPatch* patch = 
        (FaceMapSubPatch*) face_map->patches()[first_closest_point.parent_patch]; 
    Point3 first_closest_point_position;
    patch->xy_to_patch_coords(
            first_closest_point.parametric_coordinates,
            PatchSamples::EVAL_VL,
            (double*) &first_closest_point_position);

    // Check that the distance from the closest point on patch is the
    // same as the actual distance from the true closest point
    CHECK((target-true_position).length()- 
            (target - first_closest_point_position).length() <= 1e-14);

    patch = (FaceMapSubPatch*) face_map->patches()[second_closest_point.parent_patch]; 
    Point3 second_closest_point_position;
    patch->xy_to_patch_coords(
            second_closest_point.parametric_coordinates,
            PatchSamples::EVAL_VL,
            (double*) &second_closest_point_position);

    // Check that the distance from the closest point on patch is the
    // same as the actual distance from the true closest point
    CHECK((target-true_position).length() -
            (target - second_closest_point_position).length() <= 1e-14);

    // Check that the closest point on each patch is actually the same
    // point in R^3
    CHECK((first_closest_point_position - second_closest_point_position).length() <= 1e-14);
    } 
}


void check_normal_projection(DblNumMat targets, PatchSurfFaceMap* face_map){
    NearFieldMap closest_points;
    Markgrid::mark_near_field(targets, face_map, closest_points);
    REQUIRE(closest_points.size() == size_t(targets.n()));
    
    for(auto const& points_near_target: closest_points){
        int target_id = points_near_target.first;
        Point3 target(targets.clmdata(target_id));

        auto closest_on_surface_points = points_near_target.second;
        
        for(OnSurfacePoint on_surface_point : closest_on_surface_points){
            Point2 uv = on_surface_point.parametric_coordinates;
            int patch_containing_point = on_surface_point.parent_patch;

            // compute normal at and position of closest point 
            auto patch = face_map->subpatch(patch_containing_point);
            
            Point3 position;
            patch->xy_to_patch_coords(uv.array(), PatchSamples::EVAL_VL, position.array());
            Point3 normal = patch->normal(uv.array());
            double projection_magnitude = fabs(dot(normal, target - position));
            double x_minus_y_length = (target - position).length();
            CHECK(fabs(projection_magnitude - x_minus_y_length) <= x_minus_y_length*1e-14 + 1e-14);
        }

    }

}

void test_two_closest_points_are_equidistance(DblNumMat targets, 
        PatchSurfFaceMap* face_map){
    NearFieldMap closest_points;
    Markgrid::mark_near_field(targets, face_map, closest_points);
    REQUIRE(closest_points.size() == size_t(targets.n()));

    for(auto const& points_near_target: closest_points){
        int target_id = points_near_target.first;
        Point3 target(targets.clmdata(target_id));

        auto closest_on_surface_points = points_near_target.second;
        //assert(closest_on_surface_points.size() == 2);
        OnSurfacePoint first_closest_point = closest_on_surface_points[0]; 
        OnSurfacePoint second_closest_point = closest_on_surface_points[1]; 

        // Check that the distance from the closest point on each patch is
        // the same
        CHECK( fabs(first_closest_point.distance_from_target  - 
                    second_closest_point.distance_from_target) <= 1e-14);
    }

}

void test_in_out_marking_parallel_copy_of_targets(size_t num_samples, 
        PatchSurfFaceMap* face_map, double r,
        DomainMembership in_or_out){

    // sample the face-map
    PatchSamples* patch_samples = new PatchSamples("", "");
    patch_samples->bdry() = face_map;
    vector<int> patches(face_map->patches().size(), 0);
    patch_samples->patch_partition() = patches;
    patch_samples->setup(); 

    // scale the points to be either inside or outside
    Vec points;
    VecDuplicate(patch_samples->sample_point_3d_position(), &points);
    VecCopy(patch_samples->sample_point_3d_position(), points);
    if(in_or_out == INSIDE){
        VecAXPY(points, -(1.-r), patch_samples->sample_point_normal());
    } else if(in_or_out == OUTSIDE){
        VecAXPY(points, (1.-r), patch_samples->sample_point_normal());

    }
    //VecScale(points, r);
    DblNumMat points_on_sphere = get_local_vector(3, patch_samples->local_num_sample_points(), points);

    auto samples_on_surface_points = patch_samples->sample_point_as_on_surface_point();
    DblNumMat true_points= get_local_vector(3, patch_samples->local_num_sample_points(),
            patch_samples->sample_point_3d_position());

    Vec e;
    VecCreateMPI(MPI_COMM_WORLD, patch_samples->local_num_sample_points(), PETSC_DETERMINE, &e);
    write_general_points_to_vtk(patch_samples->sample_point_3d_position(), 
            1, "true_points.vtp", e, "output/");
    write_general_points_to_vtk(points, 1, "points_to_mark.vtp", e, "output/");
    // find the closest on surface point and verify the points are in/out
    //NumVec<OnSurfacePoint> on_surface_points = 
    //    Markgrid::compute_closest_on_surface_points(points_on_sphere,face_map);
    NumVec<OnSurfacePoint> on_surface_points = 
        Markgrid::mark_target_points(points_on_sphere, face_map,false);
    Vec computed_closest;
    VecDuplicate(points, &computed_closest);
    VecCopy(points, computed_closest);
    DblNumMat computed_pts(3, computed_closest);
    for(int pi = 0; pi < points_on_sphere.n(); pi++){
        auto p = on_surface_points(pi);
        REQUIRE(p.target_index == pi);
        CHECK( p.inside_domain == in_or_out);
        auto patch = face_map->patch(p.parent_patch);
        Point3 ret;
        patch->xy_to_patch_coords(p.parametric_coordinates, PatchSamples::EVAL_VL, ret.array());
        Point3 true_point(true_points.clmdata(pi));
        Point2 true_uv = samples_on_surface_points(pi).parametric_coordinates;
        
        for (int i = 0; i < 3; i++) {
            computed_pts(i,pi) = ret(i);
        }

        CHECK((ret - true_point).l2() <=1e-6);
    }
    write_general_points_to_vtk(computed_closest, 1, "computed_closest.vtp", e, "output/");

    points_on_sphere.restore_local_vector();
    VecDestroy(&points);
    delete patch_samples;
}

// create list of factors to scale the on-surface samples by, in such a fashion
// that in/out is known
vector<double> create_scale_factors(){
    vector<double> scale_factors;
    // scale_factor  = .1, .2, ..., .9
    //int num_interior = 10;
    
    int num_interior = 8;
    for(int i = 1; i < num_interior; i++){
        double scale_factor = double(i)/double(num_interior);
        scale_factors.push_back(scale_factor);
    }
    // scale_factor= .9, .99, .999, .9999
    //int num_near_int = 5;
    int num_near_int = 2;
    for(int i = 1; i < num_near_int; i++){
        double scale_factor = 1. - pow(10, -(i) );
        scale_factors.push_back(scale_factor);
    }
    // scale_factor = 1.000001, 1.00001, 1.0001, 1.001
    //int num_near_ext = 4;
    /*
    int num_near_ext = 1;
    for(int i = num_near_ext; i > 0; i--){
        double scale_factor = 1. + pow(10., -(i) );
        scale_factors.push_back(scale_factor);
    }*/
    // scale_factor = 1.025, 1.05, 1.075, 1.1
    //int num_exterior= 4;
    int num_exterior= 3;
    for(int i = 1; i < num_exterior ; i++){
        double scale_factor = 1. + 1e-1*double(i)/double(num_exterior);
        scale_factors.push_back(scale_factor);
    }
    return scale_factors;

}





TEST_CASE("Test closest point on single patch", 
        "[markgrid][geom][closest-point][critical]"){
    //stats.add_result("closest point total iterations", 0.);
    //stats.add_result("closest point total calls", 0.);
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    SECTION("test closest point to single flat patch"){
        // Single patch defined on [-1,1] x [-1, 1] x 0

        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        // Load patch + polynomial from file
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_meshes/wrl/flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_meshes/wrl/flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_meshes/poly/flat_patch.poly");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();


        assert(face_map->patches().size() == 1); // only a single flat patch
        auto patch = FaceMapSubPatch::as_subpatch(face_map->patches()[0]);
        {
            DblNumMat targets(DIM, 2*2*3);
            DblNumMat closest_points(targets.m(), targets.n());
            // test the closest point on patch to points in the xy plane: 
            // (\pm 2, \pm 2, 0); solution = (\pm 1, \pm 1, 0)
            int index =0;
            for(int i = 0; i < 2; i++){
                for(int j = 0; j < 2; j++){
                    double x_sign = pow(-1., i);
                    double y_sign = pow(-1., j);
                    for(int z = -1; z < 2; z++){
                        Point3 target(targets.clmdata(index));
                        target = Point3(x_sign*2., y_sign*2., 2.*double(z));
                        //Point3 target(x_sign*2., y_sign*2., 2.*double(z));
                        Point3 closest_point(closest_points.clmdata(index));
                        closest_point = Point3(x_sign, y_sign, 0.);
                        //Point3 actual_closest_point(x_sign, y_sign, 0.);
                        index++;
                    }
                }
            }
            test_closest_point_single_patch(targets, face_map, closest_points);
        }
        {
            int num_samples = 10;
            DblNumMat targets(DIM, num_samples*num_samples*3);
            DblNumMat closest_points(targets.m(), targets.n());
            // (-1 + 2i/N,-1 + 2j/N, \pm 1); solution = (-1 + 2i/N,-1 + 2j/N, 0)
            int index =0;
            double h = 2./double(num_samples);
            for(int i = 0; i <= num_samples; i++){
                for(int j = 0; j <= num_samples; j++){
                    for(int sign = -1; sign < 2; sign+=2){
                        Point3 target(targets.clmdata(index));
                        target = Point3(-1. + i*h, -1. + j*h, double(sign)*1e-12);
                        
                        Point3 closest_point(closest_points.clmdata(index));
                        closest_point = Point3(-1. + i*h, -1. + j*h, 0.);

                        index++;

                    }
                }
            }
            test_closest_point_single_patch(targets, face_map, closest_points);
        }
        {
            int num_samples = 10;
            DblNumMat targets(DIM, 2*2);
            DblNumMat closest_points(targets.m(), targets.n());
            // closest point to the corners should be themselves
            int index = 0;
            for(int i = 0; i < 2; i++){
                for(int j = 0; j < 2; j++){

                    double x_sign = pow(-1., i);
                    double y_sign = pow(-1., j);
                    Point3 target(targets.clmdata(index));
                    target = Point3(x_sign, y_sign, 0.);

                    Point3 closest_point(closest_points.clmdata(index));
                    closest_point = Point3(x_sign, y_sign, 0.);
                    index++;


                }
            }
            test_closest_point_single_patch(targets, face_map, closest_points);
        }
        
        {
            // closest point to an interior on-surface point should be itself
            DblNumMat target(3,1);
            DblNumMat actual_closest_point(3,1);
            for(int d = 0; d < DIM; d++){
                target(d,0) = 0.;
                actual_closest_point(d,0) = 0.;
            }
            test_closest_point_single_patch(target, face_map, actual_closest_point);
        }

        delete face_map;
    }
    stats.print_results();
}

TEST_CASE("Test closest points on two patches", "[markgrid][critical]"){
    SECTION("test closest point to two flat patches sharing an edge"){
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        // Load patch + polynomial from file
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_meshes/wrl/two_patches_shared_edge.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_meshes/wrl/two_patches_shared_edge.wrl");
        PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_meshes/poly/two_patches_shared_edge.poly");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        assert(face_map->patches().size() == 2);

        int num_samples = 10;
        double h = .2/(num_samples);
        // Target points along line above the shared edge
        // solution is the projection of the line onto the shared edge
        {
            DblNumMat targets(DIM, num_samples+1);
            DblNumMat closest_points(targets.m(), targets.n());
            int index =0;
            for(int i = 0; i <= num_samples; i++){
                Point3 target(0., .2, -.1 + i*h);;
                Point3 actual_closest_point(0., .1, -.1 + i*h);
                for(int d= 0; d < DIM; d++){
                    targets(d,i) = target(d);
                    closest_points(d,i) = actual_closest_point(d);
                }
            }
            test_closest_point_two_patches(targets, face_map, closest_points);
        }
        {
            DblNumMat targets(DIM, pow(num_samples+1, 2));
            DblNumMat closest_points(targets.m(), targets.n());
            int index =0;
            for(int i = 0; i <= num_samples; i++){
                for(int j = 0; j <= num_samples; j++){
                    Point3 target(0., -.1 + j*h, -.1 + i*h);
                    for(int d= 0; d < DIM; d++){
                        targets(d,index) = target(d);
                    }
                    index++;
                }
            }
            check_normal_projection(targets, face_map);
            test_two_closest_points_are_equidistance(targets, face_map);

        }
    stats.print_results();
    }
    SECTION("test closest point to two flat patches sharing an edge"){

        /**
         * Two flat patches with a shared vertex at (0,1,0). Medial axis is
         * given by the line x=z=0.
         */
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        // Load patch + polynomial from file
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_meshes/wrl/two_patches_shared_vertex.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_meshes/wrl/two_patches_shared_vertex.wrl");
        PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_meshes/poly/two_patches_shared_vertex.poly");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        assert(face_map->patches().size() == 2);
        int num_samples = 10;

        DblNumMat targets(DIM, 2*(num_samples+1));
        double h = .2/num_samples;
        for(int i = 0; i <= 2*num_samples; i++){
            Point3 target(targets.clmdata(i));
            target = Point3(0., -.1 + i*h, 0.);
        }
        check_normal_projection(targets, face_map);
        test_two_closest_points_are_equidistance(targets, face_map);

    stats.print_results();
    }
    SECTION("test closest point to two flat patches intersecting along an edge in interior"){
        /**
         * Two flat patches with a shared edge at along x = 0, y = 0. 
         * Medial axis is the plane x =0;
         */
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        // Load patch + polynomial from file
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_meshes/wrl/two_patches_nonmanifold.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_meshes/wrl/two_patches_nonmanifold.wrl");
        PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_meshes/poly/two_patches_nonmanifold.poly");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        assert(face_map->patches().size() == 2);
        int num_samples = 10;
        double h = .2/double(num_samples);
        {
            // check that points on the intersection return the same physical
            // closest point in space on both patches
            DblNumMat targets(DIM, (num_samples+1));
            DblNumMat closest_points(targets.m(), targets.n());
            int index = 0;

            for(int i = 0; i <= num_samples; i++){
                Point3 target(0., 0., -.1 + i*h );
                Point3 actual_closest_point(0.,0.,-.1 + i*h );
                    for(int d= 0; d < DIM; d++){
                        targets(d,index) = target(d);
                        closest_points(d,index) = actual_closest_point(d);
                    }
                    index++;


            }
            test_closest_point_two_patches(targets, face_map, closest_points);
            check_normal_projection(targets, face_map);
            test_two_closest_points_are_equidistance(targets, face_map);
        }
        {
            // Check samples from the plane x =0, and ensure the closest points
            // meet the optimality condition and the closest point on each patch
            // is the same distance.
            DblNumMat targets(DIM, pow(num_samples+1,2));
            int index = 0;
            for(int i = 0; i <= num_samples; i++){
                for(int j = 0; j <= num_samples; j++){
                    Point3 target(0., -.1 + j*h, -.1 + i*h );
                    for(int d= 0; d < DIM; d++)
                        targets(d,index) = target(d);
                    index++;
                }
            }
            check_normal_projection(targets, face_map);
            test_two_closest_points_are_equidistance(targets, face_map);

        }
    stats.print_results();
    }
    SECTION("Test find_closest_on_surface_point_in_list()") {
         //
         // Check that find_closest_on_surface_point_in_list() find the closest
         // one in a list of points that are all O(\eps) from the true solution
         // for \eps = 1e-2, ... , 1e-11
         //
        int n = 20;
        OnSurfacePoint closest_point;
        closest_point.distance_from_target = 1e-2;
        
        vector<OnSurfacePoint> on_surface_points(n, closest_point);
        srand(time(NULL)); 
        for(double noise_scale = 1e-2; noise_scale >= 1e-11; noise_scale *= 1e-1){
            for(int k = 0; k < n; k+=2){
                for(int i = 0; i < n; i++){
                    double random_double =  ((double) rand() / (RAND_MAX));
                    on_surface_points[i].distance_from_target += noise_scale*random_double;
                }
                on_surface_points[k] = closest_point;

                OnSurfacePoint computed_closest_point = 
                    Markgrid::find_closest_on_surface_point_in_list(on_surface_points);

                CHECK(fabs(computed_closest_point.distance_from_target - closest_point.distance_from_target) <=1e-14);
            }
        }
         //
         // Check that find_closest_on_surface_point_in_list() find the closest
         // one (only in terms of distance)  in a list of points that are all 
         // O(\eps) from the true solution
         // for \eps = 1e-2, ... , 1e-11 and when there are multiple closest points
         // (say 4)
         //
        for(double noise_scale = 1e-2; noise_scale >= 1e-11; noise_scale *= 1e-1){
            for(int k = 0; k < n; k+=2){
                for(int i = 0; i < n; i++){
                    double random_double =  ((double) rand() / (RAND_MAX));
                    on_surface_points[i].distance_from_target += noise_scale*random_double;
                }
                for(int j = 0; j < 4; j++){
                    // multiple correct values in random places
                    int random_offset = rand()/(RAND_MAX)*(n-1);
                    on_surface_points[(k+random_offset)%n] = closest_point;
                }

                OnSurfacePoint computed_closest_point = 
                    Markgrid::find_closest_on_surface_point_in_list(on_surface_points);
                // make sure the distance is correct
                CHECK(fabs(computed_closest_point.distance_from_target - closest_point.distance_from_target) <=1e-14);
            }
        }

        
    }

    SECTION("Test closest point to two patches via grid search"){
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        // Load patch + polynomial from file
        // Two flat patches of .2 x .2 centered at (0, 0, -.05) and (0, 0, .05)
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/two_small_flat_patch.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/two_small_flat_patch.wrl");
        Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_meshes/poly/two_small_flat_patch.poly");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();
        
        SpatialGrid* grid  = new SpatialGrid(face_map);
        int num_samples = 10;
        double h = .2/(num_samples);

        // generate target points along plane in between two patches
        DblNumMat target_points(DIM, (num_samples+1)*(num_samples+1));
        for(int i = 0; i <= num_samples; i++){
            for(int j = 0; j <= num_samples; j++){
                int index = i*(num_samples+1)+j;
                target_points(0,index) = -.1 + i*h;
                target_points(1,index) = -.1 + j*h;
                target_points(2,index) = 0;
            }
        }
        
        Markgrid::NearFieldMap closest_point_map;
        Markgrid::mark_target_points(target_points,face_map, closest_point_map);

        // map should be the same size as the number of target points 
        REQUIRE(closest_point_map.size() == target_points.n());
        for(auto target_and_closest_points : closest_point_map){

            int target_id = target_and_closest_points.first;
            auto closest_points = target_and_closest_points.second ;
            // should have one closest point per patch
            //REQUIRE(closest_points.size() == 2);
            for(auto on_surface_point : closest_points){
                // both closest points should be .05 away from target
               CHECK(fabs(on_surface_point.distance_from_target -.05) <=1e-14*.1 + 1e-14);
            }
        }
        
        delete grid;
        delete face_map;
    }
}
TEST_CASE("test", "[valgrind]"){
    SECTION("Test closest point to two patches uniform refinement"){
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        // Load patch + polynomial from file
        // Two flat patches of .2 x .2 centered at (0, 0, -.05) and (0, 0, .05)
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/two_small_flat_patch.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/two_small_flat_patch.wrl");
        Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_meshes/poly/two_small_flat_patch.poly");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_uniform(3);
        //face_map->refine_test();
        
        //SpatialGrid* grid  = new SpatialGrid(face_map);
        int num_samples = 10;
        double h = .2/(num_samples);

        // generate target points along plane in between two patches
        DblNumMat target_points(DIM, (num_samples+1)*(num_samples+1));
        for(int i = 0; i <= num_samples; i++){
            for(int j = 0; j <= num_samples; j++){
                int index = i*(num_samples+1)+j;
                target_points(0,index) = -.1 + i*h;
                target_points(1,index) = -.1 + j*h;
                target_points(2,index) = 0;
            }
        }
        
     //auto closest_points =
            //Markgrid::mark_target_points(target_points,face_map);

        Markgrid::NearFieldMap closest_point_map;
        Markgrid::mark_target_points(target_points,face_map, closest_point_map);
        // map should be the same size as the number of target points 
        for(auto target_and_closest_points : closest_point_map){

            int target_id = target_and_closest_points.first;
            auto closest_points = target_and_closest_points.second ;
            
            OnSurfacePoint closest_point = 
                Markgrid::find_closest_on_surface_point_in_list(closest_points);

            // should have one closest point per patch
                // both closest points should be .05 away from target
               CHECK(fabs(closest_point.distance_from_target -.05) <=1e-14*.1 + 1e-14);
        }
        /*REQUIRE(closest_points.m() == target_points.n());
        for(int i = 0; i < closest_points.m(); i++){
            auto closest_point =  closest_points(i);

            // should have one closest point per patch
                // both closest points should be .05 away from target
               CHECK(fabs(closest_point.distance_from_target -.05) <=1e-14*.1 + 1e-14);
        }*/
        
        delete face_map;
        cout << DBL_MAX << endl;
    }



}


TEST_CASE("Markgrid slow", "[markgrid][geom][closest-point]"){
    //stats.add_result("closest point total iterations", 0.);
    //stats.add_result("closest point total calls", 0.);


    SECTION("test near inside/outside point marking for a blob"){
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;

        // Load cube blob thing
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_meshes/wrl/cube.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_facemap_refinement_factor", "1");
        PetscOptionsSetValue(NULL, "-bis3d_spacing", ".1");
        PetscOptionsSetValue(NULL, "-bd3d_facemap_patch_order", "8");
        PetscOptionsSetValue(NULL, "-bdsurf_interpolant_spacing", "1.");
        PetscOptionsSetValue(NULL, "-bdsurf_interpolate", "0");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        size_t num_samples = 10;
        //double r_int = .999995;
        double r_int = .99;
        double r_ext = 1.000005;
        //double r_int = .75;
        //double r_ext = 1.2;

        // Generate samples on two  parallel copies of the surface, 
        // one inside and one outside of the boundary and
        // check that they're marked as inside and outside.
        //test_in_out_marking_parallel_copy_of_targets(num_samples, face_map, r_ext, OUTSIDE);
        test_in_out_marking_parallel_copy_of_targets(num_samples, face_map, r_int, INSIDE);
        

    }
    
    SECTION("Test collect_nearby_on_surface_points"){
        
        // Initialize 
        int n = 20;
        OnSurfacePoint closest_point;
        closest_point.distance_from_target = 1e-2;
        
        vector<OnSurfacePoint> on_surface_points(n, closest_point);
        srand(time(NULL)); 
        vector<int> random_offsets(n);
        for(int i = 0; i < n; i++)
            random_offsets[i] = i;
         //
         // Check that collect_nearby_on_surface_points() find the all the closest
         // ones (w.r.t distance_from_target) in a list of points that are all 
         // O(\eps) from the true solution
         // for \eps = 1e-2, ... , 1e-11          
         // 
        for(int num_near_points = 1; num_near_points < n/2; num_near_points++){
            for(double noise_scale = 1e-2; noise_scale >= 1e-11; noise_scale *= 1e-1){
                for(int k = 0; k < n; k+=2){

                    // Generate random noise for other closest points
                    for(int i = 0; i < n; i++){
                        double random_double =  ((double) rand() / (RAND_MAX));
                        on_surface_points[i].distance_from_target += noise_scale*random_double;
                    }
                    
                    // add num_near_points copies of the closest point
                    random_shuffle(random_offsets.begin(), random_offsets.end());
                    for(int j = 0; j < num_near_points; j++){
                        // multiple correct values in random places
                        int random_offset = random_offsets[j];
                        int index = (k+random_offset)%n;
                        on_surface_points[index] = closest_point;
                    }

                    vector<OnSurfacePoint> computed_nearby_points = 
                        Markgrid::collect_nearby_on_surface_points(on_surface_points, closest_point, 1e-12);
                    // make sure we found as many nearby on-surface points as
                    // we generated...
                    int num_computed_nearby_points =computed_nearby_points.size();
                    CHECK(num_computed_nearby_points == num_near_points); 
                    // make sure the distance is correct
                    for(int i =0; i < num_computed_nearby_points; i++){
                        OnSurfacePoint computed_nearby_point = computed_nearby_points[i];
                        CHECK(fabs(computed_nearby_point.distance_from_target - closest_point.distance_from_target) <=1e-14);
                    }
                }
            }
        }

    }
    

/*
    SECTION("test near inside/outside point marking for a pipe"){
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;

        // Load cube blob thing
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_meshes/wrl/pipe.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_meshes/wrl/pipe.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_facemap_refinement_factor", "1");
        PetscOptionsSetValue(NULL, "-bis3d_spacing", "1.");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        size_t num_samples = 10;
        double r_int = .995;
        double r_ext = 1.005;

        // Generate two parallel copies of the surface, 
        // one inside and one outside of the boundary and
        // check that they're marked as inside and outside.
        //test_in_out_marking_parallel_copy_of_targets(num_samples, face_map, r_ext, OUTSIDE);
        //test_in_out_marking_parallel_copy_of_targets(num_samples, face_map, r_int, INSIDE);
    }

    SECTION("Test qbkix point generation"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        Options::set_value_petsc_opts("-near_interpolation_num_samples", "4");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".5");
        //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        Options::set_value_petsc_opts("-bis3d_spacing", "1");
        
        int qbkix_order = int(Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
        double boundary_spacing = Options::get_double_from_petsc_opts("-boundary_distance_ratio");

        // initialize surface
        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();
        int num_patches = face_map->patches().size();
        int num_samples_per_patch = int(ceil(1./Options::get_double_from_petsc_opts("-bis3d_spacing")))+1;
        num_samples_per_patch *= num_samples_per_patch;
        int num_qbkix_points = num_samples_per_patch*num_patches;

        // sample surface
        PatchSamples* patch_samples = new PatchSamples("", "" );
        patch_samples->bdry() = face_map;
        vector<int> patch_partition(face_map->patches().size(), 0);
        patch_samples->patch_partition() = patch_partition;
        patch_samples->setup();
        
        
        vector<int> patch_index(num_patches, 0);
        for(int i = 0; i < num_patches; i++){
            patch_index[i] = i;
        }
        
        Vec qbkix_points = patch_samples->generate_qbkix_points_from_sample_points(vector<int>(1, floor((qbkix_order-1)/2)), patch_index);
        
        DblNumMat qbkix_points_local  = get_local_vector(DIM, num_qbkix_points, qbkix_points);
        DblNumMat sample_points_local = get_local_vector(DIM, num_qbkix_points, patch_samples->sample_point_3d_position());
        for(int pi = 0; pi < num_patches; pi++){
            FaceMapSubPatch* patch = dynamic_cast<FaceMapSubPatch*>(face_map->patches()[pi]);
            for(int qi = 0; qi < num_samples_per_patch; qi++){
                int index = pi*num_samples_per_patch + qi;
                Point3 sample_point(sample_points_local.clmdata(index));
                Point3 qbkix_point(qbkix_points_local.clmdata(index));
                CHECK( fabs((sample_point - qbkix_point).length() - 
                            double(qbkix_order/2)*patch->characteristic_length()*boundary_spacing) <=1e-14);

            }
        }
        qbkix_points_local.restore_local_vector();
        sample_points_local.restore_local_vector();
        VecDestroy(&qbkix_points);


    }*/
    
  
/*
SECTION("Test phase 2 refinement"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-bis3d_spacing", "1");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");
        
        int qbkix_order = int(Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
        double boundary_spacing = Options::get_double_from_petsc_opts("-boundary_distance_ratio");

        // initialize surface
        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();
        //FaceMapSubPatch* patch = dynamic_cast<FaceMapSubPatch*>(face_map->patches()[0]);
        //cout << patch->_parent_id << endl;
        int num_patches = face_map->patches().size();
        int num_samples_per_patch = int(ceil(1./Options::get_double_from_petsc_opts("-bis3d_spacing")))+1;
        num_samples_per_patch *= num_samples_per_patch;
        int num_qbkix_points = num_samples_per_patch*num_patches;

        
        // Refine until points are inside
        
        //p4est_connectivity_t* connectivity = build_connectivity_from_face_map(face_map);
        //p4est_t* p4est = p4est_new(MPI_COMM_WORLD, connectivity, sizeof(RefinementData<FaceMapSubPatch>), NULL, NULL);
p4est_t* p4est = face_map->_p4est;

        //refine_patches_for_qbkix_point_location(p4est, face_map);
        vector<int> stage_two_indices(1, (qbkix_order-1)/2);

        //refine_patches_point_location(p4est, face_map, stage_two_ref_data);
        refine_patches_midpoint_near_medial_axis(p4est, face_map, stage_two_indices);

        vector<Patch*> subpatches = p4est_to_face_map_subpatches(p4est, face_map);
        face_map->patches() = subpatches;

        num_patches = face_map->patches().size();
        // sample surface
        PatchSamples* patch_samples = new PatchSamples("", "" );
        patch_samples->bdry() = face_map;
        vector<int> patch_partition(face_map->patches().size(), 0);
        patch_samples->patch_partition() = patch_partition;
        patch_samples->setup();

        vector<int> patch_ids(face_map->patches().size(),-1);
        for(int pi = 0; pi < face_map->patches().size(); pi++)
            patch_ids[pi] = pi;
       
        vector<int> qbkix_indices(1,(qbkix_order-1)/2);

        Vec qbkix_points = patch_samples->generate_qbkix_points_from_sample_points(
                qbkix_indices, patch_ids);


        DblNumMat qbkix_points_local = get_local_vector(DIM, num_patches*num_samples_per_patch*qbkix_indices.size(), qbkix_points);


        NumVec<OnSurfacePoint> on_surface_points = Markgrid::mark_target_points(qbkix_points_local, face_map);

        for(int qi =0; qi < on_surface_points.m(); qi++){
            OnSurfacePoint o = on_surface_points(qi);

            switch(o.region){
                case FAR:
                    CHECK(o.parent_patch == -1);
                    break;
                case NEAR:
                    CHECK(o.parent_patch == int(qi/num_samples_per_patch));
                    FaceMapSubPatch* patch = dynamic_cast<FaceMapSubPatch*>(face_map->patches()[o.parent_patch]);

                    // TODO CHECK THAT THIS IS A SANE TEST I DONT THINK THIS IS
                    // RELAVENT ANY LONGER
                    double min_distance = patch->characteristic_length()/2.;
                    CHECK(o.distance_from_target > min_distance);
                    break;

            }
            CHECK(o.inside_domain == INSIDE);
        }

        

    }

  SECTION("Test phase 3 refinement"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/pipe.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/pipe.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".5");
        //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        Options::set_value_petsc_opts("-bis3d_spacing", "1");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");
        
        int qbkix_order = int(Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
        double boundary_spacing = Options::get_double_from_petsc_opts("-boundary_distance_ratio");

        // initialize surface
        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        //face_map->refine_test();
        int num_patches = face_map->patches().size();
        int num_samples_per_patch = int(ceil(1./Options::get_double_from_petsc_opts("-bis3d_spacing")))+1;
        num_samples_per_patch *= num_samples_per_patch;
        int num_qbkix_points = num_samples_per_patch*num_patches;

        
        // Refine until points are inside
        
        p4est_connectivity_t* connectivity = build_connectivity_from_face_map(face_map);
        p4est_t* p4est = p4est_new(MPI_COMM_WORLD, connectivity, sizeof(RefinementData), NULL, NULL);


        //refine_patches_for_qbkix_point_location(p4est, face_map);
        refine_patches_for_fixed_qbkix_points(p4est, face_map);
        
        
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(p4est, face_map);
        face_map->patches() = subpatches;

        num_patches = face_map->patches().size();
        // sample surface
        PatchSamples* patch_samples = new PatchSamples("", "" );
        patch_samples->bdry() = face_map;
        vector<int> patch_partition(face_map->patches().size(), 0);
        patch_samples->patch_partition() = patch_partition;
        patch_samples->setup();

        vector<int> patch_ids(face_map->patches().size(),-1);
        for(int pi = 0; pi < face_map->patches().size(); pi++)
            patch_ids[pi] = pi;
       
        vector<int> qbkix_indices(qbkix_order,0);
        for(int qi = 0; qi < qbkix_order; qi++)
            qbkix_indices[qi] = qi;


        Vec qbkix_points = patch_samples->generate_qbkix_points_from_sample_points(
                qbkix_indices, patch_ids);

        DblNumMat qbkix_points_local = get_local_vector(DIM, num_patches*num_samples_per_patch*qbkix_order, qbkix_points);


        NumVec<OnSurfacePoint> on_surface_points = Markgrid::mark_target_points(qbkix_points_local, face_map);

        for(int qi =0; qi < on_surface_points.m(); qi++){
            CHECK(on_surface_points(qi).region == FAR);
        }

        

    }
    SECTION("Test point full marking, in/out and near/far"){

        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".25");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".03125");
        Options::set_value_petsc_opts("-bis3d_np", "12");

        // initialize surface
        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        // sample surface
        PatchSamples* patch_samples = new PatchSamples("", "" );
        patch_samples->bdry() = face_map;
        vector<int> patch_partition(face_map->patches().size(), 0);
        patch_samples->patch_partition() = patch_partition;
        patch_samples->setup();
        
        vector<double> scale_factors;//create_scale_factors();
        scale_factors.push_back(.8);
        scale_factors.push_back(.999);
        scale_factors.push_back(1.0000001);
        scale_factors.push_back(1.1);

        Vec targets = Test::create_scaled_surface_target_points(scale_factors, patch_samples->sample_point_3d_position());
        int num_targets = Petsc::get_vec_size(targets)/DIM;
        DblNumMat targets_local = get_local_vector(DIM, num_targets, targets);

        // mark the points as near/far/ in/out
        NumVec<OnSurfacePoint> on_surface_points = Markgrid::mark_target_points(targets_local, face_map);

        int num_target_levels = scale_factors.size();
        int num_targets_per_level = num_targets/num_target_levels;
        int cnt = 0;
        
        assert(scale_factors.size()  == 4);

        // inside + far check
        for(int i = 0; i < num_targets_per_level; i++){
            OnSurfacePoint on_surface_point = on_surface_points(cnt++);
            CHECK(on_surface_point.inside_domain == INSIDE);
            CHECK(on_surface_point.region== FAR);
        }
        assert(cnt == num_targets_per_level);
        // inside + near check
        for(int i = 0; i < num_targets_per_level; i++){
            OnSurfacePoint on_surface_point = on_surface_points(cnt++);
            CHECK(on_surface_point.inside_domain == INSIDE);
            CHECK(on_surface_point.region== NEAR);
        }
        assert(cnt == 2*num_targets_per_level);

        // outside + near check
        for(int i = 0; i < num_targets_per_level; i++){
            OnSurfacePoint on_surface_point = on_surface_points(cnt++);
            CHECK(on_surface_point.inside_domain == OUTSIDE);
            CHECK(on_surface_point.region== NEAR);
        }

        assert(cnt == 3*num_targets_per_level);
        for(int i = 0; i < num_targets_per_level; i++){
            OnSurfacePoint on_surface_point = on_surface_points(cnt++);
            CHECK(on_surface_point.inside_domain == OUTSIDE);
            CHECK(on_surface_point.region== FAR);
        }
        assert(cnt == 4*num_targets_per_level);

        VecDestroy(&targets); 
    }*/


    stats.print_results();
    stats.clear();
}

