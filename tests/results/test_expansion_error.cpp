#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "common/utils.hpp"
#include "deps/petsc/vector.hpp"
#include "../utils/evaluation_utils.hpp"
#include "common/nummat.hpp"
#include "deps/petsc/vector.hpp"
#include "bie3d/evaluator_qbkix.hpp"
#include <tuple>


double extrapolation_eval_point_alt(double distance_to_closest_sample, double h){
    return -h;
}
using namespace hedgehog;
using namespace Petsc;
//using Petsc::Vector;
typedef tuple<double, double, double, double> Result;
TEST_CASE("Test qbkix evaluation error", "[results][qbkix-error]"){
    
    SECTION(""){
        Options::set_value_petsc_opts("-near_interpolation_num_samples", "12");
        int qbkix_order = Options::get_double_from_petsc_opts("-near_interpolation_num_samples");

        Petsc::Vector target(1, 3, Petsc::Vector::VectorType::mpi, MPI_COMM_WORLD);
        Petsc::Vector closest_point_mock(1, 3, Petsc::Vector::VectorType::mpi, MPI_COMM_WORLD);
        Vec expansion_distance_from_target;
        Vec expansion_point_spacing;

        Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, expansion_point_spacing);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, expansion_distance_from_target);


        Vec singularity_position;
        Vec singularity_strength;
        Petsc::create_mpi_vec(MPI_COMM_WORLD, 1*3, singularity_position);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, singularity_strength);
        DblNumMat singularity_position_local= get_local_vector(1,3,singularity_position);
        singularity_position_local(0,0) = 1.;
        singularity_position_local.restore_local_vector();
        VecSet(singularity_strength, 4*M_PI);
        Vec true_potential= 
            Test::compute_dirichlet_boundary_data(
                    MPI_COMM_WORLD,
                    Kernel3d(111, vector<double>(2,1.)),
                    singularity_position,
                    singularity_strength,
                    target.v_,
                    NULL);
        Vec interp_dir;
        Petsc::create_mpi_vec(MPI_COMM_WORLD, 3, interp_dir);

        DblNumMat interp_dir_local= get_local_vector(1,3,interp_dir);
        interp_dir_local(0,0) = -1.;
        interp_dir_local.restore_local_vector();

        vector<int> qbkix_points_to_generate;
        for(int i =0; i < qbkix_order; i++){
            qbkix_points_to_generate.push_back(i);
        }

        int num_iter = 500;
        double min_error = DBL_MAX;

        vector<Result> error_results;
        error_results.reserve(num_iter*num_iter*num_iter);

        for(int spacing_i = 0; spacing_i < num_iter; spacing_i++){
            for(int dist_i = 0; dist_i < num_iter; dist_i++){
                double singularity_dist = 1.;
                double R_over_rho = double(dist_i)/num_iter*.25;
                double r_over_R = double(spacing_i)/num_iter*2.;
                double distance_from_target = R_over_rho*singularity_dist;
                double expansion_spacing = r_over_R*distance_from_target;

                VecSet(expansion_point_spacing, expansion_spacing);
                VecSet(expansion_distance_from_target, distance_from_target);
                // target point at the origin
                target.set_all(0.);

                // make qbkix points
                Vec qbkix_points =  generate_interior_qbkix_points(target.v_, 
                        closest_point_mock.v_,qbkix_points_to_generate,interp_dir,
                        expansion_distance_from_target,expansion_point_spacing);

                // compute function values at interpolation nodes
                Vec interp_node_values  = 
                    Test::compute_dirichlet_boundary_data(
                            MPI_COMM_WORLD,
                            Kernel3d(111, vector<double>(2,1.)),
                            singularity_position,
                            singularity_strength,
                            qbkix_points,
                            NULL);
                int num_qbkix_points = Petsc::get_vec_size(qbkix_points)/3;
                DblNumMat interp_node_values_local = get_local_vector(1,num_qbkix_points, interp_node_values);
                vector<double> interpolation_nodes = get_interpolation_nodes(EXTRAPOLATE_ONE_SIDE);
                assert(int(interpolation_nodes.size()) == num_qbkix_points);
                DblNumMat interpolated_value(1,1);
                target.get_local(3,1);

            DblNumMat expansion_distance_from_target(1,1);
            expansion_distance_from_target(0,0) = distance_from_target;
            DblNumMat node_spacing(1,1);
            node_spacing(0,0) = expansion_spacing;
                cout << "R:" <<  distance_from_target << ", r: " << expansion_spacing << endl;
                lagrange_extrapolation_bary(target.local(),
                        target.local(),
                        interpolation_nodes,
                        interp_node_values_local,
                        node_spacing,
                        expansion_distance_from_target,
                        1,
                        &extrapolation_eval_point_qbkix,
                        interpolated_value);

                
               
                target.restore_local();
                DblNumMat true_potential_local = get_local_vector(1,1,true_potential);
                double abs_error = fabs( interpolated_value(0,0) - true_potential_local(0,0));
              
                    cout << "computed: " << interpolated_value(0,0) << ", actual: " <<  true_potential_local(0,0) << endl;
                if(abs_error > 1e12 && distance_from_target > .5){

                    min_error = abs_error;
                    cout << endl << endl;
                    cout << "HIGH ERROR" << endl;
                    cout << interpolated_value(0,0) << "," <<  true_potential_local(0,0) << endl;
                    cout << "distance_from_target: " << distance_from_target << endl;
                    cout << "expansion_spacing: " << expansion_spacing << endl;
                    cout << "singularity_dist: " << singularity_dist << endl;
                    cout << min_error << endl;
                    cout << endl << endl;
                    exit(0);
                }
                Result result = make_tuple(abs_error, R_over_rho, 
                        r_over_R, singularity_dist);
                error_results.push_back(result);


                VecDestroy(&interp_node_values);
                VecDestroy(&qbkix_points);
            }
        }
        DblNumMat error_and_params(4,error_results.size());
        DblNumVec index(error_results.size());
        int num_values = error_results.size();
        for(int i = 0; i < num_values; i++){
            index(i) = i;
            
            Result result = error_results[i];
            error_and_params(0,i) = get<0>(result);
            error_and_params(1,i) = get<1>(result);
            error_and_params(2,i) = get<2>(result);
            error_and_params(3,i) = get<3>(result);

        }
        write_to_file("qbkix_error_parameter_sweep.txt",
                error_and_params,
                index);

    }
}

TEST_CASE("Test qbkix extrapolation error w.r.t patch size ", "[results][qbkix-bbox-size]"){
    
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", ".1");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".066666");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "12");
        Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");

    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    Options::set_value_petsc_opts("-kt", "111");
    

    TestConfig test;
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type = SingularityType::AXIS_ALIGNED;
    
    //test.target_type   = TargetType::SINGLE;
    test.target_type   = TargetType::LINE;
    test.evaluation_scheme = EvaluationScheme::GREENS_IDENTITY;
    test.single_target_point = Point3(0.353553, -0.25, 0.353553);
    test.solution_scheme   = SolutionScheme::NONE;
    test.solver_matvec_type = NONE;
    test.compute_far_marking =false;
    
    int max_refinement_depth = 1;
    int num_iters = 10;

    vector<Result> error_results;
    error_results.reserve(max_refinement_depth*num_iters);
    double r = Options::get_double_from_petsc_opts("-interpolation_spacing_ratio");
    for(int ref_level = 0; ref_level < max_refinement_depth; ref_level++){

    //Options::set_value_petsc_opts("-boundary_distance_ratio", std::to_string(r));
        unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
        surface->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        surface->_coarse = true;
        surface->setFromOptions();
        surface->setup();
        surface->refine_uniform(ref_level);

        Vec true_potential;
        Vec computed_potential;
        Vec targets;
        Vec boundary_data;
        Vec solved_density;

        for(int Ri =0; Ri < num_iters; Ri++){
            // iterate over R from 0,..., r
            double inflation_factor = num_iters <= 1 ? 0 : double(Ri)/double(num_iters-1);

            Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", 
                    std::to_string(inflation_factor));
            int num_coarse_patches = surface->num_patches();
            string prefix = string("output/experiments/inflation_factor_");
            stats.add_result("num coarse patches", num_coarse_patches);
            stats.add_result("inflation factor", inflation_factor);
            string filename = prefix+to_string(inflation_factor)+string(".vtp");
            stats._file_prefix = prefix+to_string(inflation_factor);
            test.error_filename = filename;
            test.dump_values = true;

            auto solver = setup_solver(surface.get(), test.solver_matvec_type);

            NumVec<OnSurfacePoint> closest_points;
            cout << "SETTING TARGET DATA" << endl;
            setup_target_data(solver, test, targets, closest_points, 
                    computed_potential, true_potential);
            setup_problem_data(solver, test, boundary_data, solved_density);

            cout << "solution_scheme: " << int(test.solution_scheme) << endl;


            solve_and_evaluate(solver,
                    boundary_data,
                    targets,
                    closest_points, 
                    solved_density, 
                    computed_potential,
                    test);
            Vec error = compute_error_petsc_vec(true_potential, computed_potential);
            DblNumMat error_local = get_local_vector(1,Petsc::get_vec_size(error),error);
            cout << error_local << endl;
            tear_down(targets, true_potential, computed_potential, test, 1, 1e-6, 1e-6, test.dump_values);
            DblNumMat t = get_local_vector(DIM, Petsc::get_vec_size(targets)/DIM, targets);
            for (int i = 0; i < closest_points.m(); i++) {
                auto p = closest_points(i);
                if(i == 26){
                //if(p.region == FAR){
                    cout << "point " << i  << ": " << endl;
                    cout << Point3(t.clmdata(p.target_index))<< endl;
                    cout << "distance to target: " << p.distance_from_target << endl;
                    cout << "region: " << p.region<< endl;

                }
                
            }
            t.restore_local_vector();
            double abs_error = DBL_MAX;
            for(int i =0; i < Petsc::get_vec_size(error); i++){
                abs_error = max(abs_error, error_local(0,i));
            }
            auto patch = surface->patch(closest_points(0).parent_patch);
            double L = patch->characteristic_length();
            error_local.restore_local_vector();
            VecDestroy(&error) ;  

                Result result = make_tuple(abs_error, inflation_factor, 
                        L, double(ref_level));
            error_results.push_back(result);

        }
        VecDestroy(&true_potential);
        VecDestroy(&computed_potential);
        VecDestroy(&solved_density);
        VecDestroy(&boundary_data);
        VecDestroy(&targets);

    }
        DblNumMat error_and_params(4,error_results.size());
        DblNumVec index(error_results.size());
        int num_values = error_results.size();
        cout << error_results.size() << endl;
        for(int i = 0; i < num_values; i++){
            index(i) = i;
            
            Result result = error_results[i];
            error_and_params(0,i) = get<0>(result);
            error_and_params(1,i) = get<1>(result);
            error_and_params(2,i) = get<2>(result);
            error_and_params(3,i) = get<3>(result);

        }
    write_to_file("qbkix_error_patch_size_parameter_sweep.txt",
                error_and_params,
                index);

}

TEST_CASE("Random test", "[random]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");

    TestConfig test;
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type = SingularityType::AXIS_ALIGNED;

    test.target_type   = TargetType::SINGLE;
    test.evaluation_scheme = EvaluationScheme::GREENS_IDENTITY;
    test.single_target_point = Point3(0.353553, -0.25, 0.353553);
    test.solution_scheme   = SolutionScheme::NONE;
    test.compute_far_marking =false;



    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine_uniform(0);
    auto solver = setup_solver(surface.get(), test.solver_matvec_type);

}
