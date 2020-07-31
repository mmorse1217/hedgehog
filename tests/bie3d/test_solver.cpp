#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "../utils/evaluation_utils.hpp"
#include "common/vtk_writer.hpp"
using namespace hedgehog;

TEST_CASE("Test solver", "[solver][eval]"){
    SECTION("Test smooth quadrature with single layer kernel"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".5");
        //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".25");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
        Options::set_value_petsc_opts("-dom", "0");
        Options::set_value_petsc_opts("-bdtype", "2");


        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        vector<int> patch_partition(face_map->patches().size(), 0);  //All in one processor

        SolverGMRESDoubleLayer* solver = new SolverGMRESDoubleLayer("BIS3D_", "bis3d_");
        solver->bdry() = (PatchSurf*) face_map;
        solver->patch_partition() = patch_partition;
        solver->dom() = Options::get_int_from_petsc_opts("-dom");
        solver->set_evaluation_type(EXTRAPOLATION_AVERAGE);
        //solver->_compute_refined_surface = false;
        solver->eqcoefs() = vector<double>(2, 1.);
        solver->setFromOptions();
        solver->setup();

        Vec singularity_positions;
        Vec singularity_normals;
        Vec singularity_densities; 
        Vec boundary_data; 
        Kernel3d single_layer_kernel(solver->equation_type() + SINGLE_LAYER + VAR_U, solver->eqcoefs());
        Petsc::create_mpi_vec(solver->mpiComm(), 
                solver->patch_samples()->local_num_sample_points()*single_layer_kernel.get_sdof(), 
                boundary_data);

        axis_aligned_singularities(
                solver->mpiComm(),
                single_layer_kernel,
                singularity_positions,
                singularity_normals,
                singularity_densities);

        evaluate_singularity_solution(
                single_layer_kernel,
                singularity_positions,      
                singularity_normals,        
                singularity_densities,      
                solver->patch_samples()->sample_point_3d_position(),
                boundary_data);


        Vec targets;
        Petsc::create_mpi_vec(solver->mpiComm(), 1*DIM, targets);
        VecSet(targets, 0.);

        Vec true_potential = evaluate_singularities_along_basis(
                solver->mpiComm(),
                single_layer_kernel, 
                targets);

        //Vec computed_potential;
        //Petsc::create_mpi_vec(solver->mpiComm(), 1*single_layer_kernel.get_tdof(), computed_potential);

    Vec computed_potential = 
        greens_identity(solver->mpiComm(), 
                single_layer_kernel,
                singularity_positions,
                singularity_densities,
                targets,
                solver);
        //solver->fareval(targets, -9, boundary_data, computed_potential);

        VecView(true_potential, PETSC_VIEWER_STDOUT_SELF);
        VecView(computed_potential, PETSC_VIEWER_STDOUT_SELF);
        Petsc::destroy_vec(targets);
        Petsc::destroy_vec(computed_potential);
        Petsc::destroy_vec(true_potential);
    }

/*
    SECTION("Test general point evaluator constant density, parallel copies of surface"){
        
        // override options for particular test
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        //Options::set_value_petsc_opts("-bis3d_spacing", ".125");
        //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".06125");
        Options::set_value_petsc_opts("-near_interpolation_num_samples", "8");
        Options::set_value_petsc_opts("-bis3d_spacing", ".25");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
        Options::set_value_petsc_opts("-bis3d_np", "10");
    Options::set_value_petsc_opts("-dump_qbkix_points", "0");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
        Options::set_value_petsc_opts("-kt", "111"); // laplace

        // set up surface
        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test(); // skip refinement...

        // set up solver
        vector<int> patch_partition(face_map->patches().size(), 0);  //All in one processor

        SolverGMRESDoubleLayer* bis = new SolverGMRESDoubleLayer("BIS3D_", "bis3d_");
        bis->bdry() = (PatchSurf*) face_map;
        bis->patch_partition() = patch_partition;
        bis->dom() = Options::get_int_from_petsc_opts("-dom");
        bis->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    //bis->_compute_refined_surface = false;
        bis->eqcoefs() = vector<double>(2, 1.);
        bis->setFromOptions();
        bis->setup();

        vector<double> scale_factors;
        for(int i =1; i < 5; i++)
            scale_factors.push_back(1. - double(i)*.05);

        int num_samples = 20;
        double step = 1./double(num_samples-1);
        Vec sphere;
        Petsc::create_mpi_vec(bis->mpiComm(), num_samples*num_samples*DIM, sphere);
        DblNumMat sphere_local = get_local_vector(DIM, num_samples*num_samples, sphere);
        for(int i =0; i < num_samples; i++){
            for(int j =0; j < num_samples; j++){
                int index = i*num_samples +j;
                Point2 uv(i*step, j*step);
                double u = uv.x();
                double v = uv.y();
                double r = 1.;
                Point3 sphere_sample(
                        r*cos(u/(2.*M_PI))*sin(v/M_PI),
                        r*sin(u/(2.*M_PI))*sin(v/M_PI),
                        r*cos(v/M_PI));
                for(int d =0; d < DIM; d++)
                    sphere_local(d,index) = sphere_sample(d);
            }
        }
        sphere_local.restore_local_vector();

        Vec targets = 
            Test::create_scaled_surface_target_points(
                    scale_factors, 
                    sphere);
                    //bis->patch_samples()->sample_point_3d_position());
        //cout << "target points" << endl;
        //VecView(targets, PETSC_VIEWER_STDOUT_WORLD);
        int num_targets = Petsc::get_vec_size(targets)/DIM;

        Kernel3d kernel(bis->equation_type() + DOUBLE_LAYER + VAR_U, bis->eqcoefs());
        
        Vec density;
        Vec potential;
        VecCreateMPI(
                PETSC_COMM_WORLD,
                kernel.srcDOF()*bis->patch_samples()->local_num_sample_points(),
                PETSC_DETERMINE,
                &density);
        
        VecCreateMPI(
                PETSC_COMM_WORLD,
                kernel.trgDOF()*(Petsc::get_vec_size(targets)/DIM),
                PETSC_DETERMINE,
                &potential);


        // TODO fails when there are mixed near/far in a contiguous vector...
        NumVec<OnSurfacePoint> closest_on_surface_points(num_targets);
        cout << Petsc::get_vec_size(potential)/kernel.trgDOF() << ", " << Petsc::get_vec_size(targets)/DIM << endl; 
        cout << Petsc::get_vec_size(density)/kernel.srcDOF() << ", " << Petsc::get_vec_size(bis->patch_samples()->sample_point_blend_func_value()) << endl; 
        VecSet(density, 1.);
        bis->evaluate(targets,
                density,
                potential,
                closest_on_surface_points,
                VAR_U);

        DblNumMat potential_local = get_local_vector(kernel.trgDOF(),  num_targets, potential);
        for(int i =0; i < potential_local.n(); i++){
            OnSurfacePoint p = closest_on_surface_points(i);
            cout << potential_local(0,i) << endl;
            if(p.inside_domain == INSIDE){
                CHECK(fabs(potential_local(0,i) - 1.) <= 1e-3);
            } else if(p.inside_domain == OUTSIDE){
                CHECK(fabs(potential_local(0,i)) <= 1e-3);
            }
        }

        VecDestroy(&density);
        VecDestroy(&potential);
        delete bis;
    }*/
}

TEST_CASE("Test solver evaluation", "[green][eval]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "3");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "6");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".5");
    //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-3");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    double eval_time = omp_get_wtime();
    
    PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    //face_map->refine();
    stats.print_results(); 
    //face_map->refine_uniform(1);
    face_map->refine_test();

    
    SolverGMRESDoubleLayer* solver = new SolverGMRESDoubleLayer("BIS3D_", "bis3d_");
    solver->bdry() = face_map;
    //solver->_compute_refined_surface = false;
    int num_patches = face_map->patches().size();
    vector<int> patch_partition(num_patches, 0);  //All in one processor
    
    solver->patch_partition() = patch_partition;
    solver->dom() = Options::get_int_from_petsc_opts("-dom");
    solver->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    solver->eqcoefs() = vector<double>(2,1.0); 
    solver->_compute_refined_surface = true;
    
    solver->setFromOptions();
    solver->setup();

    Vec singularity_positions;
    Vec singularity_normals;
    Vec singularity_densities; 
    
    Kernel3d single_layer_kernel(solver->equation_type() + SINGLE_LAYER + VAR_U, solver->eqcoefs());
    axis_aligned_singularities(
        solver->mpiComm(),
        single_layer_kernel,
        singularity_positions,
        singularity_normals,
        singularity_densities);

    /*
    VecCreateMPI(solver->mpiComm(),
            solver->patch_samples()->local_num_sample_points()*single_layer_kernel.get_sdof(),
            PETSC_DETERMINE,
            &singularity_densities);
    VecCreateMPI(solver->mpiComm(),
            solver->patch_samples()->local_num_sample_points()*DIM,
            PETSC_DETERMINE,
            &singularity_normals);
    VecSet(singularity_normals, 0.);
    VecSet(singularity_densities, (4*M_PI));
    VecDuplicate(solver->patch_samples()->sample_point_3d_position(), &singularity_positions);
    VecCopy(solver->patch_samples()->sample_point_3d_position(),singularity_positions);
    VecScale(singularity_positions, 5.);
    */
    
    
    int target_dof = single_layer_kernel.get_tdof();
    Vec targets;
    VecDuplicate(solver->patch_samples()->sample_point_3d_position(), &targets);
    VecCopy(solver->patch_samples()->sample_point_3d_position(),targets);
    VecScale(targets, .999);
    int num_targets = Petsc::get_vec_size(targets)/DIM;


    //VecView(singularity_positions, PETSC_VIEWER_STDOUT_SELF);
    //cout << endl;
    //VecView(singularity_densities, PETSC_VIEWER_STDOUT_SELF);
    Vec computed_potential = 
        greens_identity(solver->mpiComm(), 
                single_layer_kernel,
                singularity_positions,
                singularity_densities,
                targets,
                solver);
    //VecView(computed_potential, PETSC_VIEWER_STDOUT_SELF);
 
    eval_time = omp_get_wtime() - eval_time;
        stats.add_result("total running time", eval_time);
        /*Vec true_potential = evaluate_singularities_along_basis(
                solver->mpiComm(),
            single_layer_kernel, 
            targets);*/
    Vec true_potential = Test::compute_dirichlet_boundary_data(
            solver->mpiComm(),
            single_layer_kernel,
            singularity_positions,
            singularity_densities,
            targets,
            NULL);
        /*
    Vec true_potential;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 
            num_targets*target_dof, 
            true_potential);
    
    evaluate_singularity_solution(
            single_layer_kernel,
            singularity_positions,      
            singularity_normals,        
            singularity_densities,     
            targets,              
            true_potential);          
            */
    DblNumMat computed_potential_local = get_local_vector(target_dof, num_targets, computed_potential);
    DblNumMat true_potential_local = get_local_vector(target_dof, num_targets, true_potential);
    for(int i=0; i < num_targets; i++){
        for(int d=0; d < target_dof; d++){
            double accurate = true_potential_local(d,i);
            double computed = computed_potential_local(d,i);
            cout << accurate << " ~= " << computed << endl;
            CHECK(fabs(accurate - computed)/fabs(accurate) <= 1e-5);
        }
    }
    true_potential_local.restore_local_vector();
    computed_potential_local.restore_local_vector();
        Vec error = compute_error_petsc_vec(true_potential, computed_potential);
        write_general_points_to_vtk(targets, target_dof, "laplace_err.vtp", error);
        VecDestroy(&error);
    stats.print_results();
    stats.clear();
    



}


TEST_CASE("weird parallel pvfmm test", "[mpi-pvfmm]"){

    int n = 6750;
    DblNumMat check_points(DIM, 6*n);
    DblNumMat refined_sample_points(DIM,4*n) ;
    DblNumMat refined_sample_points_normal(DIM,4*n) ;
    DblNumMat refined_density(1, 4*n);
    DblNumMat check_potential(1,6*n);
    
    Debug::load_mat(check_points, "output/check_point_positions_single_proc.mat");
    Debug::load_mat(refined_sample_points, "output/refined_sample_positions_single_proc.mat");
    Debug::load_mat(refined_sample_points_normal, "output/refined_sample_normals_single_proc.mat");
    Debug::load_mat(refined_density, "output/refined_density_single_proc.mat");
    Debug::load_mat(check_potential, "output/check_potential_single_proc.mat");


    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //assert(size ==3);
    Vec check_points_vec;
    Vec refined_sample_points_vec;
    Vec refined_sample_points_normals_vec;
    Vec refined_density_vec;
    Vec check_potential_vec;
    
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 6*n*DIM/size, check_points_vec);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 4*n*DIM/size, refined_sample_points_vec);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 4*n*DIM/size, refined_sample_points_normals_vec);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 4*n*1/size, refined_density_vec);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 6*n*1/size, check_potential_vec);
    DblNumMat check_points_local(DIM, check_points_vec);
    DblNumMat check_potential_local(1, check_potential_vec);
    DblNumMat refined_sample_points_local(DIM, refined_sample_points_vec) ;
    DblNumMat refined_sample_points_normals_local(DIM, refined_sample_points_normals_vec) ;
    DblNumMat refined_density_local(1, refined_density_vec);
    
    for(int d =0; d < DIM; d++){
        for(int i =0; i < check_points_local.n(); i++){
            int stride = 6*n/size;
            check_points_local(d,i) = check_points(d, rank*stride +i);
            if(d==0)
                check_potential_local(d,i) = check_potential(d, rank*stride +i);
        }
        for(int i =0; i < refined_sample_points_local.n(); i++){
            int stride = 4*n/size;
           refined_sample_points_local(d,i) = refined_sample_points(d, rank*stride +i);
           refined_sample_points_normals_local(d,i) = refined_sample_points_normal(d, rank*stride +i);
           if(d==0)
               refined_density_local(d,i) = refined_density(d, rank*stride +i);
        }
    }

    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");
    Kernel3d k(121, vector<double>(2,1.));
    auto fmm = unique_ptr<FMM>(new PvFMM(refined_sample_points_vec,
            refined_sample_points_normals_vec,
            check_points_vec, k));
    Vec parallel_check_potential_vec;
    VecDuplicate(check_potential_vec, &parallel_check_potential_vec);
    VecSet(parallel_check_potential_vec,0.);
    fmm->evaluate(refined_density_vec, parallel_check_potential_vec);

    Vec error;
    VecDuplicate(check_potential_vec, &error);
    VecCopy(check_potential_vec, error);
    VecAXPY(error, -1., parallel_check_potential_vec);
    double norm;
    VecAbs(error);
    VecNorm(error, NORM_INFINITY, &norm);

    cout << "max error: " << norm << endl;


}

