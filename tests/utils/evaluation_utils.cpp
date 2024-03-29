#include "evaluation_utils.hpp"
#include "common/vtk_writer.hpp"
#include "common/nummat.hpp"
#include "common/stats.hpp"
#include <petscvec.h>
#include <petscviewer.h>
#include <sampling.hpp>
#include "bie3d/markgrid.hpp"
#include "bie3d/solver_utils.hpp"
#include "../catch.hpp"
#include "bie3d/evaluator_qbkix.hpp"
using Sampling::sample_3d;
using Sampling::sample_2d;
using Sampling::sample_1d;
using Sampling::equispaced;

BEGIN_EBI_NAMESPACE

void setup_problem_data_constant_density(
        unique_ptr<SolverGMRESDoubleLayer>& solver, 
        Vec& boundary_data, Vec& solved_density){
    int sample_dof; // # dof from density
    int pole_dof; // # dof from multiply connected doms + poles
    int total_num_dof; // sample_dof+ pole_dof
    solver->localSize(sample_dof,pole_dof, total_num_dof);
    
    Kernel3d kernel = solver->problem_kernel();
    int source_dof = kernel.get_sdof();
    int target_dof= kernel.get_tdof();
    int num_samples = sample_dof/source_dof;
    
    Petsc::create_mpi_vec(solver->mpiComm(),
            //target_dof*total_num_dof,
            total_num_dof,
            boundary_data);
    VecSet(boundary_data, 1.);
    
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof, 
            solved_density);
    VecSet(solved_density, 1.);

}

void setup_target_data_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver, 
        double true_potential_value,
        Vec& targets, Vec& computed_potential, Vec& true_potential){
    if(targets == nullptr){
        targets = generate_targets_in_cube(25,
                Cube(Interval(-.2, .2),Interval(-.2, .2),Interval(-.2, .2)));
    }
    
    int num_target_points = Petsc::get_vec_local_size(targets)/DIM;
    Kernel3d kernel = solver->problem_kernel();
    int target_dof= kernel.get_tdof();
    Petsc::create_mpi_vec(solver->mpiComm(),
            num_target_points*target_dof,
            computed_potential);
    
    VecDuplicate(computed_potential, &true_potential);
    VecSet(true_potential, true_potential_value);

}


void tear_down(Vec targets, Vec true_potential, Vec computed_potential, TestConfig test,
        int target_dof,
        double eps_abs, double eps_rel, bool dump_values){
    int num_target_points = Petsc::get_vec_local_size(true_potential)/target_dof;
    DblNumMat true_potential_local = get_local_vector(target_dof, num_target_points, true_potential);
    DblNumMat computed_potential_local = get_local_vector(target_dof, num_target_points, computed_potential);
    
    check_error(true_potential_local, computed_potential_local, eps_abs, eps_rel);

    true_potential_local.restore_local_vector();
    computed_potential_local.restore_local_vector();

    if(dump_values){
        Vec error = compute_error_petsc_vec(true_potential, computed_potential);
        double max_error = 1e-16;
        VecMax(error, NULL, &max_error);
        Vec potential_abs;
        VecDuplicate(true_potential, &potential_abs);
        VecCopy(true_potential, potential_abs);
        VecAbs(potential_abs);
        double max_potential = 1e-16;
        VecMax(potential_abs, NULL, &max_potential);
        stats.add_result("max absolute error", max_error);
        stats.add_result("max relative error", max_error/max_potential);
        cout << "print to file : " << test.error_filename << endl;
        write_general_points_to_vtk(targets, target_dof, test.error_filename, error);
        if(test.target_type == TargetType::PLANE){
            mesh_and_save_points(DblNumMat(DIM, targets), DblNumMat(target_dof, computed_potential), stats._file_prefix +"_solution");
            mesh_and_save_points(DblNumMat(DIM, targets), DblNumMat(1, error), stats._file_prefix +"_error");

        }
        VecDestroy(&error);
        VecDestroy(&potential_abs);

    }
    


}

void check_error(DblNumMat true_values, DblNumMat computed_values, double eps_abs, double eps_rel){
    assert(true_values.m() == computed_values.m());
    assert(true_values.n() == computed_values.n());

    for(int ti =0; ti < true_values.n(); ti++){
        for(int d =0; d < true_values.m(); d++){
            double true_value = true_values(d,ti);
            double computed_value = computed_values(d,ti);
            cout << true_value <<  ", " << computed_value << ", abs: "
            << fabs(true_value- computed_value) << ", rel: " 
            << fabs(true_value- computed_value) /fabs(true_value) << endl;
            CHECK(fabs(true_value- computed_value) <= fabs(true_value)*eps_rel + eps_abs);
            if(fabs(true_value- computed_value) >=  1e-7){
                cout << "index: " << ti << endl;
            }
        }
    }
}

Vec compute_error_petsc_vec(Vec true_potential, Vec computed_potential){
    Vec error;
    VecDuplicate(true_potential, &error);
    VecCopy(true_potential, error);
    int minus_one = -1.;
    VecAXPY( error, minus_one,  computed_potential);
    VecAbs(error);
    DblNumMat error_local = get_local_vector(1, Petsc::get_vec_local_size(error), error);
    for(int i = 0; i < Petsc::get_vec_local_size(error); i++){
        error_local(0,i) = error_local(0,i) < 1e-16 ? 1e-16 : error_local(0,i);
    }
    error_local.restore_local_vector(); 
    return error;
}

Vec generate_targets_on_sphere(int num_samples_1d, double sphere_radius){
    int num_samples = num_samples_1d*num_samples_1d;
    Vec target_points;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_samples*DIM, target_points);
    
    DblNumMat target_points_local = 
        get_local_vector(DIM, num_samples, target_points);
    
    Rectangle theta_phi_domain(Interval(0,2*M_PI), Interval(0, M_PI));
    DblNumMat theta_phi_samples(2, num_samples, true,
        sample_2d<equispaced>(num_samples_1d, theta_phi_domain).data());

    for(int i =0; i < num_samples; i++){
        double theta = theta_phi_samples(0,i);
        double phi = theta_phi_samples(1,i);
        Point3 sphere_sample(
                sin(theta)*cos(phi),
                sin(theta)*sin(phi),
                cos(theta));
        sphere_sample *= sphere_radius;

        for(int d =0; d < DIM; d++)
            target_points_local(d,i) = sphere_sample(d);
    }
    target_points_local.restore_local_vector();
    return target_points;
}

void setup_singularities(TestConfig test, unique_ptr<SolverGMRESDoubleLayer>& solver,
        Vec& positions, Vec& strengths){

    Vec normals;
    Kernel3d k = solver->problem_kernel();
    DblNumMat positions_local;
    switch(test.singularity_type){
        case SingularityType::AXIS_ALIGNED:
             axis_aligned_singularities(
                    solver->mpiComm(),
                    k,
                    positions,
                    normals,
                    strengths);
             VecDestroy(&strengths);
             strengths = Test::generate_random_vector(
                     6*k.get_sdof(), 0., 1.);
             break;
        case SingularityType::SPHERE:
             positions = generate_targets_on_sphere(15, test.sphere_radius_bc);
             strengths = Test::generate_random_vector(
                     Petsc::get_vec_local_size(positions)/DIM*k.get_sdof(), 0., 1.);
             break;

        case SingularityType::SINGLE:
             Petsc::create_mpi_vec(solver->mpiComm(), 3*1, positions);
             Petsc::create_mpi_vec(solver->mpiComm(), k.get_sdof()*1, strengths);
             VecSet(strengths, 4*M_PI);
             positions_local = get_local_vector(1,3,positions);
             
             positions_local(0,0) = test.single_singularity_location(0);
             positions_local(0,1) = test.single_singularity_location(1);
             positions_local(0,2) = test.single_singularity_location(2);

             positions_local.restore_local_vector();
             break;
        default:
             assert(0);
    }
    
}


Vec generate_targets_in_cube(int num_samples_1d, Cube target_domain){
    int num_samples = num_samples_1d*num_samples_1d*num_samples_1d;
    Vec target_points;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_samples*DIM, target_points);
    
    DblNumMat target_points_local = 
        get_local_vector(DIM, num_samples, target_points);
    
    DblNumMat volume_samples(DIM, num_samples, true, 
            sample_3d<equispaced>(num_samples_1d, target_domain).data());
    for(int i =0; i < num_samples; i++){
        for(int d =0; d < DIM; d++){
            target_points_local(d,i) = volume_samples(d,i);
        }

    }
    target_points_local.restore_local_vector();

    return target_points;
}


unique_ptr<SolverGMRESDoubleLayer> setup_solver(PatchSurf* surface, 
        EvaluationType eval_type, bool compute_refined_surface){

    
    unique_ptr<SolverGMRESDoubleLayer> solver( new SolverGMRESDoubleLayer(surface));
    solver->set_evaluation_type(eval_type);
    solver->_compute_refined_surface = compute_refined_surface;
    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4;
    
    solver->setFromOptions();
    solver->setup(); // memory leak <<==========!!!!!
    return solver;
}

Vec setup_plane_targets(TestConfig test, MPI_Comm comm){
    Vec targets;
    int n = test.num_targets;
    // Sample a plane centered at point p spanned by vectors v1 and v2
    // parameter values are sampled from [-1,1]^2. only handle "square" sampling
    // of plane; for [-a, a], rescale vectors: v1 *= a, v2 *= b
    Petsc::create_mpi_vec(comm, n*n*DIM, targets);
    DblNumMat targets_local(DIM, targets);
    Point3 p = test.target_plane_point;
    Point3 v1 = test.target_plane_vec1;
    Point3 v2 = test.target_plane_vec2;
    double step = 2./double(n-1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double ti = i*step-1;
            double tj = j*step-1;

            Point3 t = p + ti*v1 + tj*v2;
            for(int d =0; d < DIM; d++) {
                targets_local(d,i*n+j) = t(d);
            }
        }
    }

    targets_local.restore_local_vector();
    return targets;
}

Vec setup_line_targets(TestConfig test, MPI_Comm comm){
    Vec targets;
    int num_targets =200;
    Petsc::create_mpi_vec(comm, num_targets*DIM, targets);
    DblNumMat targets_local = get_local_vector(DIM, num_targets, targets);
    Point3 init(-.9, 0., 0.);
    double step = 1.8/double(num_targets-1);
    for (int i = 0; i < num_targets; i++) {
        Point3 t =init+Point3(i*step, 0., 0.);
        for(int d =0; d < DIM; d++) {
            targets_local(d,i) = t(d);
        }
    }

    targets_local.restore_local_vector();
    VecView(targets, PETSC_VIEWER_STDOUT_SELF);
    return targets;
}

Vec set_single_target(TestConfig test, MPI_Comm comm){
    Vec targets;
    Petsc::create_mpi_vec(comm, 1*DIM, targets);
    VecView(targets, PETSC_VIEWER_STDOUT_SELF);
    DblNumMat targets_local = get_local_vector(DIM, 1, targets);
    for(int d =0; d < DIM; d++) {
        targets_local(d,0) = test.single_target_point(d);
    }

    targets_local.restore_local_vector();
    return targets;
}

Vec setup_ring_targets(int num_samples_1d, double ring_radius){
    int num_samples = num_samples_1d;
    Vec target_points;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_samples*DIM, target_points);
    
    DblNumMat target_points_local = 
        get_local_vector(DIM, num_samples, target_points);
    
    Rectangle theta_domain(Interval(0,2*M_PI), Interval(0, 2*M_PI));
    DblNumMat theta_samples(1, num_samples, true,
        sample_1d<equispaced>(num_samples_1d, theta_domain).data());

    for(int i =0; i < num_samples; i++){
        double theta = theta_samples(0,i);
        Point3 ring_sample(
                sin(theta),
                0.,
                cos(theta));
        ring_sample *= ring_radius;

        for(int d =0; d < DIM; d++)
            target_points_local(d,i) = ring_sample(d);
    }
    target_points_local.restore_local_vector();
    return target_points;
}

Vec set_plane_targets(TestConfig test, MPI_Comm comm){
    Vec targets;
    Petsc::create_mpi_vec(comm, test.num_targets*DIM, targets);
    DblNumMat targets_local = get_local_vector(DIM, test.num_targets, targets);
    for(int i = 0; i < test.num_targets; i++){
        for(int j = 0; j < test.num_targets; j++){
            int idx = i*test.num_targets +j ;
            double step = 2./(test.num_targets-1);
            
            Point3 p = test.target_plane_point;
            Point3 v1 = test.target_plane_vec1;
            Point3 v2 = test.target_plane_vec2;

            Point3 t = p + (i*step - 1.)*v1 +  (j*step - 1.)*v2;

            for(int d =0; d < DIM; d++) {
                targets_local(d,i) = t(d);
            }
        }
    }

    targets_local.restore_local_vector();
    return targets;

}


Vec compute_analytic_boundary_data(MPI_Comm comm,
        AnalyticSolutionType analytic_solution_type,
        Vec target_positions){
    Vec boundary_data;
    return boundary_data;
}
Vec setup_on_surface_targets(unique_ptr<SolverGMRESDoubleLayer>& solver){
    Vec targets;
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    // coarse targets that are too coarse to bother solving on
    Options::set_value_petsc_opts("-bis3d_spacing",
            to_string(sqrt(2.)*spacing));


    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = solver->patch_partition();
    samples->bdry() = solver->bdry();
    samples->setup();

    VecDuplicate(samples->sample_point_3d_position(), &targets);
    VecCopy(samples->sample_point_3d_position(), targets);

    // reset spacing back to original so we don't break anything
    Options::set_value_petsc_opts("-bis3d_spacing", to_string(spacing));
    return targets;
}



void set_target_positions(unique_ptr<SolverGMRESDoubleLayer>& solver,
        TestConfig& test, Vec& targets){

    //Interval side(-.3,.3);
    Interval side(.2,.8);
    int num_targets = test.num_targets;
    DblNumMat targets_local;
    cout << "target_type: " << int(test.target_type) << endl;
    switch(test.target_type){
        case TargetType::GRID:

            // generate num_targets^3 points on a uniform grid in the cubic volume
            // centered at the origin of with side as each dimension
            targets = generate_targets_in_cube(num_targets, Cube(side, side, side));

            break;
        case TargetType::COLLOCATION_POINTS:
            // targets are the same collocations points, copy them to be safe
            VecDuplicate(solver->patch_samples()->sample_point_3d_position(),
                    &targets);
            VecCopy(solver->patch_samples()->sample_point_3d_position(),
                    targets);
            //VecScale(targets, .99999);
            break;
        case TargetType::SPHERE_FAR:
            targets = generate_targets_on_sphere(num_targets, test.sphere_radius_targets);
            break;
        case TargetType::SINGLE:
            targets = set_single_target(test, solver->mpiComm());
            break;
        case TargetType::LINE:
            targets = setup_line_targets(test, solver->mpiComm());
            break;
        case TargetType::RING:
            targets = setup_ring_targets(num_targets*num_targets, test.sphere_radius_targets);
            break;
        case TargetType::ON_SURFACE_NON_COLLOCATION:
            targets = setup_on_surface_targets(solver);
            break;
        case TargetType::PLANE:
            targets = setup_plane_targets(test, solver->mpiComm());
            break;
        default:
            assert(0);
            break;

    }
}
void set_target_potential(unique_ptr<SolverGMRESDoubleLayer>& solver,
        TestConfig test_type, Vec targets, 
        NumVec<OnSurfacePoint> closest_on_surface_points,
        Vec& potential){
    double true_potential_value;
    Vec singularity_positions;
    Vec singularity_normals;
    Vec singularity_densities;
    
    int num_target_points = Petsc::get_vec_local_size(targets)/DIM;

    Kernel3d kernel = solver->problem_kernel();
    int target_dof= kernel.get_tdof();

    Petsc::create_mpi_vec(solver->mpiComm(),
            target_dof*num_target_points, 
            potential);

    DblNumMat potential_local = get_local_vector(
            target_dof, num_target_points, potential);

    BoundaryDataType bc_type= test_type.bc_type;
    cout << "targets in/out inside set_target_potential" << endl;
    switch(bc_type){
        case BoundaryDataType::CONSTANT_DENSITY:
            // For constant density, a target point x \in R^3 will have
            // potential 1 if x \in \Omega, 1/2 if x \in \partial \Omega, and 0
            // if x \not \in \Omega.

            if(Options::get_int_from_petsc_opts("-bdtype") == 2){
                assert(closest_on_surface_points.m() == num_target_points);
                bool is_interior_problem = Options::get_int_from_petsc_opts("-dom") == 0;
                for(int i =0 ; i <closest_on_surface_points.m(); i++){
                    OnSurfacePoint on_surface_point = closest_on_surface_points(i);
                    if(on_surface_point.region == NEAR)
                    assert(i == on_surface_point.target_index);
                    if(is_interior_problem){
                        if(on_surface_point.inside_domain == INSIDE || on_surface_point.inside_domain == ON_SURFACE){
                            true_potential_value = 1.;
                        } else if (on_surface_point.inside_domain == OUTSIDE){
                            true_potential_value = 0.;
                        } else {
                            assert(0);
                        }
                    } else{
                        if(on_surface_point.inside_domain == INSIDE|| on_surface_point.inside_domain == ON_SURFACE){
                            true_potential_value = 0.;
                        } else if (on_surface_point.inside_domain == OUTSIDE){
                            true_potential_value = -1.;
                        } else {
                            assert(0);
                        }
                    }

                    for(int d =0; d < target_dof; d++){
                        potential_local(d, i) = 
                            true_potential_value;
                    }

                }
            potential_local.restore_local_vector();
            } else {
                potential_local.restore_local_vector();
                VecSet(potential, 1.);
            }

            break;
        case BoundaryDataType::HARMONIC:
            setup_singularities(test_type, solver, singularity_positions, 
                    singularity_densities);

            potential  = Test::compute_dirichlet_boundary_data(
                    solver->mpiComm(),
                    kernel,
                    singularity_positions,
                    singularity_densities,
                    targets,
                    targets);
            VecDestroy(&singularity_positions);
            VecDestroy(&singularity_densities);
            
            break;
        case BoundaryDataType::RANDOM:
            assert(0); // do we need this?
            break;
        case BoundaryDataType::SPHERICAL_HARMONIC:
            assert(test_type.target_type == TargetType::COLLOCATION_POINTS );
            assert(Options::get_string_from_petsc_opts("-bd3d_meshfile") == 
                    string("wrl_meshes/wrl/sphere.wrl"));
            potential = Test::compute_spherical_harmonic_density(solver->mpiComm(), targets);
            break;
        default:
            assert(0);
            break;

    }
}

void setup_boundary_data(unique_ptr<SolverGMRESDoubleLayer>& solver,
        TestConfig test, Vec& boundary_data){

    Kernel3d kernel = solver->problem_kernel();
    int source_dof= kernel.get_sdof();
    
    Vec singularity_positions;
    Vec singularity_normals;
    Vec singularity_densities;
    Vec boundary_data_values;

    switch(test.bc_type){
        case BoundaryDataType::CONSTANT_DENSITY:
            Petsc::create_mpi_vec(solver->mpiComm(),
                    source_dof*solver->patch_samples()->local_num_sample_points(),
                    boundary_data_values);
           

            VecSet(boundary_data_values, 1.);
            break;
        case BoundaryDataType::HARMONIC:
            setup_singularities(test, solver, singularity_positions, singularity_densities);
            
            boundary_data_values = Test::compute_dirichlet_boundary_data(
                    solver->mpiComm(),
                    kernel,
                    singularity_positions,
                    singularity_densities,
                    solver->patch_samples()->sample_point_3d_position(),
                    solver->patch_samples()->sample_point_normal());
            
            break;
        case BoundaryDataType::RANDOM:
            assert(0); // do we need this?
        case BoundaryDataType::SPHERICAL_HARMONIC:
            assert(test.target_type == TargetType::COLLOCATION_POINTS );
            assert(Options::get_string_from_petsc_opts("-bd3d_meshfile") == 
                    string("wrl_meshes/wrl/sphere.wrl"));

            boundary_data_values = Test::compute_spherical_harmonic_bc(solver->mpiComm(), 
                    solver->patch_samples()->sample_point_3d_position());
            break;
        case BoundaryDataType::ANALYTIC:
            boundary_data_values = 
                compute_analytic_boundary_data(solver->mpiComm(), 
                    test.analytic_solution_type,
                    solver->patch_samples()->sample_point_3d_position());
            break;
        default:
            assert(0);
            break;

    }
    Petsc::create_mpi_vec(solver->mpiComm(), 
            solver->local_total_dof(),
            boundary_data);
    {
        DblNumVec bd(boundary_data);
        DblNumVec bd_values(boundary_data_values);
        assert(bd.m() >= bd_values.m());
        for (int i = 0; i < bd_values.m(); i++) {
            bd(i) = bd_values(i);
        }
            
    }
}

void setup_target_data(unique_ptr<SolverGMRESDoubleLayer>& solver, 
        TestConfig test_type, Vec& targets, 
    NumVec<OnSurfacePoint>& closest_points,
        Vec& computed_potential, Vec& true_potential){

    // determine target positions
    set_target_positions(solver, test_type, targets);
    

    // mark targets as in/out of problem domain
    int num_target_points = Petsc::get_vec_local_size(targets)/DIM;
    DblNumMat targets_local = get_local_vector(DIM, num_target_points, targets);
    if(auto f = dynamic_cast<PatchSurfFaceMap*>(solver->bdry())){
        cout << "marking targets" << endl;
        // all are near, and precomputed in patch samples
        if(test_type.target_type == TargetType::COLLOCATION_POINTS){
            closest_points = solver->patch_samples()->sample_point_as_on_surface_point();

        // all are far, just use smooth quad
        } else if(test_type.target_type == TargetType::SPHERE_FAR ||
                  test_type.target_type == TargetType::RING){
            closest_points.resize(num_target_points);
            for(int i = 0; i < closest_points.m(); i++){
                OnSurfacePoint p;
                p.distance_from_target = DBL_MAX;
                p.inside_domain = INSIDE;
                p.region = FAR;
                p.target_index = i;
                closest_points(i) = p;
            }

        } else if(test_type.target_type == TargetType::ON_SURFACE_NON_COLLOCATION){
            double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
            // coarse targets that are too coarse to bother solving on
            Options::set_value_petsc_opts("-bis3d_spacing",
                    to_string(sqrt(2.)*spacing));


            unique_ptr<PatchSamples> samples(new PatchSamples("",""));
            samples->patch_partition() = solver->patch_partition();
            samples->bdry() = solver->bdry();
            samples->setup();
            closest_points = samples->sample_point_as_on_surface_point();

            // reset spacing back to original so we don't break anything
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(spacing));

        // otherwise, explicitly mark
        } else {
        closest_points = 
            Markgrid::mark_target_points(targets_local,f,test_type.compute_far_marking);
        }
    } else {
        closest_points.resize(num_target_points);
    }

    targets_local.restore_local_vector();
    
    // determine true potential at targets 
    set_target_potential(solver, test_type, targets, closest_points, 
            true_potential);
    
    // computed potential should be the same size as true_potential; zero it out
    // initially
    VecDuplicate(true_potential, &computed_potential);
    VecSet(computed_potential, 0.);

}


void setup_problem_data(unique_ptr<SolverGMRESDoubleLayer>& solver, 
        TestConfig test_type, Vec& boundary_data, Vec& solved_density){

    // compute boundary conditions
    setup_boundary_data(solver, test_type, boundary_data);
    
    // initialize density vector 
    VecDuplicate(boundary_data, &solved_density);
    VecSet(solved_density, 0.);

}



void test_constant_boundary_data(PatchSurf* surface, bool dump_values){

    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    auto solver = setup_solver(surface, EXTRAPOLATION_AVERAGE);
    //test_far_eval_constant_density(solver, dump_values);
    //test_singular_quadrature_constant_density(solver, dump_values);
    //test_solver_constant_density(solver, dump_values);
    test_qbkix_constant_density(solver, dump_values);

    
}

void solve_and_evaluate(unique_ptr<SolverGMRESDoubleLayer>& solver,
        Vec boundary_data, Vec targets, NumVec<OnSurfacePoint> closest_points,
        Vec& solved_density, Vec& computed_potential,
        TestConfig test_type){

    Vec singularity_positions;
    Vec singularity_normals;
    Vec singularity_densities; 
    Vec neumann_data;

    Kernel3d single_layer_kernel(solver->equation_type() + SINGLE_LAYER + VAR_U, solver->eqcoefs());
    switch(test_type.solution_scheme){
        case SolutionScheme::EXPLICIT_DENSITY:
            // if we have constant density and don't want to solve, just evaluate
            // directly using the constant boundary data
            //solved_density = boundary_data;
            VecCopy(boundary_data, solved_density);
            break;
        case SolutionScheme::GMRES_SOLVE:
            // solve the pde given the dirichlet boundary data
            solver->solve(boundary_data, solved_density);
            break;
        default:
            break;
            //assert(0);
    } 
    
    Vec target_in_out;
    Vec interpolation_directions;
    Vec target_as_face_point;
    Vec closest_points_position;
    Vec closest_point_as_face_point;
    Vec target_far_field;
    Vec target_interpolant_spacing;
    

    Petsc::create_mpi_vec(solver->mpiComm(), 
            solver->patch_samples()->local_num_sample_points(),
            target_in_out);
        double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
        Options::set_value_petsc_opts("-bis3d_spacing", 
            to_string(sqrt(2.)*spacing));

    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = solver->patch_partition();
    samples->bdry() = solver->bdry();
    samples->setup();
    // reset spacing back to original so we don't break anything
    Options::set_value_petsc_opts("-bis3d_spacing", to_string(spacing));
    if(dynamic_cast<PatchSurfFaceMap*>(solver->bdry())){
        if(test_type.target_type == TargetType::COLLOCATION_POINTS || 
                test_type.target_type == TargetType::ON_SURFACE_NON_COLLOCATION){
            if(test_type.target_type == TargetType::ON_SURFACE_NON_COLLOCATION){
                VecDuplicate(samples->sample_as_face_point(), &target_as_face_point);
                VecCopy(samples->sample_as_face_point(), target_as_face_point);
            } else {
                VecDuplicate(solver->patch_samples()->sample_as_face_point(), &target_as_face_point);
                VecCopy(solver->patch_samples()->sample_as_face_point(), target_as_face_point);
            }
            solver->populate_qbx_data(closest_points,
                    target_in_out, 
                    closest_points_position,
                    closest_point_as_face_point,
                    interpolation_directions,
                    target_far_field,
                    target_interpolant_spacing);
            cout << "num closest points: " << closest_points.m() << endl;
            cout << "num 'targets': " << Petsc::get_vec_local_size(target_far_field) << endl;

        } 
    } else {
        if(test_type.target_type == TargetType::COLLOCATION_POINTS){
            VecDuplicate(solver->patch_samples()->sample_point_3d_position(), &closest_points_position);
            VecCopy(solver->patch_samples()->sample_point_3d_position(), closest_points_position);
            VecDuplicate(solver->patch_samples()->sample_as_face_point(), &target_as_face_point);
            VecCopy(solver->patch_samples()->sample_as_face_point(), target_as_face_point);
            VecDuplicate(solver->patch_samples()->sample_as_face_point(), &closest_point_as_face_point);
            VecCopy(solver->patch_samples()->sample_as_face_point(), closest_point_as_face_point);

            VecDuplicate(solver->patch_samples()->sample_point_normal(), &interpolation_directions);
            VecCopy(solver->patch_samples()->sample_point_normal(), interpolation_directions);
            VecDuplicate(solver->patch_samples()->sample_point_far_field(), &target_far_field);
            VecCopy(solver->patch_samples()->sample_point_far_field(), target_far_field);
            VecDuplicate(solver->patch_samples()->sample_point_interpolant_spacing(), &target_interpolant_spacing);
            VecCopy(solver->patch_samples()->sample_point_interpolant_spacing(), target_interpolant_spacing);

            VecScale(interpolation_directions, -1.);
        }else if (test_type.target_type == TargetType::ON_SURFACE_NON_COLLOCATION){
            VecDuplicate(samples->sample_point_3d_position(), &closest_points_position);
            VecCopy(samples->sample_point_3d_position(), closest_points_position);
            VecDuplicate(samples->sample_as_face_point(), &target_as_face_point);
            VecCopy(samples->sample_as_face_point(), target_as_face_point);
            VecDuplicate(samples->sample_as_face_point(), &closest_point_as_face_point);
            VecCopy(samples->sample_as_face_point(), closest_point_as_face_point);

            VecDuplicate(samples->sample_point_normal(), &interpolation_directions);
            VecCopy(samples->sample_point_normal(), interpolation_directions);
            VecDuplicate(samples->sample_point_far_field(), &target_far_field);
            VecCopy(samples->sample_point_far_field(), target_far_field);
            VecDuplicate(samples->sample_point_interpolant_spacing(), &target_interpolant_spacing);
            VecCopy(samples->sample_point_interpolant_spacing(), target_interpolant_spacing);

            VecScale(interpolation_directions, -1.);
        }
    }



    VecSet(target_in_out, 1);
Vec interpolated_density;
    DblNumMat potential_local;
    DblNumMat density_local;
    auto kernel = solver->problem_kernel();
    int num_local_targets = Petsc::get_vec_local_size(targets)/DIM;
    DblNumMat neumann_local;
    switch(test_type.evaluation_scheme){
        case EvaluationScheme::GREENS_IDENTITY:
            // TODO decouple
            cout << "GREENS ID EVAL" << endl;
            if(test_type.bc_type == BoundaryDataType::HARMONIC){

                // evaluate neumann data
                setup_singularities(
                        test_type,solver, 
                        singularity_positions,
                        singularity_densities);

                VecDuplicate(boundary_data, &neumann_data);
                neumann_local = DblNumMat(kernel.get_tdof(),neumann_data);

                kernel.neumann_bc_from_singularities(
                        DblNumMat(3,singularity_positions),
                        DblNumMat(kernel.get_sdof(),singularity_densities),
                        DblNumMat(3,solver->patch_samples()->sample_point_3d_position()), 
                        DblNumMat(3,solver->patch_samples()->sample_point_normal()), 
                        neumann_local);
                {
                    // make check points
                    int qbkix_id= Options::get_int_from_petsc_opts("-near_interpolation_num_samples");

                    vector<int> qbkix_ids;
                    for(int i = 0; i < qbkix_id; i++)
                        qbkix_ids.push_back(i);
                    
                    Vec id;
                    VecDuplicate(solver->patch_samples()->sample_point_normal(), &id);
                    VecCopy(solver->patch_samples()->sample_point_normal(), id);
                    VecScale(id,-.1);
                    Vec qbkix_points = 
                        generate_interior_qbkix_points(
                                solver->patch_samples()->sample_point_3d_position(),
                                solver->patch_samples()->sample_point_3d_position(),
                                qbkix_ids, 
                                id, 
                                solver->patch_samples()->sample_point_far_field(),
                                solver->patch_samples()->sample_point_interpolant_spacing());
                    // evaluate dirichlet data at check points
                    Vec potential_at_check_points  = Test::compute_dirichlet_boundary_data(
                            solver->mpiComm(),
                            kernel,
                            singularity_positions,
                            singularity_densities,
                            qbkix_points,
                            qbkix_points);
                    
                    VecDestroy(&qbkix_points);
                    VecDestroy(&id);
                    VecDestroy(&potential_at_check_points);

                }
               
                // pass closest points that were computed previously
                computed_potential = 
                    greens_identity(solver->mpiComm(), 
                            boundary_data,
                            neumann_data,
                            closest_points,
                            targets,
                            solver.get());
                
            } else if (test_type.bc_type == BoundaryDataType::CONSTANT_DENSITY){
                VecDuplicate(boundary_data, &neumann_data);
                VecSet(neumann_data, 0.);

                computed_potential = 
                    greens_identity(solver->mpiComm(), 
                            boundary_data,
                            neumann_data,
                            closest_points,
                            targets,
                            solver.get());

            }
            break;
            // off-surface quadratures:
            //  - smooth quadrature
        case EvaluationScheme::SMOOTH_QUAD:

            solver->fareval(targets, VAR_U, solved_density, computed_potential);
            break;

            //  - near-singular quadrature
        case EvaluationScheme::NEAR_SINGULAR:
            break;
        case EvaluationScheme::NEAR_QBKIX:
            break;

            //  - singular quadrature: on-surface
        case EvaluationScheme::ON_SINGULAR:
            cout << "num targets: " << Petsc::get_vec_local_size(targets)/DIM << endl;
            cout << "num targets fp: " << Petsc::get_vec_local_size(target_as_face_point)/(solver->bdry()->face_point_size_in_doubles()) << endl;
            cout << "num sources: " << Petsc::get_vec_local_size(solved_density)<< endl;
            cout << "num potential: " << Petsc::get_vec_local_size(computed_potential)<< endl;
            VecSet(computed_potential, 0.);
            solver->roneval(targets,
                    target_as_face_point, // wrong for general targets!
                    VAR_U, solved_density, computed_potential);
            

            // compute interior value
           VecAXPY(computed_potential, .5, solved_density);
            
            break;
        case EvaluationScheme::ON_QBKIX:

            solver->neaeval_extrapolate(
                    targets,  
                    NULL,
                    target_in_out, // whether the Omega_2 targets are in/out
                    closest_points_position,
                    closest_point_as_face_point,
                    target_far_field,
                    target_interpolant_spacing,
                    //solver->patch_samples()->sample_point_3d_position(),
                    //solver->patch_samples()->sample_as_face_point(),
                    interpolation_directions,
                    VAR_U,
                    solved_density, 
                    computed_potential);
            break;
        case EvaluationScheme::ON_QBKIX_AVERAGE:
            //assert(0); // fix for general targets
            solver->neaeval_extrapolate_average(
                    targets,  
                    NULL,
                    target_in_out, // whether the Omega_2 targets are in/out
                    closest_points_position,
                    closest_point_as_face_point,
                    target_far_field,
                    target_interpolant_spacing,
                    interpolation_directions,
                    VAR_U,
                    solved_density, 
                    computed_potential);

            // TODO interpolating from source points to non-collocation target
            // points 
            interpolated_density = interpolate_and_resample(
                    kernel.get_tdof(), solved_density,
                    solver->patch_samples(),
                    //solver->patch_samples());
                    samples.get());

            // compute interior value
           VecAXPY(computed_potential, .5, interpolated_density);
           VecDestroy(&interpolated_density);

            break;
        case EvaluationScheme::AUTOEVAL_QBKIX:
            solver->evaluate(targets,
                    solved_density, 
                    computed_potential,
                    closest_points,
                    VAR_U, true);
            break;

        case EvaluationScheme::CHECK_DENSITY:
            assert(Petsc::get_vec_local_size(solved_density) == 
                    Petsc::get_vec_local_size(computed_potential));
            VecCopy(solved_density, computed_potential);
        default:
            break;
    }

}
void run_coarse_fmm(unique_ptr<SolverGMRESDoubleLayer>& solver,
        Vec boundary_data, Vec targets){
    Vec computed_potential;
    Vec t;
    VecDuplicate(boundary_data, &computed_potential);
    VecSet(computed_potential, 0.);
    VecDuplicate(solver->patch_samples()->sample_point_3d_position(),
            &t);
    PetscRandom rctx;
    PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    VecSetRandom(t,rctx);
    PetscRandomDestroy(&rctx);
    
    Kernel3d k(solver->equation_type() + DOUBLE_LAYER + VAR_U, solver->eqcoefs());
    
    unique_ptr<PvFMM> fmm(new PvFMM(
                solver->patch_samples()->sample_point_3d_position(),
                solver->patch_samples()->sample_point_normal(),
                solver->patch_samples()->sample_point_3d_position(), k
                ));
    int64_t n;
    VecGetSize(solver->patch_samples()->sample_point_3d_position(), &n);
    n /= DIM;
    fmm->evaluate(boundary_data, computed_potential);
    double fmm_benchmark_time = omp_get_wtime();
    fmm->evaluate(boundary_data, computed_potential);
    stats.add_result("FMM benchmark time", omp_get_wtime() - fmm_benchmark_time);
    cout << "Reference FMM size: " << n << " sources, " << n << "targets: time " << omp_get_wtime() - fmm_benchmark_time  << endl;
    VecDestroy(&computed_potential);

}

void run_test(PatchSurf* surface,
        TestConfig test_type, bool dump_values){
    Vec true_potential;
    Vec computed_potential;
    Vec targets;
    Vec boundary_data;
    Vec solved_density;

    auto solver = setup_solver(surface, test_type.solver_matvec_type);
    

    NumVec<OnSurfacePoint> closest_points;
    setup_target_data(solver, test_type, targets, closest_points, 
            computed_potential, true_potential);
    setup_problem_data(solver, test_type, boundary_data, solved_density);
    
    //TODO Implement RHS refinement inside pipeline here.
    
    cout << "solution_scheme: " << int(test_type.solution_scheme) << endl;
    solve_and_evaluate(solver,
            boundary_data,
            targets,
            closest_points, solved_density, 
            computed_potential,
            test_type);
    cout << "solved density" << endl;
    VecView(solved_density, PETSC_VIEWER_STDOUT_SELF);
    if(test_type.time_coarse_fmm){
        run_coarse_fmm(solver, boundary_data, targets);
    }
    
    Kernel3d kernel = solver->problem_kernel();
    int target_dof= kernel.get_tdof();
    tear_down(targets, true_potential, computed_potential, test_type, target_dof, 1e-6, 1e-6, dump_values);

    VecDestroy(&true_potential);
    VecDestroy(&computed_potential);
    VecDestroy(&solved_density);
    VecDestroy(&boundary_data);
    VecDestroy(&targets);

}

void test_far_eval_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values){

    Vec true_potential;
    Vec computed_potential;
    Vec targets = nullptr;
    Vec boundary_data;
    Vec solved_density;
    
    
    // 1 in interior for interior problem, 0 in exterior for exterior problem
    double true_potential_value = Options::get_int_from_petsc_opts("-dom") == 0 ? 1. : 0. ;
    setup_target_data_constant_density(solver, true_potential_value, targets, 
            computed_potential, true_potential);
    setup_problem_data_constant_density(solver, boundary_data, solved_density);
    
    Kernel3d kernel = solver->problem_kernel();
    int target_dof= kernel.get_tdof();
    int num_target_points = Petsc::get_vec_local_size(targets)/DIM;
    
    solver->fareval(targets, VAR_U, boundary_data, computed_potential);

    TestConfig test;
    test.error_filename = "laplace_err.vtp";
    tear_down(targets, true_potential, computed_potential, test, target_dof, 1e-6, 1e-6, dump_values);

    VecDestroy(&boundary_data);
    VecDestroy(&solved_density);
    VecDestroy(&targets);
    VecDestroy(&true_potential);
    VecDestroy(&computed_potential);

}

void test_singular_quadrature_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values){

    Vec targets = solver->patch_samples()->sample_point_3d_position();
    Vec targets_as_face_points = solver->patch_samples()->sample_as_face_point();
    Vec true_potential;
    Vec computed_potential;
    Vec boundary_data;
    Vec solved_density;
    
    
    // on-surface values should always be 1/2
    setup_target_data_constant_density(solver, .5, targets, computed_potential, true_potential);
    
    setup_problem_data_constant_density(solver, boundary_data, solved_density);
    
    Kernel3d kernel = solver->problem_kernel();
    int target_dof= kernel.get_tdof();
    int num_target_points = Petsc::get_vec_local_size(targets)/DIM;
    
    //solver->solve(boundary_data, solved_density);
    solver->roneval(targets,targets_as_face_points,
            VAR_U, boundary_data, computed_potential);

    TestConfig test;
    test.error_filename = "laplace_err.vtp";
    tear_down(targets, true_potential, computed_potential, test, target_dof, 1e-5, 1e-5, dump_values);

    VecDestroy(&boundary_data);
    VecDestroy(&solved_density);
    VecDestroy(&true_potential);
    VecDestroy(&computed_potential);

}

void test_solver_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values){

    Vec targets = nullptr;
    Vec true_potential;
    Vec computed_potential;
    Vec boundary_data;
    Vec solved_density;
    
    
    double true_potential_value = Options::get_int_from_petsc_opts("-dom") == 0 ? 1. : 0. ;
    setup_target_data_constant_density(solver, true_potential_value, targets, 
            computed_potential, true_potential);
    
    setup_problem_data_constant_density(solver, boundary_data, solved_density);
    
    Kernel3d kernel = solver->problem_kernel();
    int target_dof= kernel.get_tdof();
    int num_target_points = Petsc::get_vec_local_size(targets)/DIM;
    
    solver->solve(boundary_data, solved_density);
    solver->fareval(targets, VAR_U, solved_density, computed_potential);

    TestConfig test;
    test.error_filename = "laplace_err.vtp";
    tear_down(targets, true_potential, computed_potential, test, target_dof, 1e-5, 1e-5, dump_values);

    VecDestroy(&boundary_data);
    VecDestroy(&solved_density);
    VecDestroy(&true_potential);
    VecDestroy(&computed_potential);

}

void test_qbkix_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values){

    Vec targets;
    Vec temp =solver->patch_samples()->sample_point_3d_position();
    VecDuplicate(temp, &targets);
    VecCopy(temp, targets);
    VecScale(targets, .999);

    
    //Vec targets;
    Vec true_potential;
    Vec computed_potential;
    Vec boundary_data;
    Vec solved_density;
    
    
    // on-surface values should always be 1/2
    double true_potential_value = Options::get_int_from_petsc_opts("-dom") == 0 ? 1. : 0. ;
    setup_target_data_constant_density(solver, true_potential_value, 
            targets, computed_potential, true_potential);
    
    setup_problem_data_constant_density(solver, boundary_data, solved_density);
    
    Kernel3d kernel = solver->problem_kernel();
    int target_dof= kernel.get_tdof();
    int num_target_points = Petsc::get_vec_local_size(targets)/DIM;
    
    NumVec<OnSurfacePoint> on_surface_points = solver->patch_samples()->sample_point_as_on_surface_point();//(num_target_points);
    cout << "on_surface_points: " << on_surface_points.m() << endl;

    solver->evaluate(targets,
            boundary_data,
            computed_potential,
            on_surface_points,
            VAR_U, 
            true);

    TestConfig test;
    test.error_filename = "laplace_err.vtp";
    tear_down(targets, true_potential, computed_potential, test, target_dof, 1e-5, 1e-5, dump_values);

            stats.print_results(); 
    VecDestroy(&boundary_data);
    VecDestroy(&solved_density);
    VecDestroy(&true_potential);
    VecDestroy(&computed_potential);

}






Vec evaluate_function(Vec targets, double (*function)(Point3)){

    int num_targets = Petsc::get_vec_local_size(targets)/DIM;
    
    Vec function_values;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets, function_values);

    DblNumMat function_values_local = get_local_vector(1,  num_targets, function_values);
    DblNumMat targets_local = get_local_vector(DIM,  num_targets, targets);
    
    for(int i =0; i < num_targets; i++){
        Point3 target(targets_local.clmdata(i));
        double function_value = function(target);
        function_values_local(0,i) = function_value;
        
    }
    
    function_values_local.restore_local_vector();
    targets_local.restore_local_vector();
    return function_values;
}



Vec compute_3d_position_of_on_surface_points(
        NumVec<OnSurfacePoint> on_surface_points,
        PatchSurfFaceMap* face_map){
    int num_targets = on_surface_points.m();

    Vec on_surface_points_3d_position;
    Petsc::create_mpi_vec(MPI_COMM_WORLD,
            num_targets*DIM,
            on_surface_points_3d_position);

    DblNumMat on_surface_points_3d_position_local = 
        get_local_vector(DIM, num_targets, on_surface_points_3d_position);
    for(int i =0; i < num_targets; i++){
        OnSurfacePoint p = on_surface_points(i);

        Point3 closest_on_surface_point(on_surface_points_3d_position_local.clmdata(i));
        auto patch = dynamic_cast<FaceMapSubPatch*>(face_map->patches()[p.parent_patch]);

        patch->xy_to_patch_coords(p.parametric_coordinates.array(), 
                PatchSamples::EVAL_VL, closest_on_surface_point.array());

        for(int d = 0; d < DIM; d++){
            on_surface_points_3d_position_local(d,i) = closest_on_surface_point(d);
        }
    }
    
    on_surface_points_3d_position_local.restore_local_vector();
    return on_surface_points_3d_position;
}

void test_qbkix_extrapolation_no_quad(PatchSurf* surface, 
        double (*true_function)(Point3),
        bool dump_values){

    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    int qbkix_order = Options::get_int_from_petsc_opts("-near_interpolation_num_samples");
    
    auto solver = setup_solver(surface,EXTRAPOLATION_AVERAGE,false);
    int num_patches = solver->bdry()->patches().size();
    
    Vec targets = solver->patch_samples()->sample_point_3d_position();
    Vec true_function_values = evaluate_function(targets, true_function);

    Vec computed_function_values;
    VecDuplicate(true_function_values, &computed_function_values);
    VecCopy(true_function_values, computed_function_values);

    int num_targets = Petsc::get_vec_local_size(targets)/DIM;

    vector<int> qbkix_indices;
    for(int i =0; i < qbkix_order; i++)
        qbkix_indices.push_back(i);
    vector<int> patch_indices;
    for(int i =0; i < num_patches; i++)
        patch_indices.push_back(i);

    Vec qbkix_points = solver->patch_samples()->generate_qbkix_points_from_sample_points(
            qbkix_indices, patch_indices);

    Vec qbkix_function_values = evaluate_function(qbkix_points, true_function);

    DblNumMat targets_local = get_local_vector(DIM, num_targets, targets);
    
    int target_dof =1; // scalar function values
    DblNumMat qbkix_function_values_local = 
        get_local_vector(target_dof, num_targets*qbkix_order, qbkix_function_values);

    DblNumMat computed_function_values_local = 
        get_local_vector(target_dof, num_targets, computed_function_values);
    
    NumVec<OnSurfacePoint> on_surface_points = 
        Markgrid::mark_target_points(targets_local, 
                dynamic_cast<PatchSurfFaceMap*>(surface), false);

    Vec on_surface_points_3d_positions =
        compute_3d_position_of_on_surface_points(on_surface_points, 
                dynamic_cast<PatchSurfFaceMap*>(surface));

    DblNumMat on_surface_points_3d_positions_local = 
        get_local_vector(DIM, num_targets, on_surface_points_3d_positions);
    vector<double> interpolation_nodes = get_interpolation_nodes(EXTRAPOLATE_ONE_SIDE);

    DblNumMat expansion_dist_from_boundary_local =
        get_local_vector(1, num_targets, solver->patch_samples()->sample_point_far_field());
    DblNumMat node_spacing_local=
        get_local_vector(1, num_targets, solver->patch_samples()->sample_point_interpolant_spacing());

    lagrange_extrapolation_bary(targets_local,
            on_surface_points_3d_positions_local,
            interpolation_nodes,
            qbkix_function_values_local,
            node_spacing_local,
            expansion_dist_from_boundary_local,
            target_dof, 
            &extrapolation_eval_point_blendsurf,
            computed_function_values_local);
    
    node_spacing_local.restore_local_vector();
    expansion_dist_from_boundary_local.restore_local_vector();
    targets_local.restore_local_vector();
    qbkix_function_values_local.restore_local_vector();
    computed_function_values_local.restore_local_vector();
    on_surface_points_3d_positions_local.restore_local_vector();

    TestConfig test;
    test.error_filename = "laplace_err.vtp";
    tear_down(targets, true_function_values, computed_function_values, test, 
            target_dof, 1e-14, 1e-14, dump_values);
    VecDestroy(&qbkix_points);
    VecDestroy(&qbkix_function_values);
    VecDestroy(&true_function_values);
    VecDestroy(&computed_function_values);
    //VecDestroy(&true_function_values);
    VecDestroy(&on_surface_points_3d_positions);
    
}
double eval_x(Point3 p){
    return p.x();
}
double eval_exp_sin_sin(Point3 p){
    double x = p.x()/(2*M_PI);
    double y = p.y()/(2*M_PI);
    double z = p.z()/(2*M_PI);
    return exp(sin(x)*sin(y));
}

Point2 grad_exp_sin_cos(Point2 s){
    return Point2(
            cos(s.x())*cos(s.y())*exp(sin(s.x())*cos(s.y())),
            sin(s.x())*sin(s.y())*exp(sin(s.x())*cos(s.y()))
            );
}
Point3 grad_exp_sin_sin_sin(Point3 s){
    double x = s.x()/(2*M_PI);
    double y = s.y()/(2*M_PI);
    double z = s.z()/(2*M_PI);
    return Point3(
            exp(sin(x)*sin(y)*sin(z)),
            exp(sin(x)*sin(y)*sin(z)),
            exp(sin(x)*sin(y)*sin(z)));
}

void test_gmres_solve_near_eval(
        PatchSurfFaceMap* surface,
        string output_folder,
        int i){
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type= SingularityType::SPHERE;
    //test.sphere_radius_bc = 2.;
    string d = Test::get_domain();
    if( d== "new_ppp" || d == "comb2" || d == "nearly_touching"){
        test.sphere_radius_bc = 2.;
    } else {
        test.sphere_radius_bc = 1.;
    }
    int is_adaptive = Options::get_int_from_petsc_opts("-adaptive");
    if(is_adaptive){
        if (Test::get_domain() == "interlocking_torii_flip"){
            test.bc_type = BoundaryDataType::CONSTANT_DENSITY;
        }
        
    }

    test.target_type   = TargetType::ON_SURFACE_NON_COLLOCATION;
    // Evaluate solution at points in the far-field via Green's Identity
    test.evaluation_scheme = EvaluationScheme::ON_QBKIX_AVERAGE;
    // no solve step need
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    //test.solver_matvec_type = INTERIOR_EXTRAPOLATION;
    test.solver_matvec_type = EXTRAPOLATION_AVERAGE;
    test.dump_values = true;
    
    string filename;
    int is_adaptive_test = Options::get_int_from_petsc_opts("-adaptive");
    if(is_adaptive_test){
        filename = string("qbx_adaptive_gmres_ref_lvl_")+to_string(i);
    } else {
        filename = string("qbx_gmres_ref_lvl_")+to_string(i);
    }
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;
    run_test(surface, test);

    int num_coarse_patches = surface->patches().size();
    stats.add_result("num coarse patches", num_coarse_patches);
    stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
    stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
    stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
    stats.dump_key_values("", stats._file_prefix);
    stats.print_results(); 
    stats.clear();

}

void test_gmres_solve_far_field_eval(
        PatchSurfFaceMap* surface,
        string output_folder,
        int i){
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type= SingularityType::SPHERE;
    test.sphere_radius_bc = 2.;
    
    // Evaluate solution at points in the far-field
    if(Test::get_domain() == "newtorus"){
        test.target_type   = TargetType::RING;
        test.sphere_radius_targets= .5;
    } else if (Test::get_domain() == "ppp"){
        test.target_type   = TargetType::SPHERE_FAR;
        test.sphere_radius_targets= .1;
    } else if (Test::get_domain() == "ttorus2" || Test::get_domain() == "squished_cube"){
        test.target_type   = TargetType::SPHERE_FAR;
        test.sphere_radius_targets= .01;
    } else {
        test.target_type   = TargetType::SPHERE_FAR;
        test.sphere_radius_targets= .2;
    }
    // ... via Green's Identity
    test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
    // no solve step need
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    test.solver_matvec_type = EXTRAPOLATION_AVERAGE;
    test.dump_values = true;
    
    string filename;
    int is_adaptive_test = Options::get_int_from_petsc_opts("-adaptive");
    if(is_adaptive_test){
        filename = string("qbx_adaptive_gmres_ref_lvl_")+to_string(i);
    } else {
        filename = string("qbx_gmres_ref_lvl_")+to_string(i);
    }
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;
    run_test(surface, test);

    int num_coarse_patches = surface->patches().size();
    stats.add_result("num coarse patches", num_coarse_patches);
    stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
    stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
    stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
    stats.dump_key_values("", stats._file_prefix);
    stats.print_results(); 
    stats.clear();

}

void gmres_test_base_options(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".045");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".045");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".04");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".006666");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".090909");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".090909");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "20");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "10000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    Options::set_value_petsc_opts("-bis3d_spacing", "0.05");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
}

void solver_test_base_options(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");/*
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".06");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".06");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");
    Options::set_value_petsc_opts("-bis3d_spacing", ".090909");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".090909");
    */

    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".135");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".135");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".125");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "20");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "20000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    Options::set_value_petsc_opts("-bis3d_spacing", ".04");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01");
    
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
    // For 7 digits but stagnating
    //Options::set_value_petsc_opts("-bis3d_spacing", ".04");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01333");
   
   
    
    //Options::set_value_petsc_opts("-bis3d_spacing", ".05");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".02");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".02");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");

    // 3/7/19: stagnating at 7 digits with one level upsample: need two with
    // lower quad order
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".0333");
    Options::set_value_petsc_opts("-bis3d_spacing", ".04");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".00666");
    
    
    Options::set_value_petsc_opts("-bis3d_ptsmax", "500");
    if(Options::get_int_from_petsc_opts("-kt") == 111)
    Options::set_value_petsc_opts("-bis3d_np", "20");
    else
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_spacing", "0.05");
    Options::set_value_petsc_opts("-bis3d_spacing", "0.0454545");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.09090909");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "2");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.0588");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".066");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0111");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");

    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".02");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".00333");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".015");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0025");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
}


unique_ptr<PatchSurfFaceMap> setup_face_map(PatchSurfFaceMap::SurfaceType surface_type, 
        int num_uniform_refinement_lvls, string ref_type){
    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type = surface_type;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    cout << "REF TYPE" << ref_type << endl;
    if(ref_type == "uniform"){
        surface->refine_uniform(num_uniform_refinement_lvls);
    } else if(ref_type == "adaptive"){
        int is_adaptive_test = Options::get_int_from_petsc_opts("-adaptive");
        if(is_adaptive_test){
            surface->refine_uniform(num_uniform_refinement_lvls);
        }


        surface->refine();

    }
    return std::move(surface);
}

void setup_and_run_face_map_convergence_test(
        int patch_order,
        int patch_refinement_factor,
        int kernel_enum,
        PatchSurfFaceMap::SurfaceType surface_type,
        string domain,
        string output_folder,
        void (*convergence_test)(PatchSurfFaceMap*, string, int),
        void (*options_init)(),
        string polynomial_patch_filename,
        int num_iterations,
        string ref_type){
    Options::set_value_petsc_opts("-bd3d_filename", string("wrl_files/") + string(domain));
    Options::set_value_petsc_opts("-bd3d_meshfile", string("wrl_files/") + string(domain));
    Options::set_value_petsc_opts("-poly_coeffs_file", 
            string("wrl_files/poly/") + string(polynomial_patch_filename));
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", to_string(patch_order));
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", to_string(patch_refinement_factor));
    Options::set_value_petsc_opts("-kt", to_string(kernel_enum));
    options_init();

    stats._file_prefix = "data/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    for(int i = 0; i < num_iterations; i++){
        Options::set_value_petsc_opts("-debug_test_iter", to_string(i));
        if(num_iterations == 1 && 
                domain == "ttorus2.wrl" &&
                kernel_enum == 311
          ){
            i = 3;
        }
            unique_ptr<PatchSurfFaceMap> surface = setup_face_map(surface_type, i, ref_type);
            /*if(//output_folder == "output/test_greens_identity/" &&
                    domain == "cube.wrl" && kernel_enum == 111 && i > 3
                    ){
                Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", to_string(patch_refinement_factor+1));
                surface = setup_face_map(surface_type, i-1, ref_type);
            } else {*/
                //surface = setup_face_map(surface_type, i, ref_type);
            //}
        convergence_test(surface.get(), output_folder, i);
        //test_gmres_solve_far_field_eval(surface.get(), output_folder, i);
        stats.print_results(); 
        stats.clear();
    }
}

void laplace_singluarity_cube(Vec samples, int dof,Vec& potential){
    DblNumMat samples_local(DIM, samples);
    DblNumMat normals(DIM, samples_local.n());
    Kernel3d laplace(111, vector<double>(2,1.));
    DblNumMat potential_local(laplace.get_tdof(), potential);
    for (int i = 0; i < potential_local.n(); i++) {
        Point3 x(20.,0.,0.);
        Point3 y(samples_local.clmdata(i));
        potential_local(0,i) = 1./pow((x-y).l2(),2);
    }

}
void laplace_singluarity_flat_patch(Vec samples, int dof,Vec& potential){
    DblNumMat samples_local(DIM, samples);
    DblNumMat normals(DIM, samples_local.n());
    Kernel3d laplace(111, vector<double>(2,1.));
    DblNumMat potential_local(laplace.get_tdof(), potential);
    for (int i = 0; i < potential_local.n(); i++) {
        Point3 x(0.,0.,.08);
        Point3 y(samples_local.clmdata(i));
        potential_local(0,i) = 1./pow((x-y).l2(),2);
        
    }
    

}
void laplace_singluarity_propeller(Vec samples, int dof,Vec& potential){
    DblNumMat samples_local(DIM, samples);
    DblNumMat normals(DIM, samples_local.n());
    Kernel3d laplace(111, vector<double>(2,1.));
    DblNumMat potential_local(laplace.get_tdof(), potential);
#pragma omp parallel for
    for (int i = 0; i < potential_local.n(); i++) {
        Point3 x(0., 0., 1.);
        Point3 y(samples_local.clmdata(i));
        potential_local(0,i) = 1./pow((x-y).l2(),2);
        
    }
}



END_EBI_NAMESPACE
