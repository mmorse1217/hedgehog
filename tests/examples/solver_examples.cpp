
#include "../catch.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/utils.hpp"
#include "common/vtk_writer.hpp"
#include <memory>
using namespace hedgehog;
TEST_CASE("Solver examples: blended surface", "[examples][solver][blended]"){
    // default options are listed in opt/morse_cases.opt. Any options explcitly
    // set here have priority and override defaults
    
    //1. Setup surface and solver
    
    // Setup surface
    // Relevant options in opt/morse_cases.opt are prefixed with "bd3d_". 
    // Best to stick with the defaults there.
    // important ones: 
    //  -bd3d_bdsurf_chttyp: type of blended patches to fit. These types
    //      are described in Ying/Zorin 2004 and Tosun/Zorin 2011. Isodistance
    //      gives best BIE results but all are feasible to solve on.
    //  -bd3d_bdsurf_ctrllvl: number of subdivision levels to compute before
    //      fitting the blended surface to the resulting subdivided mesh.
    //  -bd3d_bdsurf_pouctrl: type of partition of unity function used to 
    //      represent the blended surface
    // -bd3d_bdsurf_submatlibfile: needs to point to location of ccsubmatall.dat in blendsurf lib
    // -bd3d_bdsurf_bdulibfile: needs to point to location of bdsurf_U_ONE.dat in blendsurf lib
    // -bd3d_filename and -bd3d_meshfile: the geometry file to load. both should
    //      be the same value
    
    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();
    
    // Set up the solver
    // relevant options:
    // bis3d_spacing: for blended/analytic surfaces: corresponds to spacing of 
    //      sample points in 3D space
    // -boundary_distance_ratio: for hedgehog on blended/analytic surfaces:
    //      for the extrapolation distance a*h of the singular quadrature,
    //      corresponds to parameter a
    //      i.e. propotional to distance from surface to first check point (h is bis3d_spacing)
    // -interpolation_spacing_ratio: for hedgehog on blended/analytic surfaces:
    //      for distance between consecutive check points b*h of the singular quadrature,
    //      corresponds to parameter b
    //      i.e. propotional to distance two check points (h is bis3d_spacing)
    vector<int> patch_partition(surface->patches().size(), 0);  //All in one processor
    unique_ptr<SolverGMRESDoubleLayer> solver( new SolverGMRESDoubleLayer("BIS3D_", "bis3d_"));
    solver->bdry() = surface.get();
    solver->patch_partition() = patch_partition;
    solver->dom() = Options::get_int_from_petsc_opts("-dom"); // TODO remove

    // specify which matvec to use inside GMRES
    // SINGULAR_EVAL = the partition of unity scheme in Ying/Biros/Zorin 2006
    // INTERIOR_EXTRAPOLATION= one-sided hedgehog (evaluates interior limit on surface)
    // EXTRAPOLATION_AVERAGE = two-sided hedgehog (evalutes interior and exterior limits on surafce and averages; this is preferred approach)
    // SINGULAR_EVAL only works on blended surfaces
    solver->set_evaluation_type(SINGULAR_EVAL); 
    solver->_compute_refined_surface = false;

    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4; // choose non-trivial coefficient if running elasticity problem
    
    solver->setFromOptions();
    solver->setup(); 


    // 2. Setup of boundary condition:
    // 2a. initialize petsc vectors
    
    Vec boundary_condition;
    Vec solved_density;
    Vec point_charge_positions;
    Vec point_charge_densities;
    const int num_point_charges = 3;
    const int num_sample_points = solver->patch_samples()->local_num_sample_points();
    Kernel3d kernel = solver->problem_kernel();
    // number of degrees of freedom per sample for the... 
    const int source_dof = kernel.get_sdof(); // ... density... 
    const int target_dof = kernel.get_tdof(); // ... and the potential...
    // ... for the problem of interest.

    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_point_charges * DIM, point_charge_positions);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_point_charges * source_dof, point_charge_densities);
    // Boundary condition and solution needs to match the size of the solver discretization,
    // so we copy its size here
    Petsc::create_mpi_vec(MPI_COMM_WORLD,num_sample_points*source_dof,boundary_condition);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sample_points*source_dof, solved_density);

    // 2b.chose random points outside the domain with random charges
    {
      DblNumMat point_charge_positions_local(
          DIM, point_charge_positions); // get local writable petsc vector
                                        // (first arg is stride size)
      setvalue(point_charge_positions_local, 0.);
      // Choose axis aligned point charges for simplicity
      point_charge_positions_local(0, 0) = 1.;  // (1,0,0)
      point_charge_positions_local(1, 1) = -1.; //  (0,-1,0)
      point_charge_positions_local(2, 2) = 1.;  // (0,0,1)
      
      DblNumMat point_charge_densities_local( source_dof, point_charge_densities); 
      setvalue(point_charge_densities_local, 4*M_PI);
                                       
    } // destructor restores local vector
    
    // 2c. evaluate potential  due to point charges at sample points to create
    // boundary condition using an FMM
    unique_ptr<PvFMM> fmm(new PvFMM(
        point_charge_positions, // source positions
        NULL, // normal vectors unused for laplace single layer
        solver->patch_samples()->sample_point_3d_position(), // target positions
        kernel)); // kernel to evaluate
    fmm->evaluate(point_charge_densities, boundary_condition);

    // 3. Solve the PDE
    solver->solve(boundary_condition, solved_density);
}
TEST_CASE("Solver examples: face-map surface", "[examples][face-map][solver]"){
    // 0. Set options/configuration for solve
    // default options are listed in opt/morse_cases.opt. Any options explcitly
    // set here have priority and override defaults
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem (311=Stokes, 511=Elasticity)
    Options::set_value_petsc_opts("-dom", "0"); // 0 =interior problem , 1=exterior
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl"); // sphere-like blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
    // Quadrature related
    Options::set_value_petsc_opts("-bis3d_spacing", ".090909"); // floor(1/-bis3d_spacing) + 1 quadrature order
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6"); // hedgehog order
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");  // two levels
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".06");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01");
    
    // Surface related 
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "10"); // polynomial surface order
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive","0"); // non-adaptive , don't turn on for now
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","1");// one level polynomial splitting 
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic"); //changes hedgehog point spacing from O(h) to O(\sqrt(h))

    // pvfmm related 
    Options::set_value_petsc_opts("-bis3d_ptsmax", "1000"); // points per leaf box in pvfmm 
    Options::set_value_petsc_opts("-bis3d_np", "10");  // multipole order
    Options::set_value_petsc_opts("-pou_radius_constant", "1.1");
    
    //1. Setup surface and solver
    
    // Setup surface
    // Relevant options in opt/morse_cases.opt are prefixed with "bd3d_". 
    // Best to stick with the defaults there.
    // important ones for face-map (polynomial) surfaces: 
    // -bd3d_facemap_patch_order : polynomial/spline order for poly/spline surfaces
    // -bd3d_facemap_refinement_factor: amount of uniform refinement of the
    //      underlying surface before fitting patches to form polynomial surface
    // -bd3d_facemap_fit_accuracy: accuracy of adaptive polynomial patch fitting
    // -bd3d_filename and -bd3d_meshfile: the geometry file to load. both should
    //      be the same value
    // -poly_coeffs_file: file containing polynomial patch coefficients
    // Note! All  blendsurf parameters above are relevant for polynomials
    // surfaces. Polynomial surfaces are constructed by fitting polynomnials to
    // blended surfaces, so blendsurf parameters still control the underlying
    // surface.
    // In particular, make sure to use bd3d_bdsurf_chttyp  == isodistance and 
    // -bd3d_pouctrl == spline

    unique_ptr<PatchSurfFaceMap> surface( new PatchSurfFaceMap("BD3D_", "bd3d_"));
    // Indicate the underlying surface we are approximating with polynomials.
    // BLENDED = fit polynomials to a blended surface
    // POLYNOMIAL = explicit polynomials are given in a file in the
    // -poly_coeffs_file option
    // others are for debugging, just use these two
    surface->_surface_type = PatchSurfFaceMap::BLENDED;

    // Flag this surface as the coase patch discretization for the solver
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine_test(); // necessary for refinement compatibility

    // Set up the solver
    // relevant options (Note that their meanings are slightly different):
    // -bis3d_spacing: for polynomial/spline surface: corresponds to spacing of
    //      sample points in PARAMETER space. To increase number of samples for
    //      these surfaces, decrease this number: \floor(1/bis3d_spacing) + 1
    //      corresponds to the chebyshev discretization order.
    // -boundary_distance_ratio: for hedgehog on polynomial/spline surfaces:
    //      for the extrapolation distance a*L of the singular quadrature,
    //      corresponds to parameter a
    //      i.e. propotional to distance from surface to first check point (L is patch length)
    // -interpolation_spacing_ratio: for hedgehog on polynomial/spline surfaces:
    //      for distance between consecutive check points b*h of the singular quadrature,
    //      corresponds to parameter b
    //      i.e. propotional to distance two check points (h is bis3d_spacing)
    vector<int> patch_partition(surface->patches().size(), 0);  //All in one processor
    unique_ptr<SolverGMRESDoubleLayer> solver( new SolverGMRESDoubleLayer("BIS3D_", "bis3d_"));
    solver->bdry() = surface.get();
    solver->patch_partition() = patch_partition;
    solver->dom() = Options::get_int_from_petsc_opts("-dom"); // TODO remove

    // specify which matvec to use inside GMRES
    solver->set_evaluation_type(EXTRAPOLATION_AVERAGE); 
    solver->_compute_refined_surface = true;

    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4; // choose non-trivial coefficient if running elasticity problem
    
    solver->setFromOptions();
    solver->setup(); 


    // 2. Setup of boundary condition:
    // 2a. initialize petsc vectors
    
    Vec boundary_condition;
    Vec solved_density;
    Vec solution;
    Vec point_charge_positions;
    Vec point_charge_densities;
    const int num_point_charges = 3;
    const int num_sample_points = solver->patch_samples()->local_num_sample_points();
    Kernel3d kernel = solver->problem_kernel();
    // number of degrees of freedom per sample for the... 
    const int source_dof = kernel.get_sdof(); // ... density... 
    const int target_dof = kernel.get_tdof(); // ... and the potential...
    // ... for the problem of interest.

    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_point_charges * DIM, point_charge_positions);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_point_charges * source_dof, point_charge_densities);
    // Boundary condition and solution needs to match the size of the solver discretization,
    // so we copy its size here
    Petsc::create_mpi_vec(MPI_COMM_WORLD,num_sample_points*source_dof,boundary_condition);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sample_points*source_dof, solved_density);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sample_points*target_dof, solution);

    // 2b.chose random points outside the domain with random charges
    {
      DblNumMat point_charge_positions_local(
          DIM, point_charge_positions); // get local writable petsc vector
                                        // (first arg is stride size)
      setvalue(point_charge_positions_local, 0.);
      // Choose axis aligned point charges for simplicity
      point_charge_positions_local(0, 0) = 1.;  // (1,0,0)
      point_charge_positions_local(1, 1) = -1.; //  (0,-1,0)
      point_charge_positions_local(2, 2) = 1.;  // (0,0,1)
      
      DblNumMat point_charge_densities_local( source_dof, point_charge_densities); 
      setvalue(point_charge_densities_local, 4*M_PI);
                                       
    } // destructor restores local vector
    
    // 2c. evaluate potential  due to point charges at sample points to create
    // boundary condition using an FMM
    unique_ptr<PvFMM> fmm(new PvFMM(
        point_charge_positions, // source positions
        point_charge_positions, // normal vectors unused for laplace single layer, just pass a dummy vector
        solver->patch_samples()->sample_point_3d_position(), // target positions
        kernel)); // kernel to evaluate
    fmm->evaluate(point_charge_densities, boundary_condition);

    // 3. Solve the PDE
    solver->solve(boundary_condition, solved_density);
    
    // 4. Evaluate solution at desired target points
    // This is automatic (in theory). pass a set of 3D points at it will
    // evaluate them all accurately with smooth quadrature or near-singular
    // quadrature appropriately
    
    // Optional to speed up computation. each target point needs to be "marked"
    // as near or far from the boundary to use near-singular quadrature or
    // smooth quadrature to compute the layer potential. 
    // solver->evaluate() computes this for you: you need only pass 
    //   NumVec<OnSurfacePoint> on_surface_points(num_target_points);
    // We're evaluating at points sampled inside the sovler, so we can use
    // values that are precomputed. Uncomment the following line to skip point
    // marking.
    NumVec<OnSurfacePoint> on_surface_points = solver->patch_samples()->sample_point_as_on_surface_point();
    //NumVec<OnSurfacePoint> on_surface_points(num_sample_points);
    solver->evaluate(solver->patch_samples()->sample_point_3d_position(),
            solved_density,
            solution,
            on_surface_points,
            VAR_U, 
            true);
    // 5. Compute relative error in the solution
    Vec error;
    VecDuplicate(boundary_condition, &error);
    VecCopy(boundary_condition, error);
    int minus_one = -1.;
    VecAXPY( error, minus_one,  solution); // u_true - u_approx
    VecAbs(error);

    // no precision lower than 1e-16 for paraview
    DblNumMat error_local = get_local_vector(1, Petsc::get_vec_local_size(error), error);
    for(int i = 0; i < Petsc::get_vec_local_size(error); i++){
        error_local(0,i) = error_local(0,i) < 1e-16 ? 1e-16 : error_local(0,i);
    }
    error_local.restore_local_vector(); 

    // 5b. Compute max magnitude of solution 
    Vec solution_absolute_value;
    VecDuplicate(solution, &solution_absolute_value);
    VecCopy(solution, solution_absolute_value);
    VecAbs(solution_absolute_value);
    double max_potential = 1e-16;
    VecMax(solution_absolute_value, NULL, &max_potential);
    // 5c. Scale error by max magnitude
    VecScale(error, 1./max_potential);
    
    VecMax(error, NULL, &max_potential);
    cout << "max relative error " <<  max_potential << endl;

    // Scale
    // 6. Write data to vtk files for paraview. Sorry in advance for the insane
    // function signatures...
    
    // 6a. write surface geometry
    write_face_map_patches_to_vtk(DblNumMat(0,0), 
        vector<int>(surface->num_patches(), 0),
        surface.get(), 0, "surface_");
    
    // 6b. write solved density to a file
    write_general_points_to_vtk(
            solver->patch_samples()->sample_point_3d_position(), 
            source_dof, "density.vtp", solved_density);

    // 6b. write solution to a file
    write_general_points_to_vtk(
            solver->patch_samples()->sample_point_3d_position(), 
            target_dof, "solution.vtp", solution);

    // 6c. write error to a file
    write_general_points_to_vtk(
            solver->patch_samples()->sample_point_3d_position(), 
            1, "relative_error.vtp", error);
}
