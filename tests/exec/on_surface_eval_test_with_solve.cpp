//#include "../bdry3d/bd3dag.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bie3d/evaluator_far.hpp"
#include "bie3d/evaluator_near.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "common/nummat.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/exsol3d.hpp"
#include "common/stats.hpp"
#include "bie3d/solver_utils.hpp"
#include "petscsys.h"
#include "petscviewer.h"
#include "petscsnes.h"
#include <getopt.h>
#include <string>
#include <iomanip> //for setprecision in file streams
#include <time.h>
#include "common/utils.hpp"


// To reinitialize the stored values for regression testing, set this to true
using namespace std;
using namespace hedgehog;

int main(int argc, char** argv)
{
    ebiFunctionBegin;
    cout.precision(16);
    // Initialization steps
    ArgParse::read_command_line_args(argc, argv, "opt/morse_test_generic.opt");
/*
    // Figure out what kernel and domain we're testing
    int is_bounded; // bounded vs. unbounded domain
    int is_blended; // blended/analytic surface representation
    // options file with no kernel or domain data
    const struct option longopts[] = 
    {
        { "kernel",     required_argument,      0, 'k'},
        { "layer",      required_argument,      0, 'l'},
        { "prim-var",   required_argument,      0, 'v'},
        { "domain",     required_argument,      0, 'd'},
        { "et",         required_argument,      0, 'e'},
        { "boundary-distance", required_argument, 0, 'b'},
        { "interpolant-spacing", required_argument, 0, 'i'},
        { "bounded",    no_argument,    &is_bounded, 0},
        { "unbounded",  no_argument,    &is_bounded, 1},
        { "analytic",   no_argument,    &is_blended, 0},
        { "blended",    no_argument,    &is_blended, 1},
        { "face-map",   no_argument,    &is_blended, 2},
        {0,0,0,0},
    };
    int iarg = 0;
    int index;

    int __internal_kernel = 1; // always solving for velocity
    string __domain = "wrl_files/"; //domain to solve problem 
    string __temp_dom;
    int __equation_type = 11; // default to constant solution in Exsol3d
    //equation type; string for now (should be int but this should be fine)
    
    string boundary_distance_ratio_str; 
    string interpolation_spacing_ratio_str;
    while (iarg != -1){
        iarg = getopt_long(argc, argv, "k:", longopts, &index);
        switch (iarg){
            case 'k': 
                cout<< "kernel: " << optarg <<  endl;
                __internal_kernel += ArgParse::parse_kernel(optarg);
                break;
            case 'l': 
                cout<< "layer: " << optarg <<  endl;
                __internal_kernel += ArgParse::parse_layer(optarg);
                break;
            case 'v': 
                cout<< "prim-var: " << optarg <<  endl; //ignore
                break;
            case 'd': 
                cout<< "domain: " << optarg <<  endl;
                __temp_dom += string(optarg);
                __domain += __temp_dom;
                __domain += ".wrl";
                break;
            case 'e':
                __equation_type = atoi(optarg);
                cout << "equation type:" << optarg << endl;
                break;
            case 'b':
                boundary_distance_ratio_str = string(optarg);
                break;
            case 'i':
                interpolation_spacing_ratio_str = string(optarg);
                break;
        }
    }
    // Load & parse CLI arguments into Petsc
    string options_file = "test_optfiles/morse_test_generic.opt";
    PetscInitialize( &argc, &argv, options_file.c_str(), "Test"); eC;

    cout << "is_bounded:" << is_bounded << endl;
    stringstream ss;
    ss << __internal_kernel;
    string s = ss.str();
    char const* __internal_kernel_str = s.c_str(); 
    cout << __internal_kernel_str << endl;
    PetscOptionsSetValue(NULL, "-kt", __internal_kernel_str);

    ss.str(string()); ss.clear();
    ss << is_bounded;
    char const* is_bounded_str = ss.str().c_str(); 
    PetscOptionsSetValue(NULL, "-dom",is_bounded_str);

    ss.str(string()); ss.clear();
    ss << __equation_type;
    char const* equation_type_str = ss.str().c_str(); 
    PetscOptionsSetValue(NULL, "-et", equation_type_str);


    ss.str(string()); ss.clear();
    ss << is_blended;
    char const* surface_rep_str = ss.str().c_str(); 
    PetscOptionsSetValue(NULL, "-bdtype", surface_rep_str);


    PetscOptionsSetValue(NULL, "-bd3d_filename", __domain.c_str());
    PetscOptionsSetValue(NULL, "-bd3d_meshfile", __domain.c_str());
    if(!boundary_distance_ratio_str.empty()){
        PetscOptionsSetValue(NULL, "-boundary_distance_ratio", boundary_distance_ratio_str.c_str()); 
    }
    if(!interpolation_spacing_ratio_str.empty()){
    PetscOptionsSetValue(NULL, "-interpolation_spacing_ratio", interpolation_spacing_ratio_str.c_str()); 
    }
    */
    
    PetscBool flg = PETSC_FALSE;
    double boundary_distance_ratio;
    double interpolation_spacing_ratio;
    PetscOptionsGetReal(NULL, "", "-boundary_distance_ratio", &boundary_distance_ratio, &flg); 
    PetscOptionsGetReal(NULL, "", "-interpolation_spacing_ratio", &interpolation_spacing_ratio, &flg); 
    cout << "boundary_distance_ratio: " << boundary_distance_ratio << ",  interpolation_spacing_ratio: " << interpolation_spacing_ratio << endl;

    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    ebiInit(); eC;

    srand48( (long)time(NULL) );
    MPI_Comm comm; comm = PETSC_COMM_WORLD;
    int mpiRank;   iC( MPI_Comm_rank(comm, &mpiRank) );
    int mpiSize;   iC( MPI_Comm_size(comm, &mpiSize) );  //ebiAssert( mpiSize==1 );

    //char pbname[100];  iC( PetscOptionsGetString(NULL, "", "-pbname", pbname, 100, &flg) ); ebiAssert(flg==PETSC_TRUE);
    //double hs; iC( PetscOptionsGetReal(NULL, "",  "-bis3d_spacing", &hs,  &flg) );  ebiAssert(flg==PETSC_TRUE);
    //double sspc; iC( bis3dspacing2bdsurfspacing(hs, sspc) );

    // 1.) Initializing the patch-based surface representation based on options 
    // file chosen from above

    PatchSurf* bd = NULL;
    int64_t bdtype;

    iC( PetscOptionsGetInt(NULL, "",  "-bdtype", &bdtype,  &flg) );  ebiAssert(flg==PETSC_TRUE);
    // bdtype = 0: analytic (bd3dag.cpp)
    // bdtype = 1: blended (bd3dbd.cpp)

    if(bdtype==0) {	 
        PatchSurfAnalytic* bdag = new PatchSurfAnalytic("BD3D_", "bd3d_");
        iC( bdag->setFromOptions() );	 
        bd = bdag;
    } else if (bdtype == 1){	 
        PatchSurfBlended* bdbd = new PatchSurfBlended("BD3D_", "bd3d_");
        iC( bdbd->setFromOptions() );
        bd = bdbd;
    } else if (bdtype == 2){

        PatchSurfFaceMap* bdfm = new PatchSurfFaceMap("BD3D_", "bd3d_");
        bdfm->_surface_type = PatchSurfFaceMap::BLENDED;
        iC( bdfm->setFromOptions() );	 
        bd = bdfm;
    } else{ assert(0); }

    // but there is another input option of -bd3d_spacing?
    cout << "getting blended surface " << endl;
    bd->setup(); 
    if(bdtype == 2){
    ((PatchSurfFaceMap*) bd)->refine_test();
    //((PatchSurfFaceMap*) bd)->refine();
    }
    cout << "got blended surface " << endl;
    // bd->patches() is the list of overlapping surface patches 
    int numpch = bd->patches().size();
    cout << "number of patches: " << numpch << endl;
    vector<int> patch_partition(numpch, 0);  //All in one processor

    //2.) Initialize the boundary integral solver based on options file chosen above
    // -dom 0: unbounded domain
    // -dom 1: bounded domain
    int64_t dom;  iC( PetscOptionsGetInt(NULL, "",  "-dom", &dom,  &flg) ); ebiAssert(flg==PETSC_TRUE);
    //int64_t dnref;  iC( PetscOptionsGetInt(NULL, "",  "-dnref", &dnref,  &flg) ); ebiAssert(flg==PETSC_TRUE);

    // kernel type for current problem; possible values are in ebi/common/kernel3d.hpp
    int64_t eqn;  iC( PetscOptionsGetInt(NULL, "",  "-kt", &eqn,  &flg) ); 
    ebiAssert(flg==PETSC_TRUE); 

    // LL := number of Lagrange points used to interpolate
    int64_t LL;  iC( PetscOptionsGetInt(NULL, "",  "-LL", &LL,  &flg) ); ebiAssert(flg==PETSC_TRUE);


    SolverGMRESDoubleLayer* bis = new SolverGMRESDoubleLayer("BIS3D_", "bis3d_");
    bis->bdry() = bd;
    bis->patch_partition() = patch_partition;
    bis->dom() = dom;
    bis->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    bis->_compute_refined_surface = false;
    //bis->set_evaluation_type(INTERPOLATION_ACROSS_SURFACE);
    //bis->set_evaluation_type(INTERIOR_EXTRAPOLATION);
    bis->eqcoefs() = vector<double>(2,1.0); // 1 coef for stokes/laplace, 2 coefs for helmholtz and navier
    if (eqn == 511 || eqn == 521){
        //bis->eqcoefs()[0] = 1.0;
        //bis->eqcoefs()[1] = .5;
        bis->eqcoefs()[0] = .7;
        bis->eqcoefs()[1] = .4;
    }
    if(eqn == 211){
        bis->eqcoefs()[1] = .5;
    }
    // Load options into BIS
    iC( bis->setFromOptions() );
    cout << "bis setup" << endl;
    iC( bis->setup() );
    cout << "bis setupdone" << endl;

    //3.) Actually solve the equation
    // -et: equation type; possible values are listed in ebi/fmm3d/exsol3d.cpp/hpp
    int64_t bnd;  iC( PetscOptionsGetInt(NULL, "",  "-et", &bnd,  &flg) ); ebiAssert(flg==PETSC_TRUE);

    Exsol3d exsol(eqn, bis->eqcoefs(), bnd); 

    int la,lb,lm; 
    int ga,gb,gm; 
    iC( bis->localSize( la,lb,lm) );
    cout << la << ", " << lb << ", " << lm << endl;
    iC( bis->globalSize( ga,gb,gm) );
    
    Kernel3d problem_kernel(eqn, bis->eqcoefs());
    int target_dof = problem_kernel.get_tdof();
    int source_dof = problem_kernel.get_sdof();
    // Initialize source points: 6 points at +-1*e_i for i=1,2,3 with e_i unit
    // vectors
    
    
    Vec den_solved; iC( VecCreateMPI(comm, lm, PETSC_DETERMINE, &den_solved) );
    Vec val = evaluate_singularities_along_basis( comm,
            problem_kernel, 
            bis->patch_samples()->sample_point_3d_position());

    cout << "\n\ncalling solve..." << endl;
    stats.start_timer("Solve");
    bis->solve(val, den_solved);
    stats.stop_timer("Solve");
    cout << "called solve...\n\n" << endl;
    // Save density
    {
      PetscViewer viewer;
      PetscViewerBinaryOpen(comm, "tests/stokes_solved_density_singularities.bin", FILE_MODE_WRITE, &viewer);
      VecView(den_solved, viewer);
      PetscViewerDestroy(&viewer);
      }

    {        
        
        
        //--------------------------------------------------------------------

        
        // Compute ||solved-density - known-density||
        //double norm = 0.0;
        //PetscScalar scale = -1.0;
        Vec difference;
        iC(VecCreateMPI(comm, lm, PETSC_DETERMINE, &difference));

        /*
        iC(VecCopy(den_solved, difference));
        iC(VecAXPY( difference, scale, den_known)); //computes y = alpha*x + y (arguments in that order)

        iC(VecNorm(difference, NORM_2, &norm));
        cout << "|| den_solved - den_known||_2= " << norm << endl;
        iC(VecNorm(difference, NORM_INFINITY, &norm));
        cout << "|| den_solved - den_known||_inf= " << norm << endl;
        //assert(norm <=4.);
        iC(VecNorm(den_solved, NORM_2, &norm));
        cout << "Norm of den_solved:  " << norm << endl;
        iC(VecNorm(den_known, NORM_2, &norm));
        cout << "Norm of den known:  " << norm << endl;
        cout << "\n\n\n" << endl;
        */
        //--------------------------------------------------------------------
        // Choose 50 points on a circle in the far-field of the  domain 
        // Ring of points around the origin, lying in the (x,y) plane of
        // varying size
        
        int num_target_points = 50;
        
        Vec target_points_petsc;
        iC(VecCreateMPI(comm, 3*num_target_points, PETSC_DETERMINE, 
                    &target_points_petsc));

        double* target_points_petsc_ptr;
        VecGetArray(target_points_petsc, &target_points_petsc_ptr);
        DblNumMat target_points(3, num_target_points, false, target_points_petsc_ptr);

        Vec val;
        double* valarr;
        iC( VecCreateMPI(comm, target_dof*num_target_points,
                    PETSC_DETERMINE, &val) );
        iC( VecGetArray(val, &valarr) );

        // Choose radii sufficiently far from all domains in wrl_files/
        double r;
        double h = 2*M_PI/(num_target_points -1);
        if (dom  == 0){
            // Bounded domain
            r = 1e-6;
        } else if (dom == 1){
            //Unbounded domain
            r = 2;
        } else {
            cerr << "Domain not defined" << endl;
            ebiAssert(false);
        }

        for (int i = 0; i < num_target_points; i++){
            target_points(0, i) = r*cos(i*h);
            target_points(1, i) = r*sin(i*h);
            target_points(2, i) = 0; 
        }

        
        iC( VecRestoreArray(val, &valarr) );
        iC(VecRestoreArray(target_points_petsc, &target_points_petsc_ptr));
        VecDestroy(&val);
        val = evaluate_singularities_along_basis(comm,
            problem_kernel, 
            target_points_petsc);


        
        // Store the resulting potential at target points
        Vec target_value;
        iC(VecCreateMPI(comm, num_target_points*target_dof, PETSC_DETERMINE, &target_value));
       
        //iC(VecCreateMPI(comm, num_target_points, PETSC_DETERMINE, &target_value));
        int qt;
        qt = VAR_U;
        
        clock_t eval_time = clock(); 
        stats.start_timer("Far Evaluation");
        bis->fareval(target_points_petsc, qt, den_solved, target_value);
        stats.stop_timer("Far Evaluation");


        eval_time = clock() - eval_time;
        cout << endl << "Far evaluation time: " << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 
        cout << "in main()" << endl;
        stats.store_relative_error(val, 
                target_value,
                NORM_INFINITY,
                "Smooth quadrature (far evaluation)");
        Test::compute_relative_error(val, 
                target_value, 
                comm, 
                num_target_points, 
                target_dof,
                "Far Evaluation");

        // Compute || solved_potential - known_potential||_2
        //double exact_solution_norm = 0.0;
        Vec potential_difference;
        iC(VecCreateMPI(comm, target_dof*num_target_points, PETSC_DETERMINE,
                         &potential_difference));
        
        VecDestroy(&potential_difference);
        VecDestroy(&target_value);
        VecDestroy(&target_points_petsc);
        // Pressure evaluation
        qt = VAR_U;
        //--------------------------------------------------------------------
        // Test near-zone evaluator
        //--------------------------------------------------------------------
        Vec near_surface_samples;
        // MJM TODO: remove explicit references to three; possibly replace with
        // bis->dim() everywhere 
        int64_t patch_order;
        PetscOptionsGetInt(NULL,"", "-bd3d_facemap_patch_order", &patch_order, &flg);
        assert(flg);

        iC(VecCreateMPI(comm, bis->patch_samples()->local_num_sample_points()*3,
                    PETSC_DETERMINE, &near_surface_samples));

        VecCopy(bis->patch_samples()->sample_point_3d_position(), near_surface_samples);
        PetscScalar near_zone_scale = (1-pow(1e-2,patch_order-1));
        //PetscScalar intermed_zone_scale = .85; 
        VecScale( near_surface_samples, near_zone_scale);
    
        PatchSamples* patch_samples = bis->patch_samples();
        Vec near_potential;
        iC(VecCreateMPI(comm, target_dof*num_target_points, 
                    PETSC_DETERMINE, &near_potential));

        Vec near_target_in_out;
        iC(VecCreateMPI(comm, bis->patch_samples()->local_num_sample_points(),
                    PETSC_DETERMINE, &near_target_in_out));
        PetscScalar one = 1.0;
        VecSet( near_target_in_out, one);

        Vec near_targets_as_face_points;
        VecDuplicate(patch_samples->sample_as_face_point(), &near_targets_as_face_points);
        VecCopy(patch_samples->sample_as_face_point(), near_targets_as_face_points);
        
        Vec on_surface_samples;
        VecDuplicate(patch_samples->sample_point_3d_position(), &on_surface_samples);
        VecCopy(patch_samples->sample_point_3d_position(), on_surface_samples);

        Vec on_surface_face_points;
        VecDuplicate(patch_samples->sample_as_face_point(), &on_surface_face_points);
        VecCopy(patch_samples->sample_as_face_point(), on_surface_face_points);

        Vec near_target_value;
        iC(VecCreateMPI(comm, patch_samples->local_num_sample_points()*target_dof, 
                    PETSC_DETERMINE, &near_target_value));
        Vec intermediate_targets;
        VecCreateMPI(comm, 0, PETSC_DETERMINE, &intermediate_targets);
        
        Vec near_target_value_extrapolate;
        VecDuplicate(near_target_value, &near_target_value_extrapolate);
        // Perform near evaluation


        eval_time = clock(); 
        stats.start_timer("QBKIX Near");
        bis->neaeval_extrapolate(near_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     qt, 
                     den_solved, 
                     near_target_value_extrapolate);
        stats.stop_timer("QBKIX Near");
        eval_time = clock() - eval_time;
        cout << endl << "Near interior extrapolation evaluation time: " << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 

        eval_time = clock(); 
        if(bdtype != 2){/*
        bis->neaeval(near_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     qt, den_solved, near_target_value);*/
        }
        eval_time = clock() - eval_time;
        cout << endl << "Near interpolation + singular quad. evaluation time: " << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 


        // Compute exact potential
        double* near_potential_array;
        double* near_surface_samples_array;

        Vec exact_near_potential;
        iC(VecCreateMPI(comm, patch_samples->local_num_sample_points()*target_dof, 
                    PETSC_DETERMINE, &exact_near_potential));
        VecGetArray(exact_near_potential, &near_potential_array);
        DblNumVec near_potential_mat(patch_samples->local_num_sample_points()*target_dof,
                false, near_potential_array);

        VecGetArray(near_surface_samples, &near_surface_samples_array);
        DblNumMat near_surface_samples_mat(target_dof, patch_samples->local_num_sample_points(), false, 
                                        near_surface_samples_array);
        
        iC( VecRestoreArray(exact_near_potential, &near_potential_array) );
        iC( VecRestoreArray(near_surface_samples, &near_surface_samples_array) );
        VecDestroy(&val);
        exact_near_potential = evaluate_singularities_along_basis(comm,
            problem_kernel, 
            near_surface_samples);
        
        stats.store_relative_error(exact_near_potential, 
                near_target_value_extrapolate,
                NORM_INFINITY,
                "QBKIX Near evaluation");

        num_target_points = patch_samples->local_num_sample_points();
        Test::compute_relative_error(
                exact_near_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "Near Zone Evaluation via Interior Extrapolation");
        //Compare normwise

        //--------------------------------------------------------------------
        // Test on-surface evaluator
        //--------------------------------------------------------------------
        int num_local_targets = bis->patch_samples()->local_num_sample_points();
        VecDestroy(&near_target_value_extrapolate);
        VecCreateMPI(comm,
                num_local_targets*target_dof,
                PETSC_DETERMINE,
                &near_target_value_extrapolate);

        // Test 3 "qbkix"-like schemes on the same grid used to solve for
        // density
        Vec exact_on_surface_potential;
        iC(VecCreateMPI(comm, 
                    num_local_targets*target_dof, 
                    PETSC_DETERMINE, &exact_on_surface_potential));
        VecDestroy(&exact_near_potential);

        exact_on_surface_potential = evaluate_singularities_along_basis(comm,
            problem_kernel, 
            on_surface_samples);
        
        Vec on_surface_target_value;
        iC(VecCreateMPI(comm, num_local_targets*target_dof,
                    PETSC_DETERMINE, &on_surface_target_value));

        eval_time = clock(); 

        if(bdtype != 2){
        bis->roneval(on_surface_samples,
                     on_surface_face_points,
                     qt, den_solved, on_surface_target_value);
        } 
        eval_time = clock() - eval_time;
        cout << endl << "On-surface Eval via Singular Quad. evaluation time: " << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 

        // Recover interior limit of solution
        VecAXPY(on_surface_target_value, .5*one, den_solved);

        Test::compute_relative_error(
                exact_on_surface_potential, 
                on_surface_target_value,
                comm, 
                num_target_points, 
                target_dof,
                "On-Surface Evaluation via Singular Quadrature");
        
        //Compare normwise
        //double eval_relative_error;
        Vec eval_difference;
        iC(VecCreateMPI(comm, num_local_targets*target_dof,
                    PETSC_DETERMINE, &eval_difference));
        double zero = 0.;
        VecSet(near_target_value_extrapolate, zero);
        VecSet( near_target_in_out, one);
        // Compute inward facing normals to interpolation along for new near
        // eval
        Vec interpolation_directions;
        VecDuplicate(bis->patch_samples()->sample_point_normal(), &interpolation_directions);
        VecCopy(bis->patch_samples()->sample_point_normal(), interpolation_directions);
        double minus_one = -1.;
        VecScale(interpolation_directions, minus_one);
        // --------------------- Interior Extrapolation ----------------------- 
        eval_time = clock(); 
        stats.start_timer("QBKIX On");
        bis->neaeval_extrapolate(on_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     interpolation_directions,
                     qt,
                     den_solved, 
                     near_target_value_extrapolate);
        stats.stop_timer("QBKIX On");
        stats.store_relative_error(exact_on_surface_potential, 
                near_target_value_extrapolate,
                NORM_INFINITY,
                "QBKIX On-surface evaluation");
        
        Vec near_value_interior_limit;
        VecDuplicate(near_target_value_extrapolate, &near_value_interior_limit);
        VecCopy(near_target_value_extrapolate, near_value_interior_limit);

        eval_time = clock() - eval_time;
        cout << endl << "On-Surface Evaluation via Interior Extrapolation evaluation time: " 
            << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 

        Test::compute_relative_error(
                exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "On-Surface Evaluation via Interior Extrapolation");
        
        

        // ------------ Averaging Interior/Exterior Extrapolation -------------
        VecSet(near_target_value_extrapolate, zero);
        
        eval_time = clock(); /*
        bis->neaeval_extrapolate_average(
                     on_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     interpolation_directions,
                     qt,
                     den_solved, 
                     near_target_value_extrapolate);*/
        eval_time = clock() - eval_time;
        cout << endl << "On-Surface Evaluation via Two-sided Average Extrapolation evaluation time: "
            << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 
        
        //VecAXPY(near_target_value_extrapolate, .5*one, den_solved);
        Test::compute_relative_error(
                exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "On-Surface Evaluation via Two-sided Average Extrapolation");

        // -----------------  Interpolation across surface --------------------
        VecSet(near_target_value_extrapolate, zero);

        eval_time = clock();
        /*
        bis->neaeval_interpolate(
                     on_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     interpolation_directions,
                     qt,
                     den_solved, 
                     near_target_value_extrapolate);
                     */
        //VecAXPY(near_target_value_extrapolate, -.5*one, den_solved);
        eval_time = clock() - eval_time;
        cout << endl << "On-Surface Evaluation via Interpolation across Surface time: " 
            << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 
        Vec near_value_interpolate;
        VecDuplicate(near_target_value_extrapolate, &near_value_interpolate);
        VecCopy(near_target_value_extrapolate, near_value_interpolate);
        
        Test::compute_relative_error(
                exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "On-Surface Evaluation via Interpolation across Surface");
        

    }
    stats.add_result("boundary_distance_ratio", boundary_distance_ratio);
    stats.add_result("interpolation_spacing_ratio", interpolation_spacing_ratio);
    stats.print_results();
    stats.append = true;
    stats.dump();
    iC( VecDestroy(&den_solved) );
    delete bis;
    delete bd;

    // Calling this seg faults the solver at the end
    //iC( PetscFinalize() );
    ebiFunctionReturn(0);
}
