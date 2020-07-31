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
// To reinitialize the stored values for regression testing, set this to true
using namespace std;
using namespace hedgehog;

// ---------------------------------------------------------------------- 
// A test that compares the stored mat-vec of a random vector with the 
// current computed mat-vec of the system.
// ---------------------------------------------------------------------- 

int parse_kernel(char*  arg){
    string args = string(arg);
    if( args == "laplace"){
        return 100;
    } else if (args == "modhel"){
        return 200;
    } else if (args == "stokes"){
        return 300;
    } else if (args == "navier"){
        return 500;
    } else { 
        cerr << "Kernel not yet implemented" << endl;
        ebiAssert(0);
    }
}
int parse_layer(char* arg){
    string args = string(arg);
    if( args == "single"){
        return 10;
    } else if (args == "double"){
        return 20;
    } else { 
        cerr << "Only single/double layer kernels are valid" << endl;
        ebiAssert(0);
    }
}

void compute_relative_error(Vec true_potential, Vec computed_potential, 
        MPI_Comm comm, int num_local_targets, int target_dof,
        string eval_type){
    cout << endl << eval_type << endl;
    double true_minus_computed_norm = 0.0;
    Vec eval_difference;
    iC(VecCreateMPI(comm, num_local_targets*target_dof,
                PETSC_DETERMINE, &eval_difference));

    VecCopy(true_potential, eval_difference);
    double computed_norm = 0;
    VecNorm(computed_potential, NORM_2, &computed_norm);
    int minus_one = -1.;
    VecAXPY( eval_difference, minus_one,  computed_potential);
    double inf_difference_norm;
    VecNorm(eval_difference, NORM_INFINITY, &inf_difference_norm);
    VecNorm(eval_difference, NORM_2, &true_minus_computed_norm);
    cout << "||computed_potential - exact_potential ||_inf = " << inf_difference_norm << endl;


    double exact_norm = 0.0;
    VecNorm(true_potential, NORM_2, &exact_norm);

    double eval_relative_error = true_minus_computed_norm/exact_norm;
    cout << "||computed_potential - exact_potential ||_2 / ||exact_potential||_2 = "
        << eval_relative_error << endl;
    cout << "||computed_potential||_2 = " << computed_norm << endl;
    cout << "||exact_potential||_2 = " << exact_norm << endl;
    cout << endl  << endl;
    VecDestroy(&eval_difference);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char** argv)
{
    ebiFunctionBegin;
    cout.precision(16);
    // Initialization steps
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
    while (iarg != -1){
        iarg = getopt_long(argc, argv, "k:", longopts, &index);
        switch (iarg){
            case 'k': 
                cout<< "kernel: " << optarg <<  endl;
                __internal_kernel += parse_kernel(optarg);
                break;
            case 'l': 
                cout<< "layer: " << optarg <<  endl;
                __internal_kernel += parse_layer(optarg);
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

        }
    }
    // Load & parse CLI arguments into Petsc

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
    */
    string options_file = argv[1];
    PetscInitialize( &argc, &argv, options_file.c_str(), "Test"); eC;
    
    ebiInit(); eC;

    srand48( (long)time(NULL) );
    MPI_Comm comm; comm = PETSC_COMM_WORLD;
    int mpiRank;   iC( MPI_Comm_rank(comm, &mpiRank) );
    int mpiSize;   iC( MPI_Comm_size(comm, &mpiSize) );  //ebiAssert( mpiSize==1 );
    PetscBool flg = PETSC_FALSE;

    double hs; iC( PetscOptionsGetReal(NULL, "",  "-bis3d_spacing", &hs,  &flg) );  ebiAssert(flg==PETSC_TRUE);
    char pbname[100];  iC( PetscOptionsGetString(NULL, "", "-pbname", pbname, 100, &flg) ); ebiAssert(flg==PETSC_TRUE);
    double sspc; iC( bis3dspacing2bdsurfspacing(hs, sspc) );

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
        iC( bdfm->setFromOptions() );	 
        bd = bdfm;
    } else{ assert(0); }

    // but there is another input option of -bd3d_spacing?
    cout << "getting blended surface " << endl;
    iC( bd->setup() );
    cout << "got blended surface " << endl;
    // bd->patches() is the list of overlapping surface patches 
    int numpch = bd->patches().size();
    cout << "number of patches: " << numpch << endl;
    //double jac;
    vector<int> patch_partition(numpch, 0);  //All in one processor
    /*int stride = numpch/mpiSize;
    for (int i = 0; i < stride; i ++){
        patch_partition[mpiRank*stride+i] = mpiRank;
    }*/

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
    //bis->set_evaluation_type(INTERPOLATION_ACROSS_SURFACE);
    //bis->set_evaluation_type(INTERIOR_EXTRAPOLATION);
    bis->eqcoefs() = vector<double>(2,1.0); // 1 for stokes/laplace, 2 for helmholtz
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

    //int local_point_num = num_local_points(bis->patch_samples()->sample_point_3d_position());


    int la,lb,lm; 
    int ga,gb,gm; 
    iC( bis->localSize( la,lb,lm) );
    cout << la << ", " << lb << ", " << lm << endl;
    iC( bis->globalSize( ga,gb,gm) );
    
    Kernel3d problem_kernel(eqn+10, bis->eqcoefs());
    int target_dof = problem_kernel.get_tdof();
    int source_dof = problem_kernel.get_sdof();
    // Initialize source points: 6 points at +-1*e_i for i=1,2,3 with e_i unit
    // vectors
    
    /*
    Vec sources, source_normals, source_density;
    VecCreateMPI(comm, 6*3, PETSC_DETERMINE, &sources);
    VecCreateMPI(comm, 6*3, PETSC_DETERMINE, &source_normals);
    VecCreateMPI(comm, 6*problem_kernel.get_sdof(), 
                        PETSC_DETERMINE, &source_density);

    double zero = 0.;
    VecSet(sources, zero);
    VecSet(source_normals, zero);
    VecSet(source_density, zero);

    double* src_ptr;
    double* src_normal_ptr;
    double* src_den_ptr;
    VecGetArray(sources, &src_ptr);
    VecGetArray(source_normals, &src_normal_ptr);
    VecGetArray(source_density, &src_den_ptr);
    
    double R = 5.;
    double fourpi = 4*M_PI;
    for(int i=0; i < 3; i++){
        src_ptr[3*i + i] = 1.; // positive unit vector in each dimension
        src_ptr[3*(i + 3) + i] = -1.; //negative unit vector in each dimension
        // normals pointing toward the center
        src_normal_ptr[3*i + i] = -1.; 
        src_normal_ptr[3*(i + 3) + i] = 1.; 
    }
    for(int i =0; i <6; i++){
        for(int j=0; j < source_dof; j++){
            src_den_ptr[source_dof*i + j] = 1.; 
        }
    }
    VecRestoreArray(sources, &src_ptr);
    VecRestoreArray(source_normals, &src_normal_ptr);
    VecRestoreArray(source_density, &src_den_ptr);

    VecScale(sources, R); 
    VecScale(source_density, fourpi);

    {
        PvFMM* fmm = new PvFMM();
        fmm->initialize_fmm(sources, source_normals, 
                bis->patch_samples()->sample_point_3d_position(), problem_kernel);
        fmm->evaluate(source_density, val);
        delete fmm;
    }
*/
    Vec den_solved; iC( VecCreateMPI(comm, lm, PETSC_DETERMINE, &den_solved) );
    //Exsol3d exsol(311, bis->eqcoefs(), 32);
    /*
    Vec val; iC( VecCreateMPI(comm, la, PETSC_DETERMINE, &val) );
    double* val_ptr;
    double* pos_ptr;
    VecGetArray(bis->patch_samples()->sample_point_3d_position(), &pos_ptr);
    VecGetArray(val, &val_ptr);
    int num_local_targets =  num_local_points(bis->patch_samples()->sample_point_3d_position());
    DblNumMat pos_local(DIM, num_local_targets, false, pos_ptr);
    DblNumVec val_local(la, false, val_ptr);
    exsol.quantity(QNT_U, pos_local, val_local);
    VecRestoreArray(bis->patch_samples()->sample_point_3d_position(), &pos_ptr);
    VecRestoreArray(val, &val_ptr);
    */
    Vec val = evaluate_singularities_along_basis( comm,
            problem_kernel, 
            bis->patch_samples()->sample_point_3d_position());
    /*
     Vec val = evaluate_solution_x( comm,
            problem_kernel, 
            bis->patch_samples()->sample_point_3d_position());
            */

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
        //DblNumVec valvec(num_target_points*target_dof,  false, valarr);
        //iC( exsol.quantity(QNT_U, target_points, valvec) ); 
        
        iC( VecRestoreArray(val, &valarr) );
        iC(VecRestoreArray(target_points_petsc, &target_points_petsc_ptr));
        VecDestroy(&val);
        Vec far_field_samples = target_points_petsc;
        /*
        iC(VecCreateMPI(comm, bis->patch_samples()->local_num_sample_points()*3,
                    PETSC_DETERMINE, &far_field_samples));

        VecCopy(bis->patch_samples()->sample_point_3d_position(), far_field_samples);
        PetscScalar far_zone_scale = .6;
        //PetscScalar intermed_zone_scale = .85; 
        VecScale( far_field_samples, far_zone_scale);

        // Store the resulting potential at target points
        Vec target_value;
        iC(VecCreateMPI(comm, bis->patch_samples()->local_num_sample_points()*target_dof, PETSC_DETERMINE, &target_value));
        */ 
        Vec target_value;
        iC(VecCreateMPI(comm, num_target_points*target_dof, PETSC_DETERMINE, &target_value));



        val = evaluate_singularities_along_basis(comm,
            problem_kernel, 
            far_field_samples);
        /*{
        PvFMM* fmm = new PvFMM();
            fmm->initialize_fmm(sources, source_normals, 
                    target_points_petsc, problem_kernel);

            fmm->evaluate(source_density, val);
            delete fmm;
        }*/


        // Initialize far-field evaluator
        // TODO: figure out what on earth these strings in constructors are for
        //VecRestoreArray(target_points_petsc, &target_points_petsc_ptr);
        
        //iC(VecCreateMPI(comm, num_target_points, PETSC_DETERMINE, &target_value));
        int qt;
        qt = VAR_U;
        
        clock_t eval_time = clock(); 
        stats.start_timer("Far Evaluation");
        bis->fareval(far_field_samples, qt, den_solved, target_value);
        stats.stop_timer("Far Evaluation");


        eval_time = clock() - eval_time;
        cout << endl << "Far evaluation time: " << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 
        cout << "in main()" << endl;
        stats.store_relative_error(val, 
                target_value,
                NORM_INFINITY,
                "Smooth quadrature (far evaluation)");
        /*
        compute_relative_error(val, 
                target_value, 
                comm, 
                num_target_points, 
                target_dof,
                "Far Evaluation");*/

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
        PetscScalar near_zone_scale = (1-5e-4);
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
                     qt, den_solved, near_target_value_extrapolate);
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

        //VecDestroy(&near_target_in_out);
        //VecDestroy(&near_targets_as_face_points);

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
        
        //iC( exsol.quantity(QNT_U, near_surface_samples_mat, near_potential_mat) ); 
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

        /*
        {
        PvFMM* fmm = new PvFMM();
        fmm->initialize_fmm(sources, source_normals, 
                near_surface_samples, problem_kernel);

        fmm->evaluate(source_density, exact_near_potential);
        delete fmm;
        }
            */
        num_target_points = patch_samples->local_num_sample_points();
        compute_relative_error(
                exact_near_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "Near Zone Evaluation via Interior Extrapolation");
        /*
        compute_relative_error(
                exact_near_potential, 
                near_target_value, 
                comm, 
                num_target_points, 
                target_dof,
                "Near Zone Evaluation via Interpolation + Singular On-surface Quadrature");

        compute_relative_error(
                near_target_value, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "Error in interior extrapolation compared to legacy near interpolation scheme");
                */

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
        /*
        double* ptr;
        double* pos_ptr;
        VecGetArray(exact_on_surface_potential, &ptr);
        VecGetArray(bis->patch_samples()->sample_point_3d_position(), &pos_ptr);
        DblNumVec local(target_dof* num_local_targets, false, ptr);
        DblNumMat pos_local(DIM, num_local_targets, false, pos_ptr); 
        exsol.quantity(QNT_U, pos_local, local);
        VecRestoreArray(exact_on_surface_potential, &ptr);
        VecRestoreArray(bis->patch_samples()->sample_point_3d_position(), &pos_ptr);
        */
        VecDestroy(&exact_near_potential);

        exact_on_surface_potential = evaluate_singularities_along_basis(comm,
            problem_kernel, 
            on_surface_samples);
        /*{
            PvFMM* fmm = new PvFMM();
            fmm->initialize_fmm(sources, source_normals, 
                    on_surface_samples, problem_kernel);

            fmm->evaluate(source_density, exact_on_surface_potential);
            delete fmm;
        }*/
        //double half = .5;
        //VecSet(den_solved, one);
        //VecSet(exact_on_surface_potential, half);
        
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

        compute_relative_error(
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
        /*
        compute_relative_error(
                exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "On-Surface Evaluation via Interior Extrapolation");
        */ 
        

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
        /*
        compute_relative_error(
                exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "On-Surface Evaluation via Two-sided Average Extrapolation");
                */

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
       /* 
        compute_relative_error(
                exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_target_points, 
                target_dof,
                "On-Surface Evaluation via Interpolation across Surface");
        */
/*
        Vec empty;
        VecCreateMPI(comm, 0, PETSC_DETERMINE, &empty);

        Vec interpolation_points = generate_auxiliary_interpolation_points(
                on_surface_samples,
                on_surface_samples,
                problem_kernel,
                (ExpansionType) INTERPOLATION_ACROSS_SURFACE, // MJM bug remove both enums
                interpolation_directions, NULL);
        Vec interpolation_point_potential;
            VecCreateMPI(comm, 
               8*patch_samples->local_num_sample_points()*target_dof,
               PETSC_DETERMINE,
               &interpolation_point_potential);
        interpolation_point_potential = evaluate_singularities_along_basis(comm,
            problem_kernel, 
            interpolation_points);
        {
            PvFMM* fmm = new PvFMM();
            fmm->initialize_fmm(sources, source_normals, 
                    interpolation_points, problem_kernel);

            fmm->evaluate(source_density, interpolation_point_potential);
            delete fmm;
        }
        write_to_text_file(
                "interpolation_points_true_potential.txt",
               interpolation_points,
               8*patch_samples->local_num_sample_points(),
               empty,
               0,
               interpolation_point_potential,
               target_dof,
               empty,
               0,
               empty,
               0,
               empty,
               0
               );




        write_to_text_file(
                "most_recent_test_values.txt",
               on_surface_samples, // pos
               patch_samples->local_num_sample_points(),
               den_solved, //density
               source_dof,
               on_surface_target_value, //potential
               target_dof,
               exact_on_surface_potential, // val1
               target_dof,
               //near_value_interior_limit,
               //target_dof,
               near_target_value_extrapolate, //val2
               target_dof, 
               near_value_interpolate, //val3
               target_dof);*/

        // --------------------------------------------------------------------
        /*
        VecDestroy(&on_surface_samples);
        VecDestroy(&on_surface_face_points);
        VecDestroy(&near_target_in_out);
        VecDestroy(&near_target_value_extrapolate);

        // Use refined surface as targets
        VecDuplicate(bis->refined_patch_samples()->sample_point_3d_position(),
                &on_surface_samples);
        VecCopy(bis->refined_patch_samples()->sample_point_3d_position(),
                on_surface_samples);
        
        VecDuplicate(bis->refined_patch_samples()->sample_as_face_point(),
                &on_surface_face_points);
        VecCopy(bis->refined_patch_samples()->sample_as_face_point(),
                on_surface_face_points);

        iC(VecCreateMPI(comm, bis->refined_patch_samples()->local_num_sample_points(),
                    PETSC_DETERMINE, &near_target_in_out));
        VecSet( near_target_in_out, one);

        int num_local_target_points = num_local_points(bis->refined_patch_samples()->sample_point_3d_position());
        Vec on_surface_potential;
        iC(VecCreateMPI(comm, target_dof*num_local_target_points, 
                    PETSC_DETERMINE, &on_surface_potential));

        Vec on_surface_target_value;
        iC(VecCreateMPI(comm, bis->refined_patch_samples()->local_num_sample_points()*target_dof, 
                    PETSC_DETERMINE, &on_surface_target_value));
        
        VecDuplicate(on_surface_target_value, &near_target_value_extrapolate);
    
        // Perform near evaluation
        bis->roneval(on_surface_samples,
                     on_surface_face_points,
                     qt, den_solved, on_surface_target_value);
        
        VecSet(near_target_value_extrapolate, zero);
        // Compute inward facing normals to interpolation along for new near
        // eval
        Vec interpolation_directions;
        VecDuplicate(bis->refined_patch_samples()->sample_point_normal(), &interpolation_directions);
        VecCopy(bis->refined_patch_samples()->sample_point_normal(), interpolation_directions);
        double minus_one = -1.;
        VecScale(interpolation_directions, minus_one);

        bis->neaeval_extrapolate(on_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     interpolation_directions,
                     qt,
                     den_solved, 
                     near_target_value_extrapolate);
        
        VecDestroy(&on_surface_face_points);
        // Correct on-surface potential by 1/2*\phi
        double half = .5;

        //MJM TODO call EvaluatorNear::interpolate_density_to_refined_grid in 
        //order to properly correct for the jump

        // Yank out the density at the on surface and neglect poles
        //VecAXPY(on_surface_target_value, half, den_solved);
       
        //double* tmp_arr;
        //double* tmp_den_arr;
        //VecGetArray(on_surface_target_value, &tmp_arr);
        //VecGetArray(den_solved, &tmp_den_arr);
        //for(int i = 0; i < bis->refined_patch_samples()->local_num_sample_points()*target_dof; i++){
        //   tmp_arr[i] += .5*tmp_den_arr[i];
        //}
        //VecRestoreArray(on_surface_target_value, &tmp_arr);
        //VecRestoreArray(den_solved, &tmp_den_arr);
                                
        // Compute exact potential
        double* on_surface_potential_array;
        double* on_surface_samples_array;

        Vec exact_on_surface_potential;
        iC(VecCreateMPI(comm, 
                    bis->refined_patch_samples()->local_num_sample_points()*target_dof, 
                    PETSC_DETERMINE, &exact_on_surface_potential));
        VecGetArray(exact_on_surface_potential, &on_surface_potential_array);
        DblNumVec on_surface_potential_mat(bis->refined_patch_samples()->local_num_sample_points()*target_dof,
                false, on_surface_potential_array);

        VecGetArray(on_surface_samples, &on_surface_samples_array);
        DblNumMat on_surface_samples_mat(3, bis->refined_patch_samples()->local_num_sample_points(), false, 
                                        on_surface_samples_array);
        
        //iC( exsol.quantity(QNT_U, on_surface_samples_mat, on_surface_potential_mat) ); 
        iC( VecRestoreArray(exact_on_surface_potential, &on_surface_potential_array) );
        iC( VecRestoreArray(on_surface_samples, &on_surface_samples_array) );
        {
            PvFMM* fmm = new PvFMM();
            fmm->initialize_fmm(sources, source_normals, 
                    on_surface_samples, problem_kernel);

            fmm->evaluate(source_density, exact_on_surface_potential);
            delete fmm;
        }


        VecGetSize(bis->refined_patch_samples()->sample_point_3d_position(), &size);
        cout << "bis->patch_samples()->sample_point_3d_position() size: " << size << endl;

        VecGetSize(on_surface_samples, &size);
        cout << "on_surface_samples size: " << size << endl;
        
        //Compare normwise
        double on_surface_eval_norm = 0.0;
        Vec on_surface_eval_difference;
        iC(VecCreateMPI(comm, bis->refined_patch_samples()->local_num_sample_points()*target_dof, 
                    PETSC_DETERMINE, &on_surface_eval_difference));

        VecCopy(exact_on_surface_potential, on_surface_eval_difference);
        VecAXPY( on_surface_eval_difference, scale,  on_surface_target_value);
        VecNorm(on_surface_eval_difference, NORM_2, &on_surface_eval_norm);

        exact_solution_norm = 0.0;
        VecNorm(exact_on_surface_potential, NORM_2, &exact_solution_norm);
        
        double on_surface_relative_error = on_surface_eval_norm/exact_solution_norm;

        cout << "|| on_surface_zone_potential_solved - on_surface_zone_potential_exact ||_2 / ||on_surface_zone_potential_exact||_2 = "
             << .5 -on_surface_relative_error << endl;
        cout << "||on_surface_zone_potential_solved||_2 = " << on_surface_eval_norm << endl;
        cout << "||on_surface_zone_potential_exact||_2 = " << exact_solution_norm << endl;
        near_eval_norm = 0.0;

        VecDestroy(&near_eval_difference);

        iC(VecCreateMPI(comm, bis->refined_patch_samples()->local_num_sample_points()*target_dof, 
                    PETSC_DETERMINE, &near_eval_difference));

        VecSet(near_eval_difference, zero);
        VecCopy(exact_on_surface_potential, near_eval_difference);
        VecAXPY( near_eval_difference, scale,  near_target_value_extrapolate);
        VecNorm(near_eval_difference, NORM_2, &near_eval_norm);

        exact_solution_norm = 0.0;
        VecNorm(exact_near_potential, NORM_2, &exact_solution_norm);

        near_zone_relative_error = near_eval_norm/exact_solution_norm;
        cout << "|| on_surface_potential_extrapolate - on_surface_potential_exact ||_2 / ||on_surface_potential_exact||_2 = "
             << near_zone_relative_error << endl;

        
        VecDestroy(&near_target_value);
        VecDestroy(&near_target_value_extrapolate);
        */
    }
    stats.print_results();
    stats.dump();
    iC( VecDestroy(&den_solved) );
    delete bis;
    delete bd;

    // Calling this seg faults the solver at the end
    //iC( PetscFinalize() );
    ebiFunctionReturn(0);
}
