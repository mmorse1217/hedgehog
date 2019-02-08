//#include "../bdry3d/bd3dag.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bie3d/evaluator_far.hpp"
#include "bie3d/evaluator_near.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "common/nummat.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/exsol3d.hpp"
#include "common/stats.hpp"
#include "bie3d/solver_utils.hpp"
#include <getopt.h>
#include <string>
#include <iomanip> //for setprecision in file streams
#include <time.h>
// To reinitialize the stored values for regression testing, set this to true
using namespace std;
using namespace Ebi;

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
    VecNorm(eval_difference, NORM_2, &true_minus_computed_norm);

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
        { "face-map",    no_argument,    &is_blended, 2},
        {0,0,0,0},
    };
    int iarg = 0;
    int index;

    int __internal_kernel = 1; // always solving for velocity
    string __domain = "wrl_files/"; //domain to solve problem 
    string __temp_dom;
    int __equation_type; //equation type; string for now (should be int but this should be fine)
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
    string options_file = "opt/morse_test_generic.opt";
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
    
    ebiInit(); eC;

    srand48( (long)time(NULL) );
    MPI_Comm comm; comm = PETSC_COMM_WORLD;
    int mpiRank;   iC( MPI_Comm_rank(comm, &mpiRank) );
    int mpiSize;   iC( MPI_Comm_size(comm, &mpiSize) );  //ebiAssert( mpiSize==1 );
    PetscBool flg = PETSC_FALSE;

    char pbname[100];  iC( PetscOptionsGetString(NULL, "", "-pbname", pbname, 100, &flg) ); ebiAssert(flg==PETSC_TRUE);
    double hs; iC( PetscOptionsGetReal(NULL, "",  "-bis3d_spacing", &hs,  &flg) );  ebiAssert(flg==PETSC_TRUE);
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

    // kernel type for current problem; possible values are in ebi/common/kernel3d.hpp
    int64_t eqn;  iC( PetscOptionsGetInt(NULL, "",  "-kt", &eqn,  &flg) ); 
    ebiAssert(flg==PETSC_TRUE); 

    // LL := number of Lagrange points used to interpolate
    int64_t LL;  iC( PetscOptionsGetInt(NULL, "",  "-LL", &LL,  &flg) ); ebiAssert(flg==PETSC_TRUE);


    SolverGMRESDoubleLayer* bis = new SolverGMRESDoubleLayer("BIS3D_", "bis3d_");
    bis->bdry() = bd;
    bis->patch_partition() = patch_partition;
    bis->dom() = dom;
    //bis->set_evaluation_type(INTERIOR_EXTRAPOLATION);
    bis->eqcoefs() = vector<double>(2,1.0); // 1 for stokes/laplace, 2 for helmholtz
    if (eqn == 511 || eqn == 521){
        bis->eqcoefs()[0] = 1.0;
        bis->eqcoefs()[1] = .5;

    }
    // Load options into BIS
    iC( bis->setFromOptions() );
    cout << "before bis setup" << endl;
    iC( bis->setup() );
    cout << "after bis setup" << endl;

    //3.) Actually solve the equation
    // -et: equation type; possible values are listed in ebi/fmm3d/exsol3d.cpp/hpp
    int64_t bnd;  iC( PetscOptionsGetInt(NULL, "",  "-et", &bnd,  &flg) ); ebiAssert(flg==PETSC_TRUE);

    Exsol3d exsol(eqn, bis->eqcoefs(), bnd); 

    int local_point_num = num_local_points(bis->patch_samples()->sample_point_3d_position());


    int la,lb,lm; 
    int ga,gb,gm; 
    iC( bis->localSize( la,lb,lm) );
    cout << la << ", " << lb << ", " << lm << endl;
    iC( bis->globalSize( ga,gb,gm) );
    
    Kernel3d problem_kernel(eqn+10, bis->eqcoefs());
    int target_dof = problem_kernel.get_tdof();
    int source_dof = problem_kernel.get_sdof();

    double one = 1.;
    double scale = -1.;
    double exact_solution_norm = 0;
    double zero = 0.;
    Vec on_surface_samples;
    VecDuplicate(bis->patch_samples()->sample_point_3d_position(), &on_surface_samples);
    VecCopy(bis->patch_samples()->sample_point_3d_position(), on_surface_samples);
    Vec on_surface_face_points;
    VecDuplicate(bis->patch_samples()->sample_as_face_point(), &on_surface_face_points);
    VecCopy(bis->patch_samples()->sample_as_face_point(), on_surface_face_points);





        //--------------------------------------------------------------------
        // Test on-surface evaluators
        //--------------------------------------------------------------------
        int qt = VAR_U;
        int num_local_targets = bis->patch_samples()->local_num_sample_points();
        Vec den_solved;
        VecCreateMPI(comm,
                num_local_targets*target_dof,
                PETSC_DETERMINE,
                &den_solved);

        Vec near_target_value_extrapolate;
        VecCreateMPI(comm,
                num_local_targets*target_dof,
                PETSC_DETERMINE,
                &near_target_value_extrapolate);

        Vec exact_on_surface_potential;
        iC(VecCreateMPI(comm, 
                    num_local_targets*target_dof, 
                    PETSC_DETERMINE, &exact_on_surface_potential));
        double half = .5;
        VecSet(den_solved, one);
        VecSet(exact_on_surface_potential, half);
        
        Vec on_surface_target_value;
        iC(VecCreateMPI(comm, num_local_targets*target_dof,
                    PETSC_DETERMINE, &on_surface_target_value));

        clock_t eval_time = clock();
        cout << "density size: " << Petsc::get_vec_size(den_solved) << endl;
        cout << "potential size: " << Petsc::get_vec_size(on_surface_target_value) << endl;
        bis->roneval(on_surface_samples,
                     on_surface_face_points,
                     qt, den_solved, on_surface_target_value);
        VecView(on_surface_target_value,PETSC_VIEWER_STDOUT_SELF);
        eval_time = clock() - eval_time;
        cout << "Singular quadrature time: " << ((float)eval_time)/CLOCKS_PER_SEC << endl;
        //VecAXPY(on_surface_target_value, -one, den_solved);

        //VecView(on_surface_target_value, PETSC_VIEWER_STDOUT_SELF);
        //Compare normwise
        compute_relative_error(exact_on_surface_potential, 
                on_surface_target_value, 
                comm, 
                num_local_targets, 
                target_dof,
                "Singular Quadrature");

        //---------------------------------------------------------------
        
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
        //VecDestroy(&val);
        VecSet(val, one);
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
        
        // Store the resulting potential at target points
        Vec target_value;
        iC(VecCreateMPI(comm, num_target_points*target_dof, PETSC_DETERMINE, &target_value));
       
        //iC(VecCreateMPI(comm, num_target_points, PETSC_DETERMINE, &target_value));
        //int qt;
        qt = VAR_U;
        
         eval_time = clock(); 
        bis->fareval(target_points_petsc, qt, den_solved, target_value);
        eval_time = clock() - eval_time;
        cout << endl << "Far evaluation time: " << ((double)eval_time/CLOCKS_PER_SEC) << endl << endl; 
        cout << "in main()" << endl;
        
        compute_relative_error(val, 
                target_value, 
                comm, 
                num_target_points, 
                target_dof,
                "Far Evaluation");
        stats.print_results();
        stats.dump();
        exit(0);

        // Compute || solved_potential - known_potential||_2
        //double exact_solution_norm = 0.0;
        Vec potential_difference;
        iC(VecCreateMPI(comm, target_dof*num_target_points, PETSC_DETERMINE,
                         &potential_difference));
        
        VecDestroy(&potential_difference);
        VecDestroy(&target_value);
        VecDestroy(&target_points_petsc);








        
        // -------------------------------------------------------------------
        // Test 3 "qbkix"-like schemes on the same grid used to solve for
        // density
        // -------------------------------------------------------------------

        Vec near_target_in_out;
        Vec intermediate_targets;

        VecCreateMPI(comm,
                num_local_targets,
                PETSC_DETERMINE,
                &near_target_in_out);

        VecCreateMPI(comm,
                0,
                PETSC_DETERMINE,
                &intermediate_targets);

        VecSet(near_target_value_extrapolate, zero);
        VecSet( near_target_in_out,3*one);
        // Compute inward facing normals to interpolation along for new near
        // eval
        Vec interpolation_directions;
        VecDuplicate(bis->patch_samples()->sample_point_normal(), &interpolation_directions);
        VecCopy(bis->patch_samples()->sample_point_normal(), interpolation_directions);
        double minus_one = -1.;
        VecScale(interpolation_directions, minus_one);
        
        
        // --------------------- Interior Extrapolation ----------------------- 
        eval_time = clock();
        bis->neaeval_extrapolate(on_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     interpolation_directions,
                     qt,
                     den_solved, 
                     near_target_value_extrapolate);
        eval_time = clock() - eval_time;
        cout << "Interior Extrapolation time: " << ((float)eval_time)/CLOCKS_PER_SEC << endl;
        //VecView(near_target_value_extrapolate, PETSC_VIEWER_STDOUT_SELF);
        
        compute_relative_error(exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_local_targets, 
                target_dof,
                "Interior Extrapolation");
        
        // ------------ Averaging Interior/Exterior Extrapolation -------------
        VecSet(near_target_value_extrapolate, zero);
        eval_time = clock();
        bis->neaeval_extrapolate_average(
                     on_surface_samples,
                     intermediate_targets,
                     near_target_in_out, 
                     on_surface_samples,
                     on_surface_face_points,
                     interpolation_directions,
                     qt,
                     den_solved, 
                     near_target_value_extrapolate);
        eval_time = clock() - eval_time;
        cout << "Two-sided Average Extrapolation time: " << ((float)eval_time)/CLOCKS_PER_SEC << endl;
        
        compute_relative_error(exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_local_targets, 
                target_dof,
                "Two-sided Extrapolation Averaged");

        // -----------------  Interpolation across surface --------------------
        VecSet(near_target_value_extrapolate, zero);
        eval_time = clock();

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

        eval_time = clock() - eval_time;
        cout << "Interpolation across Surface time: " << ((float)eval_time)/CLOCKS_PER_SEC << endl;
        
        compute_relative_error(exact_on_surface_potential, 
                near_target_value_extrapolate, 
                comm, 
                num_local_targets, 
                target_dof,
                "Interpolation Across Surface");

        // --------------------------------------------------------------------

    
    iC( VecDestroy(&den_solved) );
    delete bis;
    delete bd;

    // Calling this seg faults the solver at the end
    //iC( PetscFinalize() );
    ebiFunctionReturn(0);
}
