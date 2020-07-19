#include "utils.hpp"
#include "bie3d/solver_utils.hpp"
#include "common/nummat.hpp"
#include <getopt.h>
#include <algorithm>
#include <math.h>
#include "common/ebi_petsc.hpp"
#include <random>
using Ebi::DblNumMat;
using Ebi::get_local_vector;
using Ebi::DIM;
int ArgParse::parse_kernel(char*  arg){
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
int ArgParse::parse_layer(char* arg){
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
void ArgParse::read_command_line_args(int argc, char** argv, string options_file){
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
    //string options_file = "test_optfiles/morse_test_generic.opt";
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
}

Vec Test::generate_random_vector(int size, double upper, double lower){
    std::mt19937_64 mt(1);
    std::uniform_real_distribution<double> dist(upper, lower);
    Vec v;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, size, v); 
    Ebi::DblNumVec vec(v);
    for(int i =0; i < size; i++){
        vec(i) = dist(mt);
    }
    return v;
}

void Test::compute_relative_error(Vec true_potential, Vec computed_potential, 
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

Vec Test::create_scaled_surface_target_points(vector<double> scale_factors, 
        Vec sample_points){
        //PatchSamples* patch_samples){

    //int num_samples = patch_samples->local_num_sample_points();
    int num_samples = Petsc::get_vec_size(sample_points)/DIM;

    // returns: scale_factors.size() parallel copies of the surface, ith copy
    // scaled by scale_factors[i]
    Vec targets;
    VecCreateMPI(
            MPI_COMM_WORLD,
            num_samples*DIM*scale_factors.size(),
            PETSC_DETERMINE,
            &targets);

    DblNumMat targets_local = get_local_vector(DIM, num_samples*scale_factors.size(), targets);
    DblNumMat sample_points_local = get_local_vector(DIM, num_samples, sample_points);

    for(size_t si = 0; si < scale_factors.size(); si++){
        double scaling = scale_factors[si];
        for(int i = 0; i < num_samples; i++){
            int point_index = si*num_samples + i;
            for(int d =0; d < DIM; d++){
                targets_local(d, point_index) = scaling*sample_points_local(d,i);
            }
        }
    }
    targets_local.restore_local_vector();
    sample_points_local.restore_local_vector();
    return targets;
}
Vec Test::compute_dirichlet_boundary_data(
        MPI_Comm comm,
        Ebi::Kernel3d kernel,
        Vec singularity_positions,
        //Vec singularity_normals,
        Vec singularity_densities,
        Vec target_positions,
        Vec target_normals){
    
    Vec boundary_data;
    
    int target_dof = kernel.get_tdof();
    int num_targets = Petsc::get_vec_local_size(target_positions)/DIM;
    int num_sources = Petsc::get_vec_local_size(singularity_positions)/DIM;
    VecCreateMPI(comm, num_targets*target_dof, PETSC_DETERMINE, &boundary_data);

    int num_singularities = Petsc::get_vec_local_size(singularity_positions)/DIM;
    

    DblNumMat source_positions_local = get_local_vector(DIM, num_sources, singularity_positions);
    DblNumMat singularity_densities_local =get_local_vector(target_dof, num_sources, singularity_densities); 
    DblNumMat target_positions_local = get_local_vector(DIM, num_targets, target_positions);
    //DblNumMat target_normals_local = get_local_vector(DIM, num_targets, target_normals);
    DblNumMat boundary_data_local = get_local_vector(target_dof, num_targets, boundary_data);
    setvalue(boundary_data_local, 0.); 
    kernel.dirichlet_bc_from_singularities(source_positions_local,
            singularity_densities_local, target_positions_local, boundary_data_local);
    /*
    for(int ti = 0; ti < num_targets; ti++){
        Point3 x(target_positions_local.clmdata(ti));
        //Point3 n_y(source_normals_local.clmdata(ti));
        for(int si = 0; si < num_sources; si++){
            Point3 y(source_positions_local.clmdata(si));
            Point3 r = x - y;
            double norm_r = r.length();

            // if x is close to y, return 0.
            if(norm_r <= 1e-14)
                norm_r = 1.;

            double coeff = 1./(4*M_PI);
            boundary_data_local(0, ti) += coeff * 1./norm_r;

        }
    }*/
    boundary_data_local.restore_local_vector();
    target_positions_local.restore_local_vector();
    source_positions_local.restore_local_vector();
    return boundary_data;

}
Vec Test::compute_neumann_boundary_data(
        MPI_Comm comm,
        Ebi::Kernel3d kernel,
        Vec singularity_positions,
        //Vec singularity_normals,
        Vec singularity_densities,
        Vec target_positions,
        Vec target_normals){

    //assert(kernel.kernelType() == 111);
    Vec boundary_data;
    
    int source_dof = kernel.get_sdof();
    int num_sources = Petsc::get_vec_local_size(singularity_positions)/DIM;
    
    int target_dof = kernel.get_tdof();
    //assert(target_dof == 1);
    int num_targets = Petsc::get_vec_local_size(target_positions)/DIM;
    VecCreateMPI(comm, num_targets*target_dof, PETSC_DETERMINE, &boundary_data);

    DblNumMat source_positions_local = get_local_vector(DIM, num_sources, singularity_positions);
    DblNumMat singularity_densities_local =get_local_vector(target_dof, num_sources, singularity_densities); 
    DblNumMat target_positions_local = get_local_vector(DIM, num_targets, target_positions);
    DblNumMat target_normals_local = get_local_vector(DIM, num_targets, target_normals);
    DblNumMat boundary_data_local = get_local_vector(target_dof, num_targets, boundary_data);
    setvalue(boundary_data_local, 0.); 
    
    kernel.neumann_bc_from_singularities(source_positions_local,
            singularity_densities_local, 
            target_positions_local, 
            target_normals_local,
            boundary_data_local);
    //setvalue(boundary_data_local, 0.); 
    /*
    for(int ti = 0; ti < num_targets; ti++){
        Point3 x(target_positions_local.clmdata(ti));
        Point3 n_x(target_normals_local.clmdata(ti));
        for(int si = 0; si < num_sources; si++){
            Point3 y(source_positions_local.clmdata(si));
            Point3 r = x - y;
            double norm_r = r.length();

            // if x is close to y, return 0.
            if(norm_r <= 1e-14)
                norm_r = 1.;

            double r3 = norm_r*norm_r*norm_r;
            double r_dot_nx = dot(r, n_x);
            double coeff = 1./(4*M_PI);
            boundary_data_local(0, ti) += coeff * 1./r3 *r_dot_nx;

        }
    }*/
    boundary_data_local.restore_local_vector();
    source_positions_local.restore_local_vector(); 
    singularity_densities_local.restore_local_vector();
    target_positions_local.restore_local_vector();
    target_normals_local.restore_local_vector();
    return boundary_data;
}



bool Options::get_bool_from_petsc_opts(string opt_name){
    PetscBool flg = PETSC_FALSE;
    PetscBool option;
    PetscOptionsGetBool(NULL, "", opt_name.c_str(), &option, &flg);
    assert(flg);
    // replace with c type
    return true ? option :false;

}
string Options::get_string_from_petsc_opts(string opt_name){
    PetscBool flg = PETSC_FALSE;
    string option;
    char option_char[300];
    PetscOptionsGetString(NULL, "", opt_name.c_str(), option_char, 300, &flg);
    assert(flg);
    return string(option_char);

}

string Test::get_domain(){
    string domain = Options::get_string_from_petsc_opts("-bd3d_meshfile");
    domain.erase(domain.begin(), domain.begin()+15); // chop off "wrl_files/" from file name
    domain.erase(domain.end()-4, domain.end()); // chop off ".wrl" from file name
    return domain;
}

string Test::get_kernel(){
    int kernel = Options::get_int_from_petsc_opts("-kt");
    switch(kernel){
        case Ebi::KNL_LAP_S_U:
            return "laplace";
        case Ebi::KNL_STK_S_U:
            return "stokes";
        case Ebi::KNL_NAV_S_U:
            return "navier";
        default:
            assert(0);
    }
}
Point3 Test::cartesian_to_spherical(Point3 p){
    double r = p.length();
    return Point3(r,  acos(p.z()/r), atan2(p.y(), p.x()));
}

double Test::associated_legendre_function(int m, int n, double x){
    assert(m <= 3 && n <=3);
    // Sorry future me.
    // taken from http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
    if(m == -1 && n == 1){
        return -.5*sqrt(1. - x*x);
    } else if(m == 0 && n == 0){
        return 1.;

    } else if(m == 0 && n == 1){
        return x;

    } else if(m == 1 && n == 1){
        return -sqrt(1. - x*x);

    } else if(m == 0 && n == 2){
        return .5*(3.*x*x-1);

    } else if(m == 1 && n == 2){
        return -3*x*sqrt(1 - x*x);

    } else if(m == 2 && n == 2){
       return 3.*(1 -  x*x);
    } else if(m == 0 && n == 3){
        return .5*x*(5*x*x - 3.);
    } else if(m == 1 && n == 3){

        return 3./2.*(1-5*x*x)*sqrt(1 - x*x);
    } else if(m == 2 && n == 3){
        return 15*x*sqrt(1 - x*x);

    } else if(m == 3 && n == 3){
        double t = sqrt(1 - x*x);
        return -15*t*t*t;
    } else {

        assert(0);
    }
}

double Test::spherical_harmonic(int m, int n, Point3 x){
    assert(m <= 3 && n <=3);

    Point3 polar = Test::cartesian_to_spherical(x);
    double r = polar.x();
    double theta = polar.y();
    double phi = polar.z();
    
    double associated_legendre_value = 
        Test::associated_legendre_function(m, n, cos(theta));
    double factorial = 1.;
    for(int i = 0; i < 2*abs(m); i++){
        factorial *= double(n + abs(m) - i);
    }
    //factorial = sqrt(1./factorial);

    double coeff = sqrt((2.*n + 1.)/(4*M_PI)/factorial);
        return coeff * associated_legendre_value * cos(m*phi);
}



Vec Test::compute_spherical_harmonic_bc(MPI_Comm comm, 
        Vec target_positions){
    int num_targets = Petsc::get_vec_size(target_positions)/DIM;
    Vec boundary_data;
    Petsc::create_mpi_vec(comm, num_targets, boundary_data);

    DblNumMat target_positions_local =
        get_local_vector(DIM, num_targets, target_positions);
    DblNumMat boundary_data_local =
        get_local_vector(1, num_targets, boundary_data);

    for(int i =0; i < num_targets; i++){
        Point3 p(target_positions_local.clmdata(i));
        int m =1;
        int n = 1;
        boundary_data_local(0,i) = double(n+1)/double(2*n+1)*(spherical_harmonic(-m,n,p) + spherical_harmonic(m,n, p))/sqrt(2.);
    }

    target_positions_local.restore_local_vector();
    boundary_data_local.restore_local_vector();

    return boundary_data;
}
Vec Test::compute_spherical_harmonic_density(MPI_Comm comm, 
        Vec target_positions){
    int num_targets = Petsc::get_vec_size(target_positions)/DIM;
    Vec density;
    Petsc::create_mpi_vec(comm, num_targets, density);

    DblNumMat target_positions_local =
        get_local_vector(DIM, num_targets, target_positions);
    DblNumMat density_local =
        get_local_vector(1, num_targets, density);
    
    for(int i =0; i < num_targets; i++){
        Point3 p(target_positions_local.clmdata(i));
        int m =1;
        int n = 1;
        density_local(0,i) = (spherical_harmonic(-m,n,p) + spherical_harmonic(m,n, p))/sqrt(2.);
    }

    target_positions_local.restore_local_vector();
    density_local.restore_local_vector();

    return density;
}

double Options::get_double_from_petsc_opts(string opt_name){
    PetscBool flg = PETSC_FALSE;
    double option;
    PetscOptionsGetReal(NULL, "", opt_name.c_str(), &option, &flg);
    assert(flg);
    return option;
}

int Options::get_int_from_petsc_opts(string opt_name){
    PetscBool flg = PETSC_FALSE;
    int64_t option;
    PetscOptionsGetInt(NULL, "", opt_name.c_str(), &option, &flg);
    assert(flg);
    return option;
}


void Options::set_value_petsc_opts(string opt_name, string value){
    PetscOptionsSetValue(NULL, opt_name.c_str(), value.c_str());
}

int Petsc::get_vec_size(Vec v){
    int64_t size;
    VecGetSize(v, &size);
    return size;
}
int Petsc::get_vec_local_size(Vec v){
    int64_t size;
    VecGetLocalSize(v, &size);
    return size;
}

void Petsc::create_mpi_vec(MPI_Comm comm, int64_t size, Vec& v){
    VecCreateMPI(comm, size, PETSC_DETERMINE, &v);
}

void Petsc::destroy_vec(Vec& v){
    VecDestroy(&v);
}

Vec Petsc::concatenate(Vec A, Vec B){
  //---------------------------------------------------------------------------
  // Concatenate these two arrays into one larger array
  // MJM TODO make this less horrifying
  PetscInt A_size;
  PetscInt B_size;
  VecGetLocalSize(A, &A_size);
  VecGetLocalSize(B, &B_size);
  // Construct output vector of size (A.size + B.size)
  Vec C;
  VecCreateMPI(PETSC_COMM_WORLD,
          A_size + B_size,
          PETSC_DETERMINE,
          &C); 

  double* A_ptr;
  double* B_ptr;
  VecGetArray(A, &A_ptr);
  VecGetArray(B, &B_ptr);

  int64_t start_idx, stop_idx;
  VecGetOwnershipRange(C, &start_idx, &stop_idx);

  // MJM BUG this aren't global need vecownershiprange
  // Make global indices for the resulting concatenated Vec
  vector<int64_t> A_index_set;
  vector<int64_t> B_index_set;
  int index = 0;
  for(int i = 0; i < A_size; i++){
      //index = i;
      A_index_set.push_back(start_idx + index++);
  }
  assert(index == A_size);
  for(int i = 0; i < B_size; i++){
      //index = DIM*num_local_A + i;
      B_index_set.push_back(start_idx + index++);
  }

  // set values of the final intermediate target Vec
  VecSetValues(C,
          A_size,
          A_index_set.data(),
          A_ptr,
          INSERT_VALUES);

  VecSetValues(C,
          B_size,
          B_index_set.data(),
          B_ptr,
          INSERT_VALUES);
  
  VecAssemblyBegin(C);
  VecAssemblyEnd(C);

  VecRestoreArray(A, &A_ptr);
  VecRestoreArray(B, &B_ptr);
  return C;
  //---------------------------------------------------------------------------
}

void Debug::save_mat(DblNumMat m, string file){
    cout << file << endl;
    ofstream f;
    f.precision(16);

    f.open(file.c_str());
    for(int i =0; i < m.m(); i++){
        for(int j =0; j < m.n()-1; j++){
            f << m(i,j) << " ";

        }
        f << m(i,m.n()-1) << endl;
    }
    f.close();
}

string Debug::load_file_to_string(string file){
    std::ifstream in(file, std::ios::in);
    if (in){
        std::ostringstream contents;
        contents << in.rdbuf();
        in.close();
        return(contents.str());
    }
    throw(errno);
}

void Debug::load_mat(DblNumMat& m, string file){
    std::ifstream in(file, std::ios::in);
    in.precision(16);
    assert(in.is_open());
    vector<double> stored_mat;
    
    double data = 0.0;
    while(in >> data){
        stored_mat.push_back(data);
    }
    assert(stored_mat.size() == m.n()*m.m());
    
    int it = 0;
    for (int i = 0; i < m.m(); i++) {
        for (int j = 0; j < m.n(); j++) {
           m(i,j) = stored_mat[it++];
        }
    }
}

void Debug::print_mat(DblNumMat m){
    cout << "[ ";
    for(int i =0; i < m.m(); i++){
        cout << "[ ";
        for(int j =0; j < m.n(); j++){
            cout << m(i,j) << ", ";

        }
        cout << "], " << endl;
    }
        cout << "] " << endl;
}


