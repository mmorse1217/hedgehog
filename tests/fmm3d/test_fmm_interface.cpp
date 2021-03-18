#include "../catch.hpp"
#include "fmm3d/pvfmm_bis_interface.hpp"
#include "fmm3d/stkfmm_interface.hpp"
#include "common/kernel3d.hpp"
#include "common/utils.hpp"
#include <memory>
#include <petscvec.h>
#include <random>
using namespace hedgehog;
void compute_kifmm_potential(DblNumMat interaction_matrix, DblNumMat source_densities,
        DblNumMat& kifmm_potential){
    int sdof = source_densities.m();
    int tdof = kifmm_potential.m();
    int num_sources = interaction_matrix.n()/sdof;
    int num_targets = interaction_matrix.m()/tdof;
    setvalue(kifmm_potential,0.);
    // Writing this is giving me an ulcer.
    //DblNumMat kifmm_potential(tdof, num_targets);
    for(int si =0; si < num_sources; si++){
        for(int ti =0; ti < num_targets; ti++){
            for(int sd =0; sd < sdof; sd++){
                for(int td =0; td < tdof; td++){
                    int sindex = si*sdof + sd;
                    int tindex = ti*tdof + td;
                    kifmm_potential(td,ti) += 
                        interaction_matrix(tindex,sindex)*source_densities(sd,si);

                }
            }
        }
    }

}

void check_error(DblNumMat kifmm_potential, DblNumMat pvfmm_potential, double eps=1e-14){
    int num_targets = kifmm_potential.n();
    int tdof = kifmm_potential.m();
    assert(pvfmm_potential.n() == num_targets);
    assert(pvfmm_potential.m() == tdof);
    cout.precision(16);
    for(int ti =0; ti < num_targets; ti++){
        for(int td =0; td < tdof; td++){
            double kifmm = kifmm_potential(td, ti);
            double pvfmm = pvfmm_potential(td,ti);
            //cout << kifmm << ", " << pvfmm << endl;
            REQUIRE(fabs( kifmm - pvfmm)
                    <= fabs(kifmm)*eps + eps);
        }
    }

}
template <class Summation>
void setup_single_source_target_problem( Kernel3d kernel){
    Options::set_value_petsc_opts("-bis3d_np", "6");


        Vec source_positions;
        Vec source_normals;
        Vec source_densities;
        
        Vec target_positions;
        Vec target_potential;
        int num_sources = 1;
        int num_targets = 1;
        int tdof = kernel.get_tdof();
        int sdof = kernel.get_sdof();
        cout << "sdof: " << sdof << ", tdof: " << tdof << endl;

        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_normals);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*DIM, target_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*sdof, source_densities);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*tdof, target_potential);

        VecSet(source_densities, 4*M_PI);
        VecSet(source_positions, 0.);
        
        int64_t index[3] = {0, 1, 2};
        double values[3]= {0., 0., .5};
        VecSetValues(target_positions, 3, index, values, INSERT_VALUES);
        
        values[2] = 1.;
        VecSetValues(source_normals, 3, index, values, INSERT_VALUES);



        unique_ptr<Summation> fmm (new Summation(source_positions,
                source_normals,
                target_positions, 
                kernel));
        /*fmm->initialize_fmm(source_positions,
                source_normals,
                target_positions, 
                kernel);*/
        fmm->evaluate(source_densities, target_potential);
        

        // evaluate the corresponding KIFMM kernel directly 
        DblNumMat source_position_local = get_local_vector(
                DIM, num_sources, source_positions);
        
        DblNumMat source_normal_local = get_local_vector(
                DIM, num_sources, source_normals);
        
        DblNumMat source_density_local = get_local_vector(
                sdof, num_sources, source_densities);

        DblNumMat target_position_local = get_local_vector(
                DIM, num_targets, target_positions);


        DblNumMat interaction_matrix(
                num_targets*tdof,
                num_sources*sdof);
        cout << interaction_matrix.n() << "," << interaction_matrix.m() << endl;
        kernel.kernel(source_position_local,source_normal_local,target_position_local, interaction_matrix);

        DblNumMat kifmm_potential(tdof, num_targets);
       
        compute_kifmm_potential(interaction_matrix, 
                source_density_local, kifmm_potential);

        DblNumMat target_potential_local = get_local_vector(
                tdof, num_targets, target_potential);
        
        check_error(kifmm_potential,target_potential_local);

        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(source_densities);
        Petsc::destroy_vec(target_potential);

}
template void setup_single_source_target_problem<PvFMM>( Kernel3d kernel);
template void setup_single_source_target_problem<STKFMM>( Kernel3d kernel);
void setup_near_singular_problem(
        Kernel3d kernel){

        Vec source_positions;
        Vec source_normals;
        Vec source_densities;
        
        Vec target_positions;
        Vec target_potential;
        int num_sources = 1;
        int num_targets = 1;
        int tdof = kernel.get_tdof();
        int sdof = kernel.get_sdof();
        cout << "sdof: " << sdof << ", tdof: " << tdof << endl;

        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_normals);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*DIM, target_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*sdof, source_densities);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*tdof, target_potential);

        VecSet(source_densities, 4*M_PI);
        VecSet(source_positions, .1);
        
        int64_t index[3] = {0, 1, 2};
        double values[3]= {.1, .1, .1 + 1e-4};
        VecSetValues(target_positions, 3, index, values, INSERT_VALUES);
        
        values[2] = 1.;
        VecSetValues(source_normals, 3, index, values, INSERT_VALUES);



        PvFMM* fmm = new PvFMM(source_positions,
                source_normals,
                target_positions, 
                kernel);
        /*fmm->initialize_fmm(source_positions,
                source_normals,
                target_positions, 
                kernel);*/
       fmm->evaluate(source_densities, target_potential);
       // fmm->evaluate_direct(source_densities, target_potential);
        

        // evaluate the corresponding KIFMM kernel directly 
        DblNumMat source_position_local = get_local_vector(
                DIM, num_sources, source_positions);
        
        DblNumMat source_normal_local = get_local_vector(
                DIM, num_sources, source_normals);
        
        DblNumMat source_density_local = get_local_vector(
                sdof, num_sources, source_densities);

        DblNumMat target_position_local = get_local_vector(
                DIM, num_targets, target_positions);


        DblNumMat interaction_matrix(
                num_targets*tdof,
                num_sources*sdof);
       //kernel.kernel(source_position_local,source_normal_local,target_position_local, interaction_matrix);
        fmm->interaction_matrix(source_position_local,source_normal_local,target_position_local, interaction_matrix);

        DblNumMat kifmm_potential(tdof, num_targets);
        
        compute_kifmm_potential(interaction_matrix, 
                source_density_local, kifmm_potential);

        DblNumMat target_potential_local = get_local_vector(
                tdof, num_targets, target_potential);
        DblNumMat true_potential(tdof, num_targets);
        if(kernel.kernel_type() == SINGLE_LAYER){
            true_potential(0,0) = 1e4;
        } else if(kernel.kernel_type() == DOUBLE_LAYER){
            true_potential(0,0) = -1e8;
        }
        
        cout << "true - direct " << endl;
        check_error(true_potential,kifmm_potential, 1e-11);
        cout << "true - fmm" << endl;
        check_error(true_potential,target_potential_local, 1e-11);
        cout << "direct- fmm" << endl;
        check_error(kifmm_potential,target_potential_local, 1e-11);

        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(source_densities);
        Petsc::destroy_vec(target_potential);

}

void test_interaction_matrix(PvFMM* fmm,
                             Vec source_positions,
                             Vec target_positions,
                             Vec source_normals, 
                             Kernel3d kernel,
                             double eps) {
  int sdof = kernel.get_sdof();
  int tdof = kernel.get_tdof();
  int num_targets = Petsc::get_vec_size(target_positions)/DIM;
  int num_sources = Petsc::get_vec_size(source_positions)/DIM;

        DblNumMat source_position_local = get_local_vector(
                DIM, num_sources, source_positions);
        
        DblNumMat source_normal_local = get_local_vector(
                DIM, num_sources, source_normals);
        

        DblNumMat target_position_local = get_local_vector(
                DIM, num_targets, target_positions);

  // Test interaction matrix
  DblNumMat pvfmm_interaction_matrix(num_targets * tdof, num_sources * sdof);
  fmm->interaction_matrix(source_position_local, source_normal_local,
                          target_position_local, pvfmm_interaction_matrix);

  DblNumMat kifmm_interaction_matrix(num_targets * tdof, num_sources * sdof);
  kernel.kernel(source_position_local, source_normal_local,
                target_position_local, kifmm_interaction_matrix);

  // Check error in interaction matrix
  for (int si = 0; si < num_sources; si++) {
    for (int ti = 0; ti < num_targets; ti++) {
      for (int sd = 0; sd < sdof; sd++) {
        for (int td = 0; td < tdof; td++) {
          int sindex = si * sdof + sd;
          int tindex = ti * tdof + td;
          double kifmm = kifmm_interaction_matrix(tindex, sindex);
          double pvfmm = pvfmm_interaction_matrix(tindex, sindex);

          REQUIRE(fabs(kifmm - pvfmm) <= fabs(kifmm) * eps + eps);
        }
      }
    }
  }
}
void test_direct_eval_vs_fmm(PvFMM *fmm, Vec source_positions,
                             Vec source_normals, Vec target_positions,
                             Kernel3d kernel, Vec source_densities,
                             Vec target_potential, Vec target_potential_direct,
                             double eps) {
  int tdof = kernel.get_tdof();
  fmm->evaluate(source_densities, target_potential);

  Options::set_value_petsc_opts("-direct_eval", "1");
  unique_ptr<PvFMM> direct_eval_fmm(new PvFMM(source_positions, source_normals,
                                  target_positions, kernel));

  //direct_eval_fmm->initialize_fmm(source_positions, source_normals,
  //                                target_positions, kernel);
  direct_eval_fmm->evaluate_direct(source_densities, target_potential_direct);
  Options::set_value_petsc_opts("-direct_eval", "0");

  DblNumMat fmm_potential_local(tdof, target_potential);
  DblNumMat direct_potential_local(tdof, target_potential_direct);

  check_error(fmm_potential_local, direct_potential_local, eps);
}

void setup_random_test_data(Kernel3d kernel,
        Vec& source_positions,
        Vec& source_normals,
        Vec& source_densities,
        
        Vec& target_positions,
        Vec& target_potential,
        Vec& target_potential_direct){
        int num_sources = 250;
        int num_targets = 250;
        int tdof = kernel.get_tdof();
        int sdof = kernel.get_sdof();

        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_normals);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*DIM, target_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*sdof, source_densities);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*tdof, target_potential);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*tdof, target_potential_direct);

        DblNumMat source_position_local = get_local_vector(
                DIM, num_sources, source_positions);
        
        DblNumMat source_normal_local = get_local_vector(
                DIM, num_sources, source_normals);
        
        DblNumMat source_density_local = get_local_vector(
                sdof, num_sources, source_densities);

        DblNumMat target_position_local = get_local_vector(
                DIM, num_targets, target_positions);
        std::mt19937_64 mt(0);
        std::uniform_real_distribution<double> dist(-1, 1);

        for (int i = 0; i < num_sources; i++) {
            for (int d = 0; d < DIM; d++) {
                source_position_local(d,i) = dist(mt);
                source_normal_local(d,i) = dist(mt);
            }   
            for (int d = 0; d < sdof; d++) {
                source_density_local(d,i) = dist(mt);
            }   
        }
        for (int i = 0; i < num_targets; i++) {
            for (int d = 0; d < DIM; d++) {
                target_position_local(d,i) = dist(mt);
            }  
        }
        source_position_local.restore_local_vector();
        source_normal_local.restore_local_vector();   
        source_density_local.restore_local_vector();  
        target_position_local.restore_local_vector(); 
}
void random_data_pvfmm_vs_direct(
        Kernel3d kernel,double eps){

        Vec source_positions;
        Vec source_normals;
        Vec source_densities;
        
        Vec target_positions;
        Vec target_potential;
        Vec target_potential_direct;
        
        setup_random_test_data(kernel, source_positions, source_normals,
                               source_densities, target_positions,
                               target_potential, target_potential_direct);

        unique_ptr<PvFMM> fmm(new PvFMM(source_positions, source_normals,
                                        target_positions, kernel));
        //fmm->initialize_fmm(source_positions,
        //        source_normals,
        //        target_positions, 
        //        kernel);
        test_interaction_matrix(fmm.get(), source_positions,
                                target_positions, source_normals,
                                kernel, eps);


        test_direct_eval_vs_fmm(fmm.get(), source_positions, source_normals,
                                target_positions, kernel, source_densities,
                                target_potential, target_potential_direct, eps);

        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(source_densities);
        Petsc::destroy_vec(target_potential);
}




// ----------------------------------------------------------------------------
//  PVFMM tests
// ----------------------------------------------------------------------------


TEST_CASE("Test PvFMM Laplace kernel interface", "[kernel][fmm][laplace][pvfmm]"){

    SECTION("Laplace Single layer"){

        Kernel3d kernel(LAPLACE + SINGLE_LAYER+ VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem<PvFMM>(kernel);
        setup_near_singular_problem(kernel);
    }
    SECTION("Laplace Double layer"){

        Kernel3d kernel(LAPLACE + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem<PvFMM>(kernel);
        setup_near_singular_problem(kernel);
    }
}
TEST_CASE("Test PvFMM Stokes kernel interface", "[kernel][fmm][stokes][pvfmm]"){
    SECTION("Stokes Single layer velocity"){

    Options::set_value_petsc_opts("-bis3d_np","4");
        Kernel3d kernel(STOKES + SINGLE_LAYER+ VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem<PvFMM>(kernel);
    }
    SECTION("Stokes Double layer velocity"){

    Options::set_value_petsc_opts("-bis3d_np","4");
        Kernel3d kernel(STOKES + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem<PvFMM>(kernel);
    }
}

TEST_CASE("Test PvFMM Stokes pressure kernel interface", "[kernel][fmm][stokes-press][pvfmm]"){
    Options::set_value_petsc_opts("-bis3d_np","4");
    SECTION("Stokes Single layer pressure"){

        Kernel3d kernel(STOKES + SINGLE_LAYER+ VAR_P, vector<double>(2,1.));
        setup_single_source_target_problem<PvFMM>(kernel);
    }
    SECTION("Stokes Double layer pressure"){

        Kernel3d kernel(STOKES + DOUBLE_LAYER + VAR_P, vector<double>(2,1.));
        setup_single_source_target_problem<PvFMM>(kernel);
    }
}
TEST_CASE("Test PvFMM Navier  kernel interface", "[kernel][fmm][navier][pvfmm]"){
    Options::set_value_petsc_opts("-np","4");
    SECTION("Navier Single layer displacement"){

        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .7;
        equation_coeffs[1] = .4;
        Kernel3d kernel(NAVIER + SINGLE_LAYER+ VAR_U, equation_coeffs);
        setup_single_source_target_problem<PvFMM>(kernel);
    }
    SECTION("Navier Double layer displacement"){
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .7;
        equation_coeffs[1] = .4;
        Kernel3d kernel(NAVIER + DOUBLE_LAYER + VAR_U, equation_coeffs);
        setup_single_source_target_problem<PvFMM>(kernel);
    }
}
TEST_CASE("Test PvFMM Modified Helmholtz kernel interface", "[kernel][fmm][mod-helm][pvfmm]"){
    Options::set_value_petsc_opts("-np","4");
    SECTION("Modified Helmholtz Single layer"){
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .5;

        Kernel3d kernel(MOD_HELMHOLTZ + SINGLE_LAYER+ VAR_U, equation_coeffs);
        setup_single_source_target_problem<PvFMM>(kernel);
    }
    SECTION("Modified Helmholtz Double layer"){
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .5;

        Kernel3d kernel(MOD_HELMHOLTZ + DOUBLE_LAYER + VAR_U, equation_coeffs);
        setup_single_source_target_problem<PvFMM>(kernel);
    }

}
TEST_CASE("Test Pvfmm direct summation via interface", "[fmm][pvfmm][direct]"){
    Options::set_value_petsc_opts("-bis3d_np", "12");
    Kernel3d kernel;
    double eps = 1e-7;
    SECTION("Laplace Single Layer"){
        kernel = Kernel3d(LAPLACE + SINGLE_LAYER + VAR_U, vector<double>(2,1.));
    }
    SECTION("Laplace Double Layer"){
        kernel = Kernel3d(LAPLACE + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
    }
    SECTION("Stokes Single Layer"){
        Options::set_value_petsc_opts("-bis3d_np", "10");
        eps=1e-5;
        kernel = Kernel3d(STOKES + SINGLE_LAYER + VAR_U, vector<double>(2,1.));
    }
    SECTION("Stokes Double Layer"){
        Options::set_value_petsc_opts("-bis3d_np", "10");
        eps=1e-5;
        kernel = Kernel3d(STOKES + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
    }
    SECTION("Navier Single Layer"){
        Options::set_value_petsc_opts("-bis3d_np", "10");
        eps=1e-5;
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .7;
        equation_coeffs[1] = .4;
        kernel = Kernel3d(NAVIER + SINGLE_LAYER + VAR_U, equation_coeffs);
    }
    SECTION("Navier Double Layer"){
        Options::set_value_petsc_opts("-bis3d_np", "10");
        eps=1e-5;
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .7;
        equation_coeffs[1] = .4;
        kernel = Kernel3d(NAVIER + DOUBLE_LAYER + VAR_U, equation_coeffs);
    }
    // eps = 1e-11 and multipole order = 16 seems to pass
    // likely need m=20 for full precision
    random_data_pvfmm_vs_direct(kernel, eps);

}
// ----------------------------------------------------------------------------
//  STKFMM tests
// ----------------------------------------------------------------------------

TEST_CASE("Test STKFMM Laplace kernel interface", "[kernel][fmm][laplace][stkfmm]"){

    SECTION("Laplace Single layer"){

        Kernel3d kernel(LAPLACE + SINGLE_LAYER+ VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem<STKFMM>(kernel);
        setup_near_singular_problem(kernel);
    }
    SECTION("Laplace Double layer"){

        Kernel3d kernel(LAPLACE + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem<STKFMM>(kernel);
        setup_near_singular_problem(kernel);
    }
}
/*TEST_CASE("Test STKFMM Stokes kernel interface", "[kernel][fmm][stokes][stkfmm]"){
    SECTION("Stokes Single layer velocity"){

    Options::set_value_petsc_opts("-bis3d_np","4");
        Kernel3d kernel(STOKES + SINGLE_LAYER+ VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
    }
    SECTION("Stokes Double layer velocity"){

    Options::set_value_petsc_opts("-bis3d_np","4");
        Kernel3d kernel(STOKES + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
    }
}

TEST_CASE("Test STKFMM Stokes pressure kernel interface", "[kernel][fmm][stokes-press][stkfmm]"){
    Options::set_value_petsc_opts("-bis3d_np","4");
    SECTION("Stokes Single layer pressure"){

        Kernel3d kernel(STOKES + SINGLE_LAYER+ VAR_P, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
    }
    SECTION("Stokes Double layer pressure"){

        Kernel3d kernel(STOKES + DOUBLE_LAYER + VAR_P, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
    }
}*/


TEST_CASE("Debug fmm in a loop", "[fmm][loop][memory]"){
    int num_sources = 100;
    int num_targets= 100;
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5");
    Options::set_value_petsc_opts("-bis3d_np", "12");
    Kernel3d kernel(121, vector<double>(2,0.));
    PetscRandom seed;
    PetscRandomCreate(MPI_COMM_WORLD, &seed);
        Vec source_positions;
        Vec source_normals;
        Vec target_positions;
        Vec density;
        Vec potential;
        Vec potential_direct;

        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_normals);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*DIM, target_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*kernel.get_sdof(), density);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*kernel.get_tdof(), potential);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*kernel.get_tdof(), potential_direct);
        VecSetRandom(source_positions, seed);
        VecSetRandom(source_normals, seed);
        VecSetRandom(target_positions, seed);
        VecSetRandom(potential, seed);
        unique_ptr<PvFMM> fmm (new PvFMM(source_positions,
                source_normals,
                target_positions, 
                kernel));
//        fmm->initialize_fmm(source_positions,
//                source_normals,
//                target_positions, 
//                kernel);
    for (int i = 0; i < 10; i++) {
        VecSetRandom(density, seed);
        VecSet(potential, 0.);
        VecSet(potential_direct, 0.);


        fmm->evaluate(density, potential);
        Options::set_value_petsc_opts("-direct_eval", "1");
        fmm->evaluate_direct(density, potential_direct);
        Options::set_value_petsc_opts("-direct_eval", "0");
        DblNumMat pot(kernel.get_tdof(), potential);
        DblNumMat pot_direct(kernel.get_tdof(), potential_direct);
        
        check_error(pot, pot_direct, 1e-5);

        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(density);
        Petsc::destroy_vec(potential);
    }
}
    
TEST_CASE("Debug pvfmm memory leak", "[fmm][debug][memory]"){
    int num_sources = 1000;
    int num_targets= 1000;
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5");
    Kernel3d kernel(121, vector<double>(2,0.));
    PetscRandom seed;
    PetscRandomCreate(MPI_COMM_WORLD, &seed);
    for (int i = 0; i < 100; i++) {

        Vec source_positions;
        Vec source_normals;
        Vec target_positions;
        Vec density;
        Vec potential;

        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*DIM, source_normals);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*DIM, target_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_sources*kernel.get_sdof(), density);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*kernel.get_tdof(), potential);
        VecSetRandom(source_positions, seed);
        VecSetRandom(source_normals, seed);
        VecSetRandom(target_positions, seed);
        VecSetRandom(density, seed);
        VecSetRandom(potential, seed);

        unique_ptr<PvFMM> fmm (new PvFMM(
        source_positions,
                source_normals,
                target_positions, 
                kernel));
        fmm->evaluate(density, potential);

        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(density);
        Petsc::destroy_vec(potential);
    }
}
    
