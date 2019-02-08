#include "../catch.hpp"
#include "fmm3d/pvfmm_bis_interface.hpp"
#include "common/kernel3d.hpp"
#include "common/utils.hpp"
#include <random>
using namespace Ebi;
void compute_kifmm_potential(DblNumMat interaction_matrix, DblNumMat source_densities,
        DblNumMat& kifmm_potential){
    int num_sources = interaction_matrix.n();
    int num_targets = interaction_matrix.m();
    int sdof = source_densities.m();
    int tdof = kifmm_potential.m();
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
            cout << kifmm << ", " << pvfmm << endl;
            CHECK(fabs( kifmm - pvfmm)
                    <= fabs(kifmm)*eps + eps);
        }
    }

}

void setup_single_source_target_problem(
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
        VecSet(source_positions, 0.);
        
        int64_t index[3] = {0, 1, 2};
        double values[3]= {0., 0., .5};
        VecSetValues(target_positions, 3, index, values, INSERT_VALUES);
        
        values[2] = 1.;
        VecSetValues(source_normals, 3, index, values, INSERT_VALUES);



        unique_ptr<PvFMM> fmm (new PvFMM());
        fmm->initialize_fmm(source_positions,
                source_normals,
                target_positions, 
                kernel);
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



        PvFMM* fmm = new PvFMM();
        fmm->initialize_fmm(source_positions,
                source_normals,
                target_positions, 
                kernel);
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
        true_potential(0,0) = -1e8;
        
        cout << "true - direct " << endl;
        check_error(true_potential,kifmm_potential);
        cout << "true - fmm" << endl;
        check_error(true_potential,target_potential_local);
        cout << "direct- fmm" << endl;
        check_error(kifmm_potential,target_potential_local);

        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(source_densities);
        Petsc::destroy_vec(target_potential);

}
void random_data_pvfmm_vs_direct(
        Kernel3d kernel,double eps){

        Vec source_positions;
        Vec source_normals;
        Vec source_densities;
        
        Vec target_positions;
        Vec target_potential;
        Vec target_potential_direct;
        
        int num_sources = 500;
        int num_targets = 500;
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
        std::uniform_real_distribution<double> dist(-1e-6, 1e-6);

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
        source_position_local = get_local_vector(
                DIM, num_sources, source_positions);
        
        source_normal_local = get_local_vector(
                DIM, num_sources, source_normals);
        
        source_density_local = get_local_vector(
                sdof, num_sources, source_densities);

        target_position_local = get_local_vector(
                DIM, num_targets, target_positions);

        unique_ptr<PvFMM> fmm (new PvFMM());
        fmm->initialize_fmm(source_positions,
                source_normals,
                target_positions, 
                kernel);
        // Test interaction matrix
        DblNumMat pvfmm_interaction_matrix(num_targets*tdof,num_sources*sdof);
          fmm->interaction_matrix(source_position_local,
                  source_normal_local,
                  target_position_local, 
                  pvfmm_interaction_matrix);

        DblNumMat kifmm_interaction_matrix(
                num_targets*tdof,
                num_sources*sdof);
        kernel.kernel(source_position_local,
                source_normal_local,
                target_position_local, 
                kifmm_interaction_matrix);

        DblNumMat test(
                num_targets*tdof,
                num_sources*sdof);
        for (int i = 0; i < num_targets*tdof; i++) {
             for (int j = 0; j <num_sources*sdof; j++) {
                test(i,j) = pvfmm_interaction_matrix(i,j) - kifmm_interaction_matrix(i,j);         
             }
        }
        check_error(pvfmm_interaction_matrix, kifmm_interaction_matrix, eps);



        source_position_local.restore_local_vector();
        source_normal_local.restore_local_vector();   
        source_density_local.restore_local_vector();  
        target_position_local.restore_local_vector(); 

        fmm->evaluate(source_densities, target_potential);
        fmm->evaluate_direct(source_densities, target_potential_direct);

        DblNumMat fmm_potential_local(tdof, target_potential);
        DblNumMat direct_potential_local(tdof, target_potential_direct);
        
        check_error(fmm_potential_local, direct_potential_local,eps);



        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(source_densities);
        Petsc::destroy_vec(target_potential);

}






TEST_CASE("TestKIFMM/PvFMM Laplace kernel interface", "[kernel][fmm][laplace]"){

    /*SECTION("Laplace Single layer"){

        Kernel3d kernel(LAPLACE + SINGLE_LAYER+ VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
        setup_near_singular_problem(kernel);
    }*/
    SECTION("Laplace Double layer"){

        Kernel3d kernel(LAPLACE + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
        setup_near_singular_problem(kernel);
    }
}
TEST_CASE("TestKIFMM/PvFMM Stokes kernel interface", "[kernel][fmm][stokes]"){
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

TEST_CASE("TestKIFMM/PvFMM Stokes pressure kernel interface", "[kernel][fmm][stokes-press]"){
    Options::set_value_petsc_opts("-bis3d_np","4");
    SECTION("Stokes Single layer pressure"){

        Kernel3d kernel(STOKES + SINGLE_LAYER+ VAR_P, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
    }
    SECTION("Stokes Double layer pressure"){

        Kernel3d kernel(STOKES + DOUBLE_LAYER + VAR_P, vector<double>(2,1.));
        setup_single_source_target_problem(kernel);
    }
}
TEST_CASE("TestKIFMM/PvFMM Navier  kernel interface", "[kernel][fmm][navier]"){
    SECTION("Navier Single layer displacement"){

        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .7;
        equation_coeffs[1] = .4;
        Kernel3d kernel(NAVIER + SINGLE_LAYER+ VAR_U, equation_coeffs);
        setup_single_source_target_problem(kernel);
    }
    SECTION("Navier Double layer displacement"){
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .7;
        equation_coeffs[1] = .4;
        Kernel3d kernel(NAVIER + DOUBLE_LAYER + VAR_U, equation_coeffs);
        setup_single_source_target_problem(kernel);
    }
}
TEST_CASE("TestKIFMM/PvFMM Modified Helmholtz kernel interface", "[kernel][fmm][mod-helm]"){
    Options::set_value_petsc_opts("-np","4");
    SECTION("Modified Helmholtz Single layer"){
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .5;

        Kernel3d kernel(MOD_HELMHOLTZ + SINGLE_LAYER+ VAR_U, equation_coeffs);
        setup_single_source_target_problem(kernel);
    }
    SECTION("Modified Helmholtz Double layer"){
        vector<double> equation_coeffs(2,1.);
        equation_coeffs[0] = .5;

        Kernel3d kernel(MOD_HELMHOLTZ + DOUBLE_LAYER + VAR_U, equation_coeffs);
        setup_single_source_target_problem(kernel);
    }

}
TEST_CASE("Test Pvfmm direct summation via interface", "[fmm][pvfmm-direct]"){
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Kernel3d kernel(LAPLACE + DOUBLE_LAYER + VAR_U, vector<double>(2,1.));
    // eps = 1e-11 and multipole order = 16 seems to pass
    // likely need m=20 for full precision
    random_data_pvfmm_vs_direct(kernel, 1e-11);

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

        unique_ptr<PvFMM> fmm (new PvFMM());
        fmm->initialize_fmm(source_positions,
                source_normals,
                target_positions, 
                kernel);
        fmm->evaluate(density, potential);

        Petsc::destroy_vec(source_positions);
        Petsc::destroy_vec(source_normals);
        Petsc::destroy_vec(target_positions);
        Petsc::destroy_vec(density);
        Petsc::destroy_vec(potential);
    }
}
    
