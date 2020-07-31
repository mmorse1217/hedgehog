#include "../catch.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_samples.hpp"
#include "common/kernel3d.hpp"
#include "common/stats.hpp"
#include <sampling.hpp>
#include <common.hpp>
#include "../utils/evaluation_utils.hpp"
#include "common/vtk_writer.hpp"

using namespace hedgehog;
using Sampling::sample_2d;
using Sampling::equispaced;

TEST_CASE("Navier Neumman condition", "[boundary-data][neumann][navier]"){

    // integral over the boundary with no sources inside the domain should be
    // force-free, i.e. no displacement?
    Options::set_value_petsc_opts("-kt", "511"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
        //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
        //Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
        //Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.0078125");
    Options::set_value_petsc_opts("-bis3d_spacing", ".01");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");
        //Options::set_value_petsc_opts("-upsampling_type", "uniform");
        //Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "3");

    // it's worth noting here that interpolating the blended surface messes up
    // convergence; bottoms out at 1e-8
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0"); 

    Options::set_value_petsc_opts("-bis3d_np", "10");
    vector<double> equation_coefficients(2,1.);
    equation_coefficients[1] = .4;
    Kernel3d kernel(511, equation_coefficients);

    // Set up and sample a surface 
        /*unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
        surface->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        surface->_coarse = true;
        surface->setFromOptions();
        surface->setup();
        surface->refine_test();*/
    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();

    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = vector<int>(surface->patches().size(), 0);
    samples->bdry() = surface.get();
    samples->setup();
    int num_samples = samples->global_num_sample_points(); 
    DblNumMat target_positions = 
        get_local_vector(DIM, num_samples, samples->sample_point_3d_position());
    DblNumMat target_normals= 
        get_local_vector(DIM, num_samples, samples->sample_point_normal());
    DblNumMat quad_weight = 
        get_local_vector(1,num_samples, samples->sample_point_combined_weight());

    SECTION("Test force free navier"){


        int num_sources = 20;

        srand(0);
        DblNumMat source_densities(kernel.get_sdof(), num_sources);
        DblNumMat source_positions(DIM, num_sources);

        double min = -10.;
        double max = 10.;
        for (int i = 0; i < num_sources; i++) {
            for (int d = 0; d < kernel.get_sdof(); d++) {
                // value \in [0,1]
                double value = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                source_densities(d,i) = (max-min)*value + min; // rescale to [min,max]
            }
            for (int d = 0; d < DIM; d++) {
                double value = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                source_positions(d,i) = 1. + 3*value; // rescale to [1,4]
                // ensures that sources are outside the domain
            }
        }

        DblNumMat boundary_data(kernel.get_tdof(), num_samples);

        // compute neumann data
        kernel.neumann_bc_from_singularities(
                source_positions,
                source_densities,
                target_positions, 
                target_normals,
                boundary_data);

        // compute integral over boundary
        DblNumVec integral(kernel.get_tdof());
        setvalue(integral, 0.);
        for (int i = 0; i < num_samples; i++) {
            for (int d = 0; d < kernel.get_tdof(); d++) {
                integral(d) += quad_weight(0,i)*boundary_data(d,i);
            }
        }

        // should be zero 
        for (int d = 0; d < kernel.get_tdof(); d++) {
            cout << log10(fabs(integral(d))) << endl;
            CHECK(fabs(integral(d) ) <=1e-12);
        }
        target_positions.restore_local_vector();
        target_normals.restore_local_vector();
        quad_weight.restore_local_vector();
    }
    SECTION("Test single point force navier"){

        int num_sources = 1;

        srand(1);
        DblNumMat source_densities(kernel.get_sdof(), num_sources);
        DblNumMat source_positions(DIM, num_sources);
        setvalue(source_positions, 0.);

        double min = -10.;
        double max = 10.;
        for (int d = 0; d < kernel.get_sdof(); d++) {
            // value \in [0,1]
            double value = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
            source_densities(d,0) = (max-min)*value + min; // rescale to [min,max]
        }

        DblNumMat boundary_data(kernel.get_tdof(), num_samples);

        // compute neumann data
        kernel.neumann_bc_from_singularities(
                source_positions,
                source_densities,
                target_positions, 
                target_normals,
                boundary_data);

        // compute integral over boundary

        DblNumVec integral(kernel.get_tdof());
        setvalue(integral, 0.);
        for (int i = 0; i < num_samples; i++) {
            for (int d = 0; d < kernel.get_tdof(); d++) {
                integral(d) += quad_weight(0,i)*boundary_data(d,i);
            }
        }

        // should be zero 
        for (int d = 0; d < kernel.get_tdof(); d++) {
            double rel_error = fabs(integral(d) - source_densities(d,0))/fabs(source_densities(d,0));
            cout << integral(d) << ", " << source_densities(d,0) << ", " << log10(rel_error) << endl;
            CHECK(rel_error <=1e-12);
        }
        target_positions.restore_local_vector();
        target_normals.restore_local_vector();
        quad_weight.restore_local_vector();
    }
    SECTION("Test multiple point force navier"){

        int num_sources = 20;

        srand(2);
        DblNumMat source_densities(kernel.get_sdof(), num_sources);
        DblNumMat source_positions(DIM, num_sources);
        setvalue(source_positions, 0.);

        DblNumVec true_integral(kernel.get_tdof());
        setvalue(true_integral, 0.);
        for (int i = 0; i < num_sources; i++) {
            for (int d = 0; d < DIM; d++) {
                double min = -.1;
                double max = .1;
                
                double value = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                source_positions(d,i) = (max-min)*value + min; // rescale to [min,max]
            }

            for (int d = 0; d < kernel.get_sdof(); d++) {
                double min = -10.;
                double max = 10.;
                // value \in [0,1]
                double value = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                source_densities(d,i) = (max-min)*value + min; // rescale to [min,max]
                true_integral(d) += source_densities(d,i);
            }
        }

        DblNumMat boundary_data(kernel.get_tdof(), num_samples);

        // compute neumann data
        kernel.neumann_bc_from_singularities(
                source_positions,
                source_densities,
                target_positions, 
                target_normals,
                boundary_data);

        // compute integral over boundary
        DblNumMat quad_weight = 
            get_local_vector(1, samples->local_num_sample_points(), samples->sample_point_combined_weight());

        DblNumVec integral(kernel.get_tdof());
        setvalue(integral, 0.);
        for (int i = 0; i < num_samples; i++) {
            for (int d = 0; d < kernel.get_tdof(); d++) {
                integral(d) += quad_weight(0,i)*boundary_data(d,i);
            }
        }

        // should be zero 
        for (int d = 0; d < kernel.get_tdof(); d++) {
            double rel_error = fabs(integral(d) - true_integral(d))/fabs(true_integral(d)) ;
            cout << integral(d) << ", " << true_integral(d) << ", " 
                << log10(rel_error) << endl;
            CHECK(rel_error <=1e-12);
        }
        target_positions.restore_local_vector();
        target_normals.restore_local_vector();
        quad_weight.restore_local_vector();
    }
}
TEST_CASE("Dirichlet Neumman condition", "[boundary-data][navier][dirichlet]"){

    // integral over the boundary with no sources inside the domain should be
    // force-free, i.e. no displacement?
    Options::set_value_petsc_opts("-kt", "511"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".05");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");

    Options::set_value_petsc_opts("-bis3d_np", "10");
    vector<double> equation_coefficients(2,1.);
    equation_coefficients[1] = .4;
    Kernel3d kernel(511, equation_coefficients);

    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();

    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = vector<int>(surface->patches().size(), 0);
    samples->bdry() = surface.get();
    samples->setup();
    int num_samples = samples->global_num_sample_points(); 
    
    DblNumMat target_positions = 
        get_local_vector(DIM, num_samples, samples->sample_point_3d_position());

    SECTION("fmm kernel matching dirichlet data"){
        int num_sources = 20;

        srand(2);
        Vec source_densities;
        Vec source_positions;
        Vec source_normals;
        Petsc::create_mpi_vec(MPI_COMM_WORLD, kernel.get_sdof()*num_sources, source_densities);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, DIM*num_sources, source_positions);
        Petsc::create_mpi_vec(MPI_COMM_WORLD, DIM*num_sources, source_normals);
        DblNumMat computed_boundary_data(kernel.get_tdof(), num_samples);
        {
            DblNumMat source_densities_local(kernel.get_sdof(), source_densities);
            DblNumMat source_positions_local(DIM, source_positions);
            setvalue(source_positions_local, 0.);

            for (int i = 0; i < num_sources; i++) {
                for (int d = 0; d < DIM; d++) {
                    double min = -.1;
                    double max = .1;

                    double value = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                    source_positions_local(d,i) = (max-min)*value + min; // rescale to [min,max]
                }

                for (int d = 0; d < kernel.get_sdof(); d++) {
                    double min = -10.;
                    double max = 10.;
                    // value \in [0,1]
                    double value = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                    source_densities_local(d,i) = (max-min)*value + min; // rescale to [min,max]
                }
            }


            // compute neumann data
            kernel.dirichlet_bc_from_singularities(
                    source_positions_local,
                    source_densities_local,
                    target_positions, 
                    computed_boundary_data);
        }
         // evaluate with pvfmm
        unique_ptr<PvFMM> fmm(new PvFMM(
                    source_positions,
                    source_normals,
                    samples->sample_point_3d_position(), kernel));
         Vec true_boundary_data;
         Petsc::create_mpi_vec(MPI_COMM_WORLD, num_samples*kernel.get_tdof(), true_boundary_data);
         VecSet(true_boundary_data,0.);
         fmm->evaluate(source_densities, true_boundary_data);
         DblNumMat true_bd_local(kernel.get_tdof(), true_boundary_data);
         for(int i = 0; i < num_samples; i++){
            for(int d =0; d < kernel.get_tdof(); d++){
                double true_value = true_bd_local(d,i);
                double computed_value = computed_boundary_data(d,i);
                CHECK(fabs(true_value - computed_value) <=1e-6);
            }
         }
    }
}
