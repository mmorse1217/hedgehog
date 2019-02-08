#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "evaluation_utils.hpp"
#include "common/vtk_writer.hpp"
using namespace Ebi;

TEST_CASE("dump qbkix points for paraview", "[qbkix-points][debug]"){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "3");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "1000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");

    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type = PatchSurfFaceMap::SurfaceType::POLYNOMIAL;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine_uniform(1);

    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = vector<int>(surface->num_patches(), 0);
    samples->bdry() = surface.get();
    samples->setup();

    vector<int> qbkix_ids(Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    vector<int> patches_to_sample(surface->num_patches());
    std::iota(qbkix_ids.begin(), qbkix_ids.end(), 0);
    std::iota(patches_to_sample.begin(), patches_to_sample.end(), 0);
    
    /*Vec qbkix_points =  samples->generate_qbkix_points_from_sample_points(
            qbkix_ids, patches_to_sample);*/
    Vec normals;
    VecDuplicate(samples->sample_point_normal(), &normals);
    VecCopy(samples->sample_point_normal(), normals);
    VecScale(normals,-1);
   Vec qbkix_points= 
        generate_interior_qbkix_points(
                samples->sample_point_3d_position(),
                samples->sample_point_3d_position(),
                qbkix_ids, 
                normals, 
                samples->sample_point_far_field(),
                samples->sample_point_interpolant_spacing()
                );

    Vec values;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, Petsc::get_vec_size(qbkix_points)/3, values);
    write_general_points_to_vtk(qbkix_points, 1, 
            "qbkix_points.vtp", values, "output/");

}
