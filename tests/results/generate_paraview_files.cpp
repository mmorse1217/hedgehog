#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "../utils/evaluation_utils.hpp"
#include "common/vtk_writer.hpp"
#include "bie3d/markgrid.hpp"

using namespace hedgehog;

TEST_CASE("Generate data for qbkix schematic", "[results][images][schematic]"){
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","1");
    Options::set_value_petsc_opts("-qbkix_convergence_type", "adaptive");
    
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();

    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = vector<int>(face_map->num_patches(),0);
    samples->bdry() = face_map.get();
    samples->setup();
    string prefix = "output/figures/data/";
    // coarse patches
    write_face_map_patches_to_vtk(DblNumMat(0,0), vector<int>(face_map->num_patches(),0),
            face_map.get(), 0, prefix+ "coarse");
    Vec values;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, samples->local_num_sample_points(), values);
    laplace_singluarity_propeller(samples->sample_point_3d_position(), 1, values);
    write_general_points_to_vtk(samples->sample_point_3d_position(),
            1, "cube_coarse_sample_values.vtp", values, prefix);
    VecDestroy(&values);
    
    // surface
    write_face_map_mesh_to_vtk(face_map.get(), 0, prefix+ "coarse_surface_mesh");

    Vec target;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 3*1, target);
    
    DblNumMat t(3,target);
    t(0,0) = 0.1;
    t(1,0) = .68;
    t(2,0) = .53;
    
    vector<int> qbkix_ids;
    for (int i = 0; i < 6; i++) {
       qbkix_ids.push_back(i); 
    }
    NumVec<OnSurfacePoint> closest_on_surface_point = 
        Markgrid::mark_target_points(t, face_map.get(), false);
    
    auto on_surface_point = closest_on_surface_point(0);
    auto patch = face_map->subpatch(on_surface_point.parent_patch);
    
    Vec closest_point;
    Vec normal;
    Vec far_field;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 3*1, closest_point);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 3*1, normal);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, far_field);
    VecSet(far_field, .1);
    DblNumMat cp(3, closest_point);
    DblNumMat n(3, normal);
    
    vector<Point3> pos_and_derivs(3,Point3(0.));
    patch->xy_to_patch_coords(on_surface_point.parametric_coordinates,
            PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, (double*)pos_and_derivs.data());
    
    Point3 point_n = cross(pos_and_derivs[1], pos_and_derivs[2]);
    point_n /= point_n.length();
    for (int i = 0; i < 3; i++) {
        n(i,0) = point_n(i);
        cp(i,0) = pos_and_derivs[0](i);
    }
    n.restore_local_vector();
    cp.restore_local_vector();
    t.restore_local_vector();


    Vec check_points= 
        generate_interior_qbkix_points(
                target,
                closest_point,
                qbkix_ids, 
                normal, 
                far_field,
               far_field 
                );
    VecDestroy(&normal);
    VecDestroy(&far_field);

    face_map->refine_uniform(2);
    Vec target_value;
    Vec check_point_value;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, target_value);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 6*1, check_point_value);
    laplace_singluarity_propeller(target, 1, target_value);
    laplace_singluarity_propeller(check_points, 1, check_point_value);
    cout << "target" << endl;
    write_general_points_to_vtk(target, 1, "cube_schematic_target_value.vtp", target_value, prefix);
    write_general_points_to_vtk(closest_point, 1, "cube_schematic_closest_point.vtp", target_value, prefix);
    VecDestroy(&closest_point);
    cout << "check" << endl;
    write_general_points_to_vtk(check_points, 1, "cube_schematic_check_point_values.vtp", check_point_value, prefix);
    
    // upsampled patches
    write_face_map_patches_to_vtk(DblNumMat(0,0), vector<int>(face_map->num_patches(),0),
            face_map.get(), 0, prefix+ "upsampled");
    write_face_map_mesh_to_vtk(face_map.get(), 0, prefix+ "upsampled_surface_mesh");

    unique_ptr<PatchSamples> samples_fine(new PatchSamples("",""));
    samples_fine->patch_partition() = vector<int>(face_map->num_patches(),0);
    samples_fine->bdry() = face_map.get();
    samples_fine->setup();
    
    Petsc::create_mpi_vec(MPI_COMM_WORLD, samples_fine->local_num_sample_points(), values);
    laplace_singluarity_propeller(samples_fine->sample_point_3d_position(), 1, values);
    write_general_points_to_vtk(samples_fine->sample_point_3d_position(),
            1, "cube_upsampled_sample_values.vtp", values, prefix);
    VecDestroy(&values);

}

TEST_CASE("Generate data for patch coefficients", "[results][images][patch-coeffs]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "6");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","2");
    
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();

    int num_patches = face_map->num_patches();
    Vec coefficients;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_patches*17*17*3, coefficients);
    Vec values;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_patches*17*17, values);
    DblNumMat coeffs(3, coefficients);
    DblNumMat pids(1, values);
    int i =0;
    vector<int> ids;
    for (int pi = 0; pi < face_map->face_map()._patches.size(); pi++) {
        auto patch = face_map->face_map()._patches[pi];
        auto patch_coeffs = patch->_polynomial_patch;
        assert(patch_coeffs->length() == 17*17);
        for (int i = 0; i < patch_coeffs->length(); i++) {
            Point3 p = (*patch_coeffs)(i);
            for (int d = 0; d < 3; d++) {
                coeffs(d,17*17*pi+i) = p(d);
            }
            pids(0,17*17*pi+i) = pi;
        }
        ids.push_back(pi);
    }
    string prefix = "output/figures/data/";
    write_general_points_to_vtk(coefficients,
            1, "cube_patch_coeffs.vtp", values, prefix);
    // surface
    write_face_map_mesh_to_vtk(face_map.get(), 0, prefix+ "cube_patch_coeff_surface_mesh");
    write_face_map_patch_bounding_boxes_to_vtk(ids, face_map.get(),0, prefix + "cube_patch_");
    write_face_map_patch_bounding_boxes_to_vtk(ids, face_map.get(),0, prefix + "cube_patch_inflated", true);




}

void setup_test_case_and_singularities(){

    string prefix = "output/figures/data/";
    
    TestConfig test;
    test.bc_type = BoundaryDataType::HARMONIC;

    test.singularity_type= SingularityType::SPHERE;
    test.sphere_radius_bc = 1.;
    
    // Solve with GMRES using solver_matvec_type
    test.solver_matvec_type = EvaluationType::EXTRAPOLATION_AVERAGE;
    
    test.solution_scheme = SolutionScheme::GMRES_SOLVE;
    test.evaluation_scheme = EvaluationScheme::ON_QBKIX;

    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    if(Options::get_string_from_petsc_opts("-bd3d_filename") == "wrl_files/newtorus.wrl"){
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    } else { 
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    }
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();

    auto solver = setup_solver(face_map.get(), test.solver_matvec_type);
    Vec positions;
    Vec strengths;
    setup_singularities(test, solver, positions, strengths);

    write_face_map_mesh_to_vtk(face_map.get(), 0, prefix+ Test::get_domain() + "_test_");
    write_general_points_to_vtk(positions,
            1, Test::get_domain() + "_test_singularities.vtp", strengths, prefix);
    VecDestroy(&positions);
    VecDestroy(&strengths);
}

TEST_CASE("Generate data for test cases", "[results][images][tests]"){
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","0");
    Options::set_value_petsc_opts("-qbkix_convergence_type", "adaptive");

    /*SECTION("cube"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","2");
        setup_test_case_and_singularities();
    }*/
    SECTION("ttorus2"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/ttorus2.wrl"); // blob
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/ttorus2.wrl");
        setup_test_case_and_singularities();
    }
    /*
    SECTION("pipe"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/pipe.wrl"); // blob
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/pipe.wrl");
        setup_test_case_and_singularities();
    }
    SECTION("newtorus"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl"); // blob
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
        setup_test_case_and_singularities();
    }*/

}

TEST_CASE("Generate data for propeller rhs-refinement", "[results][images][rhs-refine]"){
}
TEST_CASE("Generate data for propeller admissibility", "[results][images][admissibility]"){
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "20");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","1");
    Options::set_value_petsc_opts("-qbkix_convergence_type", "classic");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    
    stats._file_prefix= "output/figures/data/";
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine();
}
    



