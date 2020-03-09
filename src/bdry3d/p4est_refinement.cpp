#include "bdry3d/p4est_refinement.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "common/nummat.hpp"
#include "common/numvec.hpp"
#include "bdry3d/on_surface_point.hpp"
#include "bdry3d/geogram_interface.hpp"
#include "bie3d/markgrid.hpp"
#include "common/stats.hpp"
#include "common/vtk_writer.hpp"
BEGIN_EBI_NAMESPACE


template<class PatchType>
int is_panel_inadmissible(p4est_t* p4est, p4est_topidx_t which_tree, 
        p4est_quadrant_t* quadrant){
  return static_cast<RefinementData<PatchType>*>(quadrant->p.user_data)->refine;
}

void copy_user_data_from_parent_quad(p4est_t* p4est,
        p4est_topidx_t which_tree, 
        int num_outgoing,
        p4est_quadrant_t* outgoing[],
        int num_incoming,
        p4est_quadrant_t* incoming[]){

    assert(num_incoming == 4);
    assert(num_outgoing == 1);
    p4est_quadrant_t* parent = outgoing[0];

    auto parent_ref_data = get_refinement_data(parent);
    assert(parent_ref_data);

    FaceMapSubPatch* parent_patch = parent_ref_data->patch;
    assert(parent_patch);
 
    Rectangle parent_domain(parent_patch->_x_interval,parent_patch->_y_interval);
        
    for(int qi = 0; qi < num_incoming; qi++){
        p4est_quadrant_t* quad = incoming[qi];
        auto child_ref_data = new RefinementData<FaceMapSubPatch>();

        Rectangle domain = parameter_domain_of_child_quad(
                parent_patch->_level, quad, parent_domain);

        FaceMapSubPatch* patch = 
            new FaceMapSubPatch(
                    parent_patch->_face_map_patch, 
                    domain.first, 
                    domain.second, 
                    parent_patch->_level+1,
                    0, 0,
                    parent_patch->_quadrature_weights);

        patch->_quadrant = quad; 

        patch->_coarse_parent_patch = parent_patch->_coarse_parent_patch;
        patch->_quad_id_within_p4est_tree = 0;

        child_ref_data->refine = parent_ref_data->refine;
        child_ref_data->patch = patch;
        quad->p.user_data = child_ref_data;

    }
    //cout << "deleting Patch ["  << parent_patch << "]"  << endl;
    if(!parent_ref_data->refine){
        delete parent_patch;
        // MJM TODO ADD THIS AND TEST
        delete parent_ref_data;
    }
    
}
void balance(p4est_t* p4est){
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, NULL,copy_user_data_from_parent_quad);

}

struct FunctionRefinementData {
    //FunctionHandle f;
    Vec function_values_at_children;
    Vec function_values_at_parent;
    Vec child_sample_points;
    Vec uv_coordinates_single_patch;
    double eps;
    int range_dim;
    int num_samples_1d;
    vector<int> patches_to_refine;
};

extern "C" void update_patch_validity(p4est_iter_volume_info_t * info, void *user_data){
    RefinementData<FaceMapSubPatch>* r = ((RefinementData<FaceMapSubPatch>*)info->quad->p.user_data);
    FunctionRefinementData* f_ref_data = static_cast<FunctionRefinementData*>(user_data);
    FaceMapSubPatch* patch = static_cast<FaceMapSubPatch*>(r->patch);
    auto patches = f_ref_data->patches_to_refine;
    auto it = std::find(patches.begin(), patches.end(), patch->V()) ;
    if(it !=patches.end()){
    assert(patch != NULL);
   
    int num_samples_1d = f_ref_data->num_samples_1d;
    int i = std::distance(patches.begin(), it);
    vector<PetscInt> idx = {i}; // TODO fix for new patch sampling set
    IS parent_data_is;
    ISCreateBlock(MPI_COMM_WORLD,
            f_ref_data->range_dim*num_samples_1d*num_samples_1d,
            1,
            idx.data(),
            PETSC_COPY_VALUES,
            &parent_data_is);
    IS child_data_is;
    ISCreateBlock(MPI_COMM_WORLD,
            4*f_ref_data->range_dim*num_samples_1d*num_samples_1d,
            1,
            idx.data(),
            PETSC_COPY_VALUES,
            &child_data_is);

    Vec parent_data;
    Vec child_data;
    VecGetSubVector(f_ref_data->function_values_at_parent, parent_data_is, &parent_data);
    VecGetSubVector(f_ref_data->function_values_at_children, child_data_is, &child_data);

    r->refine = !(patch->is_patch_valid(parent_data,
                child_data,
                f_ref_data->uv_coordinates_single_patch,
                f_ref_data->range_dim, f_ref_data->eps));
    VecRestoreSubVector(f_ref_data->function_values_at_parent, parent_data_is, &parent_data);
    VecRestoreSubVector(f_ref_data->function_values_at_children, child_data_is, &child_data);
    ISDestroy(&parent_data_is);
    ISDestroy(&child_data_is);
    }
}


void resolve_function(p4est_t* p4est, PatchSurfFaceMap*& face_map, FunctionHandle f, int range_dim, double eps_abs /*, double eps_rel*/){
    // mark all patches in the quadtree forest for refinement
    mark_all_patches_for_refinement<FaceMapSubPatch>(p4est);
    bool patches_still_need_refinement = true;
    // while there are still quads that need refinement
    FunctionRefinementData f_ref_data;
    f_ref_data.eps = eps_abs;
    f_ref_data.range_dim= range_dim;
    f_ref_data.num_samples_1d = floor(1./Options::get_double_from_petsc_opts("-bis3d_spacing"))+1;
    int it = 0;
    while(patches_still_need_refinement){
        unique_ptr<PatchSamples> samples(new PatchSamples("", ""));
        samples->bdry() = face_map;
        samples->patch_partition() = vector<int>(face_map->num_patches(), 0);
        samples->setFromOptions();
        samples->setup();
        // TODO use to only resample patches that need refine rather than all of
        // them...
    
        int num_samples_1d = f_ref_data.num_samples_1d;

        vector<int> invalid_patches_temp = face_map_patch_ids_to_refine(p4est);
        f_ref_data.patches_to_refine = invalid_patches_temp;
        vector<PetscInt> invalid_patches(invalid_patches_temp.begin(), invalid_patches_temp.end());

        int num_patches_to_refine = invalid_patches.size();
        Petsc::create_mpi_vec(face_map->mpiComm(), 
                num_patches_to_refine*
                num_samples_1d*num_samples_1d*range_dim,
                f_ref_data.function_values_at_parent);

        Petsc::create_mpi_vec(face_map->mpiComm(), 
                4*num_patches_to_refine*
                num_samples_1d*num_samples_1d*range_dim,
                f_ref_data.function_values_at_children);

        IS patch_samples_index_set;
        ISCreateBlock(MPI_COMM_WORLD,
                DIM*num_samples_1d*num_samples_1d,
                num_patches_to_refine,
                invalid_patches.data(),
                PETSC_COPY_VALUES,
                &patch_samples_index_set);

        Vec invalid_patches_to_sample_3d_position;
        VecGetSubVector(samples->sample_point_3d_position(), patch_samples_index_set, &invalid_patches_to_sample_3d_position);

        f(invalid_patches_to_sample_3d_position, range_dim, 
                f_ref_data.function_values_at_parent);

        generate_samples_on_child_patches(face_map->mpiComm(), 
                face_map,
                f_ref_data.num_samples_1d,
                f_ref_data.uv_coordinates_single_patch, 
                f_ref_data.child_sample_points,
                invalid_patches);


        // evaluate target function at child sample points
        f(f_ref_data.child_sample_points, range_dim, 
                f_ref_data.function_values_at_children);
        // iterate over the leaves, check if each one is valid via callbacks
    cout << "total num quads to check : " <<invalid_patches.size() << endl;
        p4est_iterate(p4est,
                NULL,
                &f_ref_data,
                update_patch_validity,
                NULL,
                NULL);



/*
        // evaluate target function at parent sample points
        Petsc::create_mpi_vec(face_map->mpiComm(), 
                face_map->num_patches()*
                f_ref_data.num_samples_1d*f_ref_data.num_samples_1d*range_dim,
                f_ref_data.function_values_at_parent);
        Petsc::create_mpi_vec(face_map->mpiComm(), 
                4*face_map->num_patches()*
                f_ref_data.num_samples_1d*f_ref_data.num_samples_1d*range_dim,
                f_ref_data.function_values_at_children);

        f(samples->sample_point_3d_position(), range_dim, 
                f_ref_data.function_values_at_parent);
        generate_samples_on_child_patches(face_map->mpiComm(), 
                face_map,
                f_ref_data.num_samples_1d,
                f_ref_data.uv_coordinates_single_patch, 
                f_ref_data.child_sample_points);
        // evaluate target function at child sample points
        f(f_ref_data.child_sample_points, range_dim, 
                f_ref_data.function_values_at_children);
        // iterate over the leaves, check if each one is valid via callbacks
        p4est_iterate(p4est,
                NULL,
                &f_ref_data,
                update_patch_validity,
                NULL,
                NULL);
*/

        // if so, continue
        //p4est_refine_ext(p4est, 0, -1, is_panel_inadmissible<FaceMapSubPatch>, NULL, copy_user_data_from_parent_quad);
refine_p4est_quads( p4est);
        
        bool dump_points = Options::get_int_from_petsc_opts("-dump_qbkix_points"); 
        if(dump_points){
            DblNumMat test(1,1);
            vector<int> pids(face_map->num_patches(),0);
            for (int i = 0; i < pids.size(); i++) {
                pids[i] = i;
            }
            write_face_map_patches_to_vtk(test, 
                    pids,
                    face_map, it++,
                    "data/rhs_refined_patches");
            string s = "rhs_interp_values_" + to_string(it) + ".vtp";
            write_general_points_to_vtk(invalid_patches_to_sample_3d_position, range_dim, 
                    s, f_ref_data.function_values_at_parent, "data/");
            /*if(it == 5){
                IS is;
                long long idxx[1] = {14};
                ISCreateBlock(MPI_COMM_WORLD,
                        DIM*num_samples_1d*num_samples_1d,
                        num_patches_to_refine,
                        idxx,
                        PETSC_COPY_VALUES,
                        &is);
                Vec tmp;
                VecGetSubVector(f_ref_data.child_sample_points, is, &tmp);

            }*/
        }
        VecRestoreSubVector(samples->sample_point_3d_position(), patch_samples_index_set, &invalid_patches_to_sample_3d_position);
        ISDestroy(&patch_samples_index_set);
        VecDestroy(&f_ref_data.function_values_at_parent);
        VecDestroy(&f_ref_data.function_values_at_children);
        VecDestroy(&f_ref_data.uv_coordinates_single_patch);
        VecDestroy(&f_ref_data.child_sample_points);

        // else, split it using split() inside of the refine_ext callback
        patches_still_need_refinement = check_refinement_criteria<FaceMapSubPatch>(p4est);
        // extract final set of patches of the leaves of the refined forest.
        vector<FaceMapSubPatch*> patches= collect_patches<FaceMapSubPatch>(p4est);
        face_map->patches().assign( patches.begin(), patches.end());
        if(it == 5){ 
            
            //exit(0);
        }
    //vector<Patch*> subpatches = p4est_to_face_map_subpatches(p4est, face_map);
    //face_map->patches() = subpatches;
    }
    //_patches = collect_patches<FaceMapSubPatch>(p4est);

}


void refine_p4est_quads(p4est_t* p4est){
    p4est_refine_ext(p4est, 0, -1, is_panel_inadmissible<FaceMapSubPatch>, NULL, copy_user_data_from_parent_quad);
    update_quad_indicies(p4est);
}


void refine_patches_uniform(int max_level, p4est_t* p4est, PatchSurfFaceMap*& face_map){
    mark_all_patches_for_refinement<FaceMapSubPatch>(p4est);
    for(int level = 0; level < max_level; level++){
        refine_p4est_quads(p4est);
        //p4est_refine_ext(p4est, 0, -1, is_panel_inadmissible<FaceMapSubPatch>, NULL, copy_user_data_from_parent_quad);
    }
}

void refine_patches_for_qbkix_point_location(p4est_t* p4est, PatchSurfFaceMap*& face_map){

    // 
    //initialize_p4est_leaves_with_subpatches(p4est, face_map);
    int qbkix_order = Options::get_int_from_petsc_opts("-near_interpolation_num_samples");
    vector<int> stage_one_qbkix_indices(1, qbkix_order - 1);

    refine_patches_point_location(p4est, face_map, stage_one_qbkix_indices);
    
    vector<int> stage_two_qbkix_indices(1, (qbkix_order-1)/2);
    refine_patches_midpoint_near_medial_axis(p4est, face_map, stage_two_qbkix_indices);

    // TODO ADD ANOTHER STEP THAT CHECK ALL QBKIX POINTS!!!!!
    //vector<int> stage_three_qbkix_indices(qbkix_order);
    //for(int i =0; i <qbkix_order; i++)
    //    stage_three_qbkix_indices[i] = i;
    //refine_patches_midpoint_near_medial_axis(p4est, face_map, stage_three_ref_data);
/*
    // each patch should be refined until the qbkix points from the previous
    // steps are far from all remaining patches
    refine_patches_for_fixed_qbkix_points(p4est, face_map);
*/
    //set_coarse_patch_ids( p4est);
}

void refine_patches_point_location(p4est_t* p4est, PatchSurfFaceMap*& face_map, 
        vector<int> qbkix_indices){
    assert(p4est != NULL);
    assert(face_map != NULL);
    
    PatchSurfFaceMap* intermediate_face_map=NULL;
    Vec qbkix_points = NULL; 
    PatchSamples* test_samples=NULL;

    mark_all_patches_for_refinement<FaceMapSubPatch>(p4est);
    vector<pair<int, int> > patches_to_refine = patches_that_need_refinement(p4est);
    
    int it=0;
    cout << "BEGIN LAST POINT INSIDE DOMAIN REFINEMENT" << endl;
    while (patches_to_refine.size() > 0){

        // updates intermediate_face_map, test_samples and qbkix_points for a
        // given p4est refinement tree
        // qbkix point gen change
        update_and_resample_face_map(p4est, face_map, intermediate_face_map,
                test_samples, qbkix_points, qbkix_indices);


        // qbkix point gen change
        int64_t num_qbkix_samples;
        VecGetSize(qbkix_points, &num_qbkix_samples);
        num_qbkix_samples /= DIM;
        int num_qbkix_samples_per_patch = num_qbkix_samples/patches_to_refine.size();
        
        DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);
        
        //   mark each as inside/outside + near/far
        // @REFACTOR needs to change for tree algo
        //NumVec<OnSurfacePoint> closest_on_surface_points = test_samples->closest_points_to_qbkix_points();
        NumVec<OnSurfacePoint> closest_on_surface_points = 
            Markgrid::mark_target_points(
                    qbkix_points_local,
                    intermediate_face_map,
                    true);

        assert(closest_on_surface_points.m() == qbkix_points_local.n());
        cout << "finished point marking" << endl;
        qbkix_points_local.restore_local_vector();
        // for each patch that needs to be refined...
        for(size_t pi = 0; pi < patches_to_refine.size(); pi++){
            auto refinement_data = 
                get_refinement_data_from_forest(p4est, patches_to_refine[pi]);
            
            if(refinement_data->refine){
                FaceMapSubPatch* patch = refinement_data->patch;
                assert(patch != NULL);
                assert(patch->_quadrant != NULL);
                //assert(patch->_quadrant == quad);
                

                // TODO BUG CHECK this variable is initialized properly...
                // BUG PENDING refactor so that we aren't resampling all patches  each time
                // and finding distances to all qbkix points for all patches...
                // should only need to check the patches that need to be
                // refined....
                int current_patch_qbkix_point_index = 
                    pi*num_qbkix_samples_per_patch;

                //   if any of the QBKIX points are outside, 
                //     mark the patch as "invalid" w.r.t in/out 
                //     i.e. needs to refined
                refinement_data->refine = false;
                for(int qbkix_i = 0; qbkix_i < num_qbkix_samples_per_patch; qbkix_i++){
                    OnSurfacePoint on_surface_point = 
                        closest_on_surface_points(current_patch_qbkix_point_index + qbkix_i);
                    
                    if(on_surface_point.inside_domain == OUTSIDE || 
                            on_surface_point.inside_domain == OUTSIDE_SPATIAL_GRID){
                        refinement_data->refine = true;
                        break;
                    }

                }

            }
        }
        refine_p4est_quads(p4est);
        //p4est_refine_ext(p4est, 0, -1, is_panel_inadmissible<FaceMapSubPatch>, NULL, copy_user_data_from_parent_quad);
        
        // update list of patches to refine and repeat until list is empty
        patches_to_refine = patches_that_need_refinement(p4est);
        qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);

        bool dump_points = Options::get_int_from_petsc_opts("-dump_qbkix_points"); 
        if(dump_points){
            vector<int> ids(intermediate_face_map->num_patches(),0);
            for(size_t i =0; i < ids.size(); i++){
                ids[i] = i;
            }
            string s = stats._file_prefix + "interior_ref_";
            dump_vtk_data_for_paraview(qbkix_points_local, closest_on_surface_points, it++, ids, intermediate_face_map, s);
        }

        qbkix_points_local.restore_local_vector();
    }
}


void refine_patches_midpoint_near_medial_axis(p4est_t*& p4est, PatchSurfFaceMap*& face_map, 
        vector<int> qbkix_indices){
    assert(p4est != NULL);
    assert(face_map != NULL);
    
    PatchSurfFaceMap* intermediate_face_map=NULL;
    Vec qbkix_points = NULL; 
    PatchSamples* test_samples=NULL;

    // mark all patches as needing refinement
    mark_all_patches_for_refinement<FaceMapSubPatch>(p4est);
    vector<pair<int, int> > patches_to_refine = patches_that_need_refinement(p4est);

    bool dump_points = Options::get_int_from_petsc_opts("-dump_qbkix_points"); 
    int it=0;
    while (patches_to_refine.size() > 0){

        // updates intermediate_face_map, test_samples and qbkix_points for a
        // given p4est refinement tree
        update_and_resample_face_map(face_map->_p4est, face_map, intermediate_face_map,
                test_samples, qbkix_points, qbkix_indices);
        

        // Yank out the local qbkix points so we can look at them
        int64_t num_qbkix_samples;
        VecGetSize(qbkix_points, &num_qbkix_samples);
        num_qbkix_samples /= DIM;
        int num_qbkix_samples_per_patch = num_qbkix_samples/patches_to_refine.size();
        
        DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);
        cout << "NUM PATCHES TO REFINE " << patches_to_refine.size() << endl;
          cout << "NUM QBKIX POINTS TO MARK " << num_qbkix_samples << endl;
        //   mark each as inside/outside + near/far
        Markgrid::NearFieldMap closest_on_surface_point_map;
        Markgrid::mark_target_points(
                qbkix_points_local,
                intermediate_face_map,
                closest_on_surface_point_map);
        

        // we should have a map from the ith qbkix point to it's list of closest
        // on-surface points for i = 0, ... num_qbkix_samples
        assert(closest_on_surface_point_map.size() == size_t(qbkix_points_local.n()));

        // Computed closest on-surface point
        NumVec<OnSurfacePoint> final_closest_on_surface_points(num_qbkix_samples);
        setvalue(final_closest_on_surface_points, OnSurfacePoint()); 

        // used in paraview; -1 = patch was not looked at on this refinment
        // pass, i=0,... k = patch is the ith patch inspected on this pass
        vector<int> patches_refined_relative_ids(intermediate_face_map->num_patches(), -1);

        // For each patch we need to refine...
        for(size_t pi = 0; pi < patches_to_refine.size(); pi++){
            
            // Get the p4est tree/quadrant associated with the patch
            auto refinement_data = 
                get_refinement_data_from_forest(p4est, patches_to_refine[pi]);

            // if we need to check the patch for refinement admissibility...
            if(refinement_data->refine){
                FaceMapSubPatch* patch = refinement_data->patch;
                
                // make sure we're looking at the same patch that we asked p4est
                // for
                assert(patch->_parent_id == patches_to_refine[pi].first);
                assert(patch->_quad_id_within_p4est_tree == patches_to_refine[pi].second);
                
                FaceMapSubPatch* p = extract_patch_from_quad(patch->_quadrant); 
                assert(p == patch);

                patches_refined_relative_ids.at(patch->_id) = int(pi);



                // update_and_resample_face_map() produces the set of qbkix
                // points that correspond to patches that need to be refined
                // ONLY. since pi=0... num_patches_to_refine, we can index the
                // QBKIX points this way.
                // TODO write a unit test for this, possible bug...
                int current_patch_qbkix_point_index = 
                    pi*num_qbkix_samples_per_patch;

                refinement_data->refine = false;
                for(int qbkix_i = 0; qbkix_i < num_qbkix_samples_per_patch; qbkix_i++){
                    Point3 target(qbkix_points_local.clmdata(current_patch_qbkix_point_index + qbkix_i));
                    vector<OnSurfacePoint> closest_on_surface_points = 
                        closest_on_surface_point_map[current_patch_qbkix_point_index + qbkix_i];
                    


                    //shuffle(closest_on_surface_points.begin(), closest_on_surface_points.end());
                    /*cout << "near patch ids to point " << qbkix_i << ": " << endl;
                    for(int kk = 0; kk < closest_on_surface_points.size(); kk++){
                        cout << closest_on_surface_points[kk].parent_patch << ", ";
                    }
                    cout << endl;
                    cout << "looking for patch " << patch->V() << endl;
                    */
                    assert(closest_on_surface_points.size() > 0);
                    // find the closest on-surface point
                    // TODO replace with
                    // Markgrid::final_closest_on_surface_points()
                    OnSurfacePoint closest_point =
                        Markgrid::select_closest_point_biased_toward_target_patch(
                            closest_on_surface_points, 
                            patch, target,
                            intermediate_face_map);
                    //cout.precision(16);
                    if(closest_point.parent_patch != patch->V()){
                        cout << "point " << qbkix_i << " on patch " << patch->V() << " caused refinement: " << endl;
                        cout << closest_point.distance_from_target << endl;
                        cout << closest_point.parent_patch << endl;
                        cout << closest_point.target_index << endl;
                        OnSurfacePoint s; s.parent_patch = -1;
                        for(auto p : closest_on_surface_points){
                            if(p.parent_patch == patch->V()){
                                s = p;
                            } 
                        }
                        if(s.parent_patch == -1){

                        cout << "parent patch not found!" << endl;
                        } else {
                            cout << "closest point on patch" << endl;
                        cout << s.distance_from_target << endl;
                        cout << s.parent_patch << endl;
                        cout << s.target_index << endl;
                        cout << "difference in distance: " << fabs(s.distance_from_target - closest_point.distance_from_target) << endl;

                        }

                        refinement_data->refine = true;
                        //closest_point.parent_patch = -200;
                    }
                    final_closest_on_surface_points(current_patch_qbkix_point_index + qbkix_i) = closest_point;
                }

            }
        }
        
        refine_p4est_quads(p4est);
        
        // update list of patches to refine and repeat until list is empty
        patches_to_refine = patches_that_need_refinement(p4est);
        
        if(dump_points){
        string s = stats._file_prefix + "medial_axis_";
            dump_vtk_data_for_paraview(qbkix_points_local,
                    final_closest_on_surface_points, it,
                    patches_refined_relative_ids, intermediate_face_map, s);
            /*
            DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);
            write_qbkix_points_to_vtk(qbkix_points_local, final_closest_on_surface_points,it);
            write_face_map_patches_to_vtk(qbkix_points_local, 
                    patches_refined_relative_ids,
                    intermediate_face_map, it);
            write_lines_from_qbkix_to_closest_point(qbkix_points_local, final_closest_on_surface_points, intermediate_face_map, it);
            qbkix_points_local.restore_local_vector();
             */           
            it++;
        }
        /*if(it == 2)
        exit(0);*/
        qbkix_points_local.restore_local_vector();
    }
}



// Stage 3 refinement
void refine_patches_for_fixed_qbkix_points(p4est_t*& p4est, PatchSurfFaceMap*& face_map){
    assert(p4est != NULL);
    assert(face_map != NULL);
    
    PatchSurfFaceMap* intermediate_face_map=NULL;
    Vec qbkix_points = NULL; 
    PatchSamples* test_samples=NULL;

    // check every patch
    mark_all_patches_for_refinement<FaceMapSubPatch>(p4est, true);

    

    // compute fixed qbkix points to check against
    int qbkix_order = Options::get_int_from_petsc_opts("-near_interpolation_num_samples");
    vector<int> qbkix_indices(qbkix_order);
    for(int i =0; i < qbkix_order; i++)
        qbkix_indices[i] = i;

    update_and_resample_face_map(p4est, face_map, intermediate_face_map,
            test_samples, qbkix_points, qbkix_indices);

    // qbkix point gen change
    int64_t num_qbkix_samples;
    VecGetSize(qbkix_points, &num_qbkix_samples);
    num_qbkix_samples /= DIM;
    //int num_qbkix_samples_per_patch = num_qbkix_samples/intermediate_face_map->num_patches();
    
    DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);
   
    vector<pair<int, int> > patches_to_refine = patches_that_need_refinement(p4est);
    bool dump_points = Options::get_int_from_petsc_opts("-dump_qbkix_points"); 
    string upsampling_algo = 
        Options::get_string_from_petsc_opts("-adaptive_upsampling");// bbox
    int upsampling_switch_iter = 
        Options::get_int_from_petsc_opts("-adaptive_upsampling_switch_iter");

    int it=0;

    while (patches_to_refine.size() > 0){

        //eval_time = omp_get_wtime() - eval_time;
        string s = "adaptive_upsampling_iter_";
        s += to_string(it);
        stats.start_timer(s);
        //stats.add_result("total running time", eval_time);
        // updates intermediate_face_map, test_samples and qbkix_points for a
        // given p4est refinement tree
        
        // make a new grid containing the current face-map
        //Markgrid::SpatialGrid* grid = new Markgrid::SpatialGrid(intermediate_face_map);
        auto subpatches_to_refine = face_map_patches_to_refine(p4est);
        vector<Patch*> parent_subpatches_to_refine(subpatches_to_refine.size());
        parent_subpatches_to_refine.assign(subpatches_to_refine.begin(),
                subpatches_to_refine.end());
        cout << "num patches to upsample: " << parent_subpatches_to_refine.size() << endl;
        unique_ptr<AABBTree> aabb_tree(new AABBTree(parent_subpatches_to_refine,true));
       /* 
        Markgrid::SpatialGrid* grid = new Markgrid::SpatialGrid(parent_subpatches_to_refine);

        // Add all qbkix points from the coarse grid
        for(int qi = 0; qi < num_qbkix_samples; qi++){
            Point3 qbkix_point(qbkix_points_local.clmdata(qi));
            if(qbkix_point < grid->max() && qbkix_point > grid->min()){
                grid->insert(qi, qbkix_point);
            } else {
                assert(0);
            }
        }

        // find the near-zone bounding boxes near to each qbkix point.
        // Below is a map from qbkix point ids to lists of patch bounding boxes
        // that share a grid box with the qi-th qbkix point.
        map<uint, vector<uint> > bounding_boxes_near_points = 
            *grid->boxes_near_points();
        */
        // For each qbkix point...
        mark_all_patches_for_refinement<FaceMapSubPatch>(p4est, false);
        NumVec<OnSurfacePoint> final_closest_on_surface_points(num_qbkix_samples);
        setvalue(final_closest_on_surface_points, OnSurfacePoint());

        bool compute_closest_points = it >= upsampling_switch_iter && upsampling_algo == "bbox_closest_point";
        // TODO make sure this is thread safe
        // TODO BUG closest point computation is not strictly necessary! At best
        // it save one level of refinement. Remove entirely and replace with 
        // grid queries only.
        cout << "checking if qbkix points are valid: compute_closest_points = " << compute_closest_points << endl;
#pragma omp parallel for 
        for(int qi = 0; qi < num_qbkix_samples; qi++){
            // and for each patch whose near-zone overlapped the qi-th qbkix
            // point's grid box...
            final_closest_on_surface_points(qi).parent_patch = -1;

        if(long(100*qi/num_qbkix_samples) % 20 == 0){
            cout << long(100*qi/num_qbkix_samples) <<"%"<< endl;
        }
            Point3 qbkix_point(qbkix_points_local.clmdata(qi));
            // 1. find all nearby refined leaf bounding boxes
            //vector<uint> patches_near_qbkix_point = bounding_boxes_near_points[qi];
            vector<uint> patches_near_qbkix_point = aabb_tree->patches_near_point(qbkix_point);
            /*uint near_patch_id = aabb_tree->closest_patch_to_point(qbkix_point);
            vector<uint> patches_near_qbkix_point;
            patches_near_qbkix_point.push_back(near_patch_id);
            vector<uint> reindexed_patches;
            for(auto pi : patches_near_qbkix_point){
                // Get the actual patch from the face-map
                int patch_id = subpatches_to_refine[pi]->V();
                reindexed_patches.push_back(patch_id);
            }*/
            // 2. find closest points patches
            size_t num_patches_near_point = patches_near_qbkix_point.size();
            vector<OnSurfacePoint> on_surface_points;
            if(compute_closest_points){
                //cout << "computing closest points!" << endl;
                on_surface_points= 
                    Markgrid::compute_closest_points_on_patches(
                            qbkix_point,
                            qi,
                            intermediate_face_map, 
                            patches_near_qbkix_point)[qi];
                assert(num_patches_near_point == on_surface_points.size());
            }
            // 3. mark the leaf boxes as valid/invalid based on leaf surface
            // area
            for(int ppi =0; ppi < patches_near_qbkix_point.size();ppi++){
                // Get the actual patch from the face-map
                int pi = patches_near_qbkix_point[ppi];
                int patch_id = subpatches_to_refine[pi]->V();
                   FaceMapSubPatch* patch = 
                           intermediate_face_map->subpatch(patch_id);
                   
                   /*Point3 bbox_min;
                   Point3 bbox_max;
                   patch->inflated_bounding_box(bbox_min, bbox_max);*/
                   // yank out the correpsonding p4est quad and mark it for
                   // refinement
                   //Markgrid::BoundingBox* bounding_box = grid->_stored_bounding_boxes[pi];
                   //if(qbkix_point > bounding_box->min() && qbkix_point < bounding_box->max()){
               //if(qbkix_point > bbox_min && qbkix_point < bbox_max){
                   int tree_id = patch->_parent_id;
                   int quad_id = patch->_quad_id_within_p4est_tree;
                   auto refinement_data = 
                       get_refinement_data_from_forest(p4est, tree_id, quad_id);
                   if(!compute_closest_points){
                       final_closest_on_surface_points(qi).parent_patch = -10;
                       refinement_data->refine = true; 
                   } else {//if (upsampling_algo == "bbox_closest_point")
                       OnSurfacePoint on_surface_point = on_surface_points[ppi];
                       //cout << "closest point region: " << on_surface_point.region << endl;
                       if(on_surface_point.region == NEAR){
                           //final_closest_on_surface_points(qi).parent_patch = on_surface_point.parent_patch;
                           //final_closest_on_surface_points(qi)= on_surface_point;
                           //cout <<"patch containing closest pt: " <<  on_surface_point.parent_patch << "; " << endl; 
                           //cout << "actual patch we're marking: " << patch_id  <<"; " << endl;
                           refinement_data->refine = true;

                       } else {
                           // far, mark as -10 for paraview
                       final_closest_on_surface_points(qi).parent_patch = -10;

                       }

                   }
                   //}

            }
            //final_closest_on_surface_points(qi) = 
            //    Markgrid::find_closest_on_surface_point_in_list( on_surface_points);
            
            
        }
        stats.stop_timer(s);
        cout << "dump_qbkix_points : " << dump_points << endl;
        if(!dump_points){
            cout << "PRINTING DATA" << endl;
            
            vector<int> patches_refined_relative_ids(intermediate_face_map->num_patches(), -1);
            for(int pi = 0; pi < intermediate_face_map->num_patches(); pi++){
                patches_refined_relative_ids[pi] = pi;

            }
        string s = stats._file_prefix + "upsample_";
            dump_vtk_data_for_paraview(qbkix_points_local,
                    final_closest_on_surface_points, it,
                    patches_refined_relative_ids, intermediate_face_map,
                    s);
            /*
            DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);
            write_qbkix_points_to_vtk(qbkix_points_local, final_closest_on_surface_points,it);
            write_face_map_patches_to_vtk(qbkix_points_local, 
                    patches_refined_relative_ids,
                    intermediate_face_map, it);
            write_lines_from_qbkix_to_closest_point(qbkix_points_local, final_closest_on_surface_points, intermediate_face_map, it);
            qbkix_points_local.restore_local_vector();
            */
            it++;
        }
        
        refine_p4est_quads(p4est);
        //update_face_map(p4est, intermediate_face_map, intermediate_face_map, test_samples);
        update_face_map(p4est, face_map, intermediate_face_map, test_samples);

        // update list of patches to refine and repeat until list is empty
        patches_to_refine = patches_that_need_refinement(p4est);
        //delete grid;
    }

        qbkix_points_local.restore_local_vector();
        VecDestroy(&qbkix_points);
        check_patches(p4est);
}


END_EBI_NAMESPACE
