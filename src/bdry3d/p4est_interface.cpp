#include "p4est_interface.hpp" 
#include "common/ebi.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bie3d/evaluator_near.hpp"
#include "bie3d/markgrid.hpp"
#include "bie3d/solver_utils.hpp"
#include "common/vtk_writer.hpp"
#include <sampling.hpp>
#include <bitset>
#include <iomanip>

extern "C" {
#include <p4est_search.h>
#include <sc_containers.h>
}

BEGIN_EBI_NAMESPACE

void initialize_user_pointer(p4est_t* p4est,
        p4est_topidx_t which_tree, 
        p4est_quadrant_t* quadrant){
        quadrant->p.user_data = new RefinementData<FaceMapSubPatch>();
}

FaceMapSubPatch* extract_patch_from_quad(p4est_quadrant_t* quad){
    auto r = static_cast<RefinementData<FaceMapSubPatch>*>(quad->p.user_data);
    return r->patch;
}

RefinementData<FaceMapSubPatch>* get_refinement_data_from_forest(p4est_t* p4est,
        int tree_id, int quad_id_within_tree){
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, tree_id);
    p4est_quadrant_t* quad = p4est_quadrant_array_index(&tree->quadrants, quad_id_within_tree);
    return static_cast<RefinementData<FaceMapSubPatch>*>(quad->p.user_data);
}

RefinementData<FaceMapSubPatch>* get_refinement_data_from_forest(p4est_t* p4est,
        pair<int,int> tree_and_quad_id){
    return get_refinement_data_from_forest(p4est, tree_and_quad_id.first, tree_and_quad_id.second);
}

RefinementData<FaceMapSubPatch>* get_refinement_data(p4est_quadrant_t* quad){
    return static_cast<RefinementData<FaceMapSubPatch>*>(quad->p.user_data);
}

void dump_vtk_data_for_paraview(DblNumMat qbkix_points,
        NumVec<OnSurfacePoint> closest_on_surface_points, int it,
        vector<int> global_patch_ids, PatchSurfFaceMap* face_map,
        string file_prefix ){
    
    // qbkix_points - qbkix points generated during refinement
    // closest_on_surface_points - closest on surface points to qbkix_points,
    //                              elementwise. thuse
    //                              len(closest_on_surface_points) ==
    //                              len(qbkix_points)
    // it - refinement iteration 
    // global_patch_ids - enumeration of patch ids on the it-th level of
    //                      refinement. This should match subpatch->_id for
    //                      subpatch in face_map
    // face_map - face_map that generated all of the associated input data
    write_qbkix_points_to_vtk(qbkix_points, closest_on_surface_points,it, file_prefix);
    write_face_map_patches_to_vtk(qbkix_points, global_patch_ids, face_map, it, file_prefix);
    //write_face_map_patch_bounding_boxes_to_vtk( global_patch_ids, face_map, it, file_prefix);
    write_lines_from_qbkix_to_closest_point(qbkix_points, closest_on_surface_points, face_map, it, file_prefix);
    write_face_map_patch_bounding_boxes_to_vtk(global_patch_ids, face_map, it, file_prefix, true);
    write_face_map_mesh_to_vtk(face_map, it, file_prefix, 11);

}

vector<FaceMapSubPatch*> collect_face_map_subpatches(p4est_t* p4est){
    vector<FaceMapSubPatch*> patches= collect_patches<FaceMapSubPatch>(p4est);
    return patches;
}


Rectangle parameter_domain_of_child_quad(int parent_patch_level, p4est_quadrant_t* quad, Rectangle parent_domain){
    int8_t level = int8_t(parent_patch_level+1);
    // shift significant bits defining quad's location within it's parent to the
    // 1's place
    int xcoord_int = quad->x >> (P4EST_MAXLEVEL - quad->level);
    int ycoord_int = quad->y >> (P4EST_MAXLEVEL - quad->level);

    // yank out the x-/y-location of the quad w.r.t its parent (i.e. 
    // (0,0), (0,1), (1,0), or (1,1)
    xcoord_int = xcoord_int & 1;
    ycoord_int = ycoord_int & 1; 

    // Need to scale integer coordinates by the width of the child subdomain 
    // and shift by theparent's x/y-param lower bounds to properly position
    // child param domain
    double x_min = double(xcoord_int)/double(pow2(level)) + parent_domain.first.first;
    double y_min = double(ycoord_int)/double(pow2(level)) + parent_domain.second.first;

    double x_max = x_min + 1./double(pow2(level));
    double y_max = y_min + 1./double(pow2(level));

    Rectangle child_domain(Interval(x_min, x_max), Interval(y_min, y_max));

    return child_domain;
}

extern "C" int point_in_quad_wrapper(p4est_t * p4est,
            p4est_topidx_t which_tree,
            p4est_quadrant_t * quadrant,
            p4est_locidx_t local_num,
            void *point){
   PointInQuad* functor = static_cast<PointInQuad*>(p4est->user_pointer);
   return (*functor)(p4est, which_tree, quadrant, local_num, point);

}

void find_leaf_quad_containing_points(p4est_t* p4est, vector<OnSurfacePoint> points_to_process){
    sc_array_t* a = sc_array_new(sizeof(pair<int, OnSurfacePoint>));
    sc_array_init(a, sizeof(pair<int, OnSurfacePoint>));
    
    for(int i =0; i < points_to_process.size(); i++){
        pair<int,OnSurfacePoint>* ptr = (pair<int,OnSurfacePoint>*) sc_array_push(a);
        *ptr = pair<int,OnSurfacePoint>(i, points_to_process[i]);

    }

   assert(0); //doesn't allow for concurrent access...
   PointInQuad* functor = static_cast<PointInQuad*>(p4est->user_pointer);
   functor->quadrant_containing_points.resize(points_to_process.size());

    p4est_search (p4est,
            NULL,
            point_in_quad_wrapper,
            a);
    for(int i =0; i < points_to_process.size(); i++){

    }
    sc_array_reset(a);
}
struct Storage{
    PatchSurfFaceMap* face_map;
    int quad_counter;

};
extern "C" void reinitialize_subpatches(p4est_iter_volume_info_t * info, void *user_data){
    p4est_quadrant_t* quad = info->quad;
    int tree_id = info->treeid;
    
    auto r = get_refinement_data(quad);
    // only call this function initially before any refinement or processing...
    
    // 0th level w.r.t to quadrature refinement. important 
    // distinction between level in the p4est tree and level of quadrature
    // refinement, which takes an existing tree as input.
    int level = 0; 
    Rectangle domain;
    quad_to_face_map_uv_bounds(quad, domain, level);
    
    Storage* s = static_cast<Storage*>(user_data); 
    int face_map_id = s->quad_counter++;

    PatchSurfFaceMap* face_map = s->face_map;
    FaceMapPatch* parent_face_map_patch = 
        static_cast<FaceMapPatch*>(face_map->patch(face_map_id));

    FaceMapSubPatch* patch = 
        new FaceMapSubPatch(
                parent_face_map_patch, 
                domain.first, 
                domain.second, 
                level,
                0,
                0,
                face_map->_quadrature_weights);
    patch->_quadrant = quad; 
    patch->_quad_id_within_p4est_tree = 0;
    r->patch = patch;
}

void initialize_p4est_leaves_with_subpatches(p4est_t* p4est, PatchSurfFaceMap* face_map){
    // TODO rewrite as a loop over tree and quads to remove the counter variable
    Storage s;
    s.quad_counter = 0;
    s.face_map =face_map;

    p4est_iterate(p4est,
            NULL,
            &s,
            reinitialize_subpatches,
            NULL,
            NULL);
    update_quad_indicies(p4est);
}
void update_quad_indicies(p4est_t* p4est){
    int global_quad_id = 0;
    for(size_t tree_id = 0; tree_id < p4est->trees->elem_count; tree_id++){
        p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, tree_id);

        for(size_t quad_id = 0; quad_id < tree->quadrants.elem_count; quad_id++){
            p4est_quadrant_t* quad = p4est_quadrant_array_index(&tree->quadrants, quad_id);
            //auto refinement_data = static_cast<RefinementData<FaceMapSubPatch>*>(quad->p.user_data);
            auto refinement_data = get_refinement_data(quad);
            FaceMapSubPatch* patch = refinement_data->patch;
            assert(patch != NULL);
            assert(patch->_quadrant != NULL);
            patch->_quadrant = quad;
            patch->_parent_id = tree_id;
            patch->_quad_id_within_p4est_tree = quad_id;
            patch->_id = global_quad_id++;
        }
    }
}

// to be called after patch admissibility refinement.
void set_coarse_patch_ids(p4est_t* p4est){
    int global_quad_id = 0;
    for(size_t tree_id = 0; tree_id < p4est->trees->elem_count; tree_id++){
        p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, tree_id);

        for(size_t quad_id = 0; quad_id < tree->quadrants.elem_count; quad_id++){
            p4est_quadrant_t* quad = p4est_quadrant_array_index(&tree->quadrants, quad_id);
            auto refinement_data = get_refinement_data(quad);
            FaceMapSubPatch* patch = refinement_data->patch;
            patch->_coarse_parent_patch = patch->_id;
        }

    }
}

void check_patches(p4est_t* p4est){
    vector<FaceMapSubPatch*> subpatches = collect_face_map_subpatches(p4est);
    for(auto patch : subpatches){
        assert(patch != NULL);
        assert(dynamic_cast<FaceMapSubPatch*>(patch));
    }

}

vector<Patch*> p4est_to_face_map_subpatches(p4est_t* p4est, 
        PatchSurfFaceMap* face_map){
    cout << "checking patches" << endl;
    check_patches(p4est);

    vector<Patch*> subpatches;

    vector<FaceMapSubPatch*> subpatches_to_convert = collect_face_map_subpatches(p4est);
    
    subpatches.assign(subpatches_to_convert.begin(), subpatches_to_convert.end());
    
    for(int i =0; i < subpatches_to_convert.size(); i++){
        assert(subpatches_to_convert[i] != NULL);
        assert(subpatches[i] != NULL);
        assert(dynamic_cast<FaceMapSubPatch*>(subpatches_to_convert[i]));
        assert(dynamic_cast<FaceMapSubPatch*>(subpatches[i]));
    }

    
    return subpatches;
}



void update_face_map(p4est_t* p4est, 
        PatchSurfFaceMap* face_map,
        PatchSurfFaceMap*& intermediate_face_map,
        PatchSamples*& intermediate_patch_samples){
    if(intermediate_patch_samples != NULL)
        delete intermediate_patch_samples;
    auto temp_face_map = intermediate_face_map;


    
    // Construct face-map patches from current p4est quads and the original face-map, 
    // and assign them to the intermediate face-map 
    intermediate_face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
    intermediate_face_map->_surface_type = face_map->_surface_type;
    intermediate_face_map->setFromOptions();
    intermediate_face_map->setup_from_existing_face_map(face_map);
    intermediate_face_map->initialize_with_existing_p4est(p4est);

    // Re-sample the new set of patches
    vector<int> partition(intermediate_face_map->num_patches(), 0);
    intermediate_patch_samples = new PatchSamples("", "");
    intermediate_patch_samples->bdry() = intermediate_face_map;
    intermediate_patch_samples->patch_partition() = partition;
    intermediate_patch_samples->setup();

}

void resample_qbkix_points(p4est_t* p4est,
        PatchSamples* patch_samples,
        Vec& qbkix_points, 
        vector<int> qbkix_indices){
    // generate qbkix targets from the new samples
    vector<int> invalid_patches = face_map_patch_ids_to_refine(p4est);
    qbkix_points = 
        patch_samples->generate_qbkix_points_from_sample_points(
                qbkix_indices, invalid_patches);


}

void update_and_resample_face_map(p4est_t* p4est, 
        PatchSurfFaceMap* face_map,
        PatchSurfFaceMap*& intermediate_face_map,
        PatchSamples*& intermediate_patch_samples,
        Vec& qbkix_points, 
        vector<int> qbkix_indices){
     
        // EVERYTHING BREAKS IF THIS DELETE IS CALLED
        // TODO BUG SCARY FIGURE OUT WHY
    /*if(intermediate_face_map != NULL)
        delete intermediate_face_map;*/
        update_face_map(p4est, face_map, 
                intermediate_face_map, intermediate_patch_samples);
        resample_qbkix_points(p4est,
                intermediate_patch_samples, qbkix_points, qbkix_indices);

}





// TODO kill and replace with functor or code from <p4est_interface.hpp>
vector<pair<int, int> > patches_that_need_refinement(p4est_t* p4est){
    vector<pair<int, int> > patches_to_refine;
    for(size_t i = 0; i < p4est->trees->elem_count; i++){
        p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, i);

        for(size_t j = 0; j < tree->quadrants.elem_count; j++){
            p4est_quadrant_t* quad = p4est_quadrant_array_index(&tree->quadrants, j);
            auto refinement_data = get_refinement_data(quad);
            if(refinement_data->refine){
                patches_to_refine.push_back(pair<int, int>(i,j));
            }
        }
    }
    return patches_to_refine;
}


vector<FaceMapSubPatch*> face_map_patches_to_refine(p4est_t* p4est){
    vector<pair<int, int> > patches_to_refine = patches_that_need_refinement(p4est);
    
    vector<FaceMapSubPatch*> face_map_patches;
    for(size_t pi = 0; pi < patches_to_refine.size(); pi++){
        auto refinement_data = 
            get_refinement_data_from_forest(p4est, patches_to_refine[pi]);

        assert(refinement_data->refine);
        FaceMapSubPatch* patch = refinement_data->patch;
        face_map_patches.push_back(patch);
    }
    return face_map_patches;
}

vector<int> face_map_patch_ids_to_refine(p4est_t* p4est){
    vector<FaceMapSubPatch*> face_map_patches = face_map_patches_to_refine(p4est);
    
    vector<int> face_map_patch_ids;
    face_map_patch_ids.reserve(face_map_patches.size());
    
    for(const auto& patch : face_map_patches){
        face_map_patch_ids.push_back(patch->V());
        
    }
    return face_map_patch_ids;
}



void store_p4est(string filename, p4est_t* p4est){
    p4est_save(filename.c_str(), p4est, 1);
}

void load_p4est(string filename, 
        MPI_Comm comm,
        //size_t data_size
        p4est_t*& p4est,
        p4est_connectivity_t** connectivity){

    p4est = p4est_load (filename.c_str(), 
            comm,
            sizeof(RefinementData<FaceMapSubPatch>),
            0,
            NULL,
            connectivity);
    p4est_reset_data(p4est, 
            sizeof(RefinementData<FaceMapSubPatch>),
            initialize_user_pointer,
            NULL);

}

// NOTE this is copied and pasted from p4est in order to not complete mess
// something up. memcpy's of user-data variables are removed.
p4est_t* object_safe_p4est_copy(p4est_t* input){
  int copy_data = 1; // always copy user_data
  const p4est_topidx_t num_trees = input->connectivity->num_trees;
  const p4est_topidx_t first_tree = input->first_local_tree;
  const p4est_topidx_t last_tree = input->last_local_tree;
  size_t              icount;
  size_t              zz;
  p4est_topidx_t      jt;
  p4est_t            *p4est;
  p4est_tree_t       *itree, *ptree;
  p4est_quadrant_t   *iq, *pq;
  sc_array_t         *iquadrants, *pquadrants;

  /* create a shallow copy and zero out dependent fields */
  p4est = P4EST_ALLOC (p4est_t, 1);
  memcpy (p4est, input, sizeof (p4est_t));
  p4est->global_first_quadrant = NULL;
  p4est->global_first_position = NULL;
  p4est->trees = NULL;
  p4est->user_data_pool = NULL;
  p4est->quadrant_pool = NULL;

  /* allocate a user data pool if necessary and a quadrant pool */
  if (copy_data && p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->data_size = 0;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* copy quadrants for each tree */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    itree = p4est_tree_array_index (input->trees, jt);
    ptree = p4est_tree_array_index (p4est->trees, jt);
    memcpy (ptree, itree, sizeof (p4est_tree_t));
    sc_array_init (&ptree->quadrants, sizeof (p4est_quadrant_t));
  }
  for (jt = first_tree; jt <= last_tree; ++jt) {
    itree = p4est_tree_array_index (input->trees, jt);
    iquadrants = &itree->quadrants;
    icount = iquadrants->elem_count;
    ptree = p4est_tree_array_index (p4est->trees, jt);
    pquadrants = &ptree->quadrants;
    sc_array_resize (pquadrants, icount);
    memcpy (pquadrants->array, iquadrants->array,
            icount * sizeof (p4est_quadrant_t));
    if (p4est->data_size > 0) {
      P4EST_ASSERT (copy_data);
      for (zz = 0; zz < icount; ++zz) {
        iq = p4est_quadrant_array_index (iquadrants, zz);
        pq = p4est_quadrant_array_index (pquadrants, zz);
        pq->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
        pq->p.user_data =  iq->p.user_data;
        //auto pr = static_cast<RefinementData<FaceMapSubPatch>*>(pq->p.user_data);
        //auto ir = static_cast<RefinementData<FaceMapSubPatch>*>(iq->p.user_data);
        //memcpy (pq->p.user_data, iq->p.user_data, p4est->data_size);
      }
    }
  }

  /* allocate and copy global quadrant count */
  p4est->global_first_quadrant =
    P4EST_ALLOC (p4est_gloidx_t, p4est->mpisize + 1);
  memcpy (p4est->global_first_quadrant, input->global_first_quadrant,
          (p4est->mpisize + 1) * sizeof (p4est_gloidx_t));

  /* allocate and copy global partition information */
  p4est->global_first_position = P4EST_ALLOC (p4est_quadrant_t,
                                              p4est->mpisize + 1);
  memcpy (p4est->global_first_position, input->global_first_position,
          (p4est->mpisize + 1) * sizeof (p4est_quadrant_t));

  /* check for valid p4est and return */
  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
}




END_EBI_NAMESPACE
