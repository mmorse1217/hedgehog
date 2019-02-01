#ifndef __COLLOCATION_PATCH_SAMPLES_HPP__
#define __COLLOCATION_PATCH_SAMPLES_HPP__

#include "solver_utils.hpp"
#include "common/ebi.hpp"
#include "common/numvec.hpp"
#include "evaluator.hpp"
#include "bdry3d/patch_samples.hpp"
#include "common/ebiobject.hpp"

BEGIN_EBI_NAMESPACE
class CollocationPatchSamples: public EbiObject{
    protected:

        //!< the list of collocation points for the current partition, 
        //  all patches together 
        vector<int> _colloc_point_global_indices;

        //!< total number of collocation points for each
        //!patch in the current partition (i.e., its own
        //!sample points, and points of overlapping
        //!patches that share a face with this patch
        vector<int> _num_colloc_point_in_patch;

        //! offsets of the list of collocation points for a patch in the global 
        //! vector of collocation points
        vector<int> _colloc_point_starting_index;

        Vec _colloc_point_3d_position;
        Vec _colloc_point_as_face_point;

        VecScatter _colloc_to_target_density_scatter;
        VecScatter _colloc_to_target_value_scatter;
        
        int get_patch_offset(int pi);
        //int get_patch_offset(int pi);
    public:
        CollocationPatchSamples(): EbiObject("", ""),
                                _colloc_point_3d_position(NULL), 
                                _colloc_point_as_face_point(NULL),
                                _colloc_to_target_density_scatter(NULL),
                                _colloc_to_target_value_scatter(NULL)
                                 {;}

        ~CollocationPatchSamples() { clear(); }

        void clear() {
            if(_colloc_point_3d_position!=NULL) {
                VecDestroy(&_colloc_point_3d_position);
                _colloc_point_3d_position=NULL;
            }
            if(_colloc_point_as_face_point!=NULL) {
                VecDestroy(&_colloc_point_as_face_point);
                _colloc_point_as_face_point=NULL;
            }
            if(_colloc_to_target_density_scatter!=NULL) { 
                VecScatterDestroy(&_colloc_to_target_density_scatter);
                _colloc_to_target_density_scatter=NULL;
            }
            if(_colloc_to_target_value_scatter!=NULL) {
                VecScatterDestroy(&_colloc_to_target_value_scatter);
                _colloc_to_target_value_scatter=NULL;
            }
        }

        // The meat and potatoes
        int distribute_collocation_points(Vec target_3d_position, 
                                          Vec target_as_face_point,
                                          PatchSamples* patch_samples, int sdof, int tdof);

        //! Collocation point on a patch P is a sample point of this patch, or a
        //! sample point of any other patch sharing a face with P
        //! Construct:
        //! colloc_point_global_indices    - the list of collocation points for
        //the current partition
        //! num_colloc_point_in_patch - number of collocation points per patch 
        //! collocation_points_starting_indexvec - offset of the list of 
        //collocation points for the
        //! current patch in the global Vec of collocation points 
        //(colloc_point_as_face_point,
        //! colloc_point_3d_position stored in CollocationPointData, attached 
        //to Be3dOv derived classes 
        //! target_as_face_point is used to get the list of sample points for 
        //this partition 

        int build_local_collocation_lists(Vec target_as_face_point, 
                vector<int>& colloc_point_global_indices,
                vector<int>& num_colloc_point_in_patch,
                vector<int>& collocation_points_starting_indexvec,
                PatchSamples* patch_samples);

        int num_collocation_points(int pi);  

        int collocation_points_starting_index(int pi); 

        DblNumMat colloc_point_3d_position(int pi);

        DblNumMat colloc_point_as_face_point(int pi, int face_point_size_in_doubles);

        DblNumMat coldat(int pi, Vec dat, int dof); 


        vector<int>& colloc_point_global_indices() {
            return _colloc_point_global_indices;
        }
        vector<int>& num_colloc_point_in_patch() {
            return _num_colloc_point_in_patch;
        }
        vector<int>& colloc_point_starting_index() {
            return _colloc_point_starting_index;
        }
        Vec& colloc_point_3d_position() {
            return _colloc_point_3d_position;
        }
        Vec& colloc_point_as_face_point() {
            return _colloc_point_as_face_point;
        }
        VecScatter& colloc_to_target_density_scatter() {
            return _colloc_to_target_density_scatter;
        }
        VecScatter& colloc_to_target_value_scatter() {
            return _colloc_to_target_value_scatter;
        }
};

END_EBI_NAMESPACE
#endif
