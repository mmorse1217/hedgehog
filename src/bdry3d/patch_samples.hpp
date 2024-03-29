#ifndef _DN3DOV_HPP_
#define _DN3DOV_HPP_

#include "patch_surf.hpp"
#include "on_surface_point.hpp"
#include <bdsurf.hpp>
#include <vector>

BEGIN_EBI_NAMESPACE

//----------------------------------------------------------
// PatchSamples handles the regular sampling of the blended surface 
// In practice, two instances are need in the solver:
// one at the target spacing, and another at a refined spacing
// with refinement determined by _refinment_factor (default 4) 
//----------------------------------------------------------
class PatchSurfFaceMap;
class PatchSamples: public EbiObject
{
public:
    enum {	 
        EVAL_VL = BdSurf::EVAL_VALUE,
        EVAL_FD = BdSurf::EVAL_1ST_DERIV,
        EVAL_2ND_DERIV = BdSurf::EVAL_2ND_DERIV 
    };
  struct Tag {
	int _gid; //group id == index of multiply connected boundary components of 
              //domain 
	bool _dmt; //dominant or not; i.e. whether a sample points face
		   //coordinates wrt its patch center are < 0.5 
  };
    /** 
     * Needed for blendsurf support
     * _patch_sampling_index:(sqrt(num samples)) x (sqrt(num samples)) matrix
     * that serves as an indexing scheme of sample points in the bounding
     * box of the ith patch. Sample points **outside** the patch are labeled
     * with a -1; sample points that fall within the patch are labeled with a 
     * unique integer. Continuing with our example, the patch would be labelled
     * as such (pardon the improper spacing):
     *
     * ______________________________________________
     * |-1 -1 -1 -1 -1 -1     / \ -1 -1 -1 -1 -1 -1  |
     * |-1 -1 -1 -1 -1      / 1  2\   -1 -1 -1 -1 -1 |
     * |------------------/3  4  5  6\---------------|
     * |  \   7  8  9 10 11 12 13 14 15 16 17 18 19 /|
     * |-1 \ 20 21 22 23 24 25 26 27 28 29 30 31 32/-1|
     * |-1 / 33 34 35 36    ^ 37 38 39 40 41 42 \ -1 |
     * |  / 43 44 45   / -1 -1 -1 \ 46 47 48 49  \   |
     * ______________________________________________
     *
     * Although the image above is rather terrifying, it should convey the general
     * idea. Sample points that land within the domain of the ith patch are
     * labelled in consecutive order from 1 to num_sample_points_in_patch[i], 
     * and all the rest are labeled as -1. This is used to select and identify
     * relavent sample points for computation.
     */
    vector<IntNumMat> _patch_sampling_index;
protected:
    PatchSurf* _bdry;

    //!< max 3d distance spacing of samples
    double _spacing;              

    double _boundary_distance_ratio;
    double _interpolation_spacing_ratio;

    //!<  Number of Lagrange interpolation points
    //int _surface_interpolation_num_samples;                 

    // ---------------------  Per patch data  ---------------------
    // Vectors below are of size # of patches and vec[i] indicates
    // some property or value corresponding to the ith patch.
    // all stored indices are global unless noted otherwise 

    //!<  patch_partition[i] = partition patch is assigned to
    vector<int> _patch_partition;

    //!<  _step_size[i] is the spacing between regular sample points in the
    //! i-th patch domain
    vector<double> _step_size;

    /** _num_sample_points: the square root of the 
     * number of total sample points taken in the bounding box around the ith
     * patch. In practice, this number is 2*(half-width of bounding box around
     * the ith patch)*(jacobian of ith patch)/_spacing. This looks roughly like
     * (of course with the bottom two points actually meeting :) )
     *
     * ________________
     * |....../\......|
     * |...../..\.....|
     * |----/....\----|
     * |\............/|
     * |.\........../.|
     * |./....^.....\.|
     * |/.../...\....\|
     * ________________
     *
     * The sample points occur at EACH ascii character, ".", "-", "\", "/",
     * "|", and "^", and the space between each character cooresponds to
     * _step_size[i] for the ith patch. _num_sample_points[i] is the total
     * number of these sample points over the entire bounding box (in the above
     * example, it's the total number of all these ascii characters), squared.
     * Or another way to to think of it is that _num_sample_points[i] is the
     * number of samples in a horizontal slice of the boundaing box.
     */
    vector<int>    _num_sample_points;
                                                            
    //. how to order samples
    
    /** _num_sample_points_in_patch: number of sample points in the bounding box
     * around the ith patch that ***fall within the bounds of the patch.***
     * This is different from _num_sample_points. The picture instead looks
     * like:
     *
     * ________________
     * |      /\      |
     * |     /..\     |
     * |----/....\----|
     * |\............/|
     * | \........../ |
     * | /....^.....\ |
     * |/.../   \....\|
     * ________________
     *
     * The sample_points that reside within the patch are only the ascii
     * characters ".", "\", "/", and "^", still occuring at spacing
     * _step_size[i] for the ith patch. So, _num_sample_points_in_patch is the
     * total number of samples inside the patch (in the above example, it's the total
     * number of ".", "\", "/", and "^" that occur).
     *
     * Note that _num_sample_points_in_patch[i] <= _num_sample_points[i], and
     * generally the inequality is strict.
     */
    vector<int> _num_sample_points_in_patch;
    
    //! the starting index of the unform sample points on the ith patch
    vector<int> _sample_point_starting_index;
  
  

  // ---------------------  per sample data   ---------------------
  //  Stored in distributed Petsc vectors; 
  /** Patch points represented by (face_idx, x,y) where x,y is coordinates 
   *  relative to the face; face_idx stored as double
   */
  Vec _sample_as_face_point;
  /**  3d positions of sample points  
   */
  Vec _sample_point_3d_position;

  //!<  Unit normal vector to the surface at patch coordinate of the sample 
  //! point
  Vec _sample_point_normal;

  //!<  Jacobian of the surface at the sample point
  Vec _sample_point_jacobian;
  //!<  Patch blending function value of the patch at the sample point 
  Vec _sample_point_blend_func_value;

  //!<  Quadrature weight for the ith sample point. Equal for all samples of 
  //! a patch i  to h_i^2 (trapezoid rule)
  Vec _sample_point_quad_weight;

  //!<  Vector of Tag structs indicating the group id and whether a sample 
  //! point is
  Vec _sample_point_props;
  Vec _sample_point_far_field;
  Vec _sample_point_interpolant_spacing;
  
  Vec _sample_point_combined_weight;  

  Vec _sample_point_parametric_preimage;
  NumVec <OnSurfacePoint> _sample_point_as_on_surface_point;
  /** Bottom of page 9 in Ying et. al. A High-order 3d BIE...
   *  wcb[k] = jacobian[k] * weight[k] * alpha[k]
   * 
   * _sample_point_combined_weight[i] is the w_k(g_k(c_k))*J_k(c_k)*dc_k in the
   * intergral at the bottom of page 8. This incorporates the weight funciton of 
   * the  partition of unity of patches, the change of variables with respect
   * to the surface parametrization (i.e. patch representation) and the
   * quadrature weight (dc_k = h_i^2 for trapezoidal rule).
   */

  // Note the following vectors have the same relative order, i.e.
  // patches[i] has center boundary_component_center[i] and orientation ortvec[i]
  //LEXING: DOUBLE CHECK THE FOLLOWING  vector<DblNumMat> _refined_datvec;
  //temporary data  Vec _sample_point_data; //used to store temporary data  int _dof;
  DblNumVec _barycentric_weights_x;
  DblNumVec _barycentric_weights_y;
  DblNumVec _interpolation_nodes_x;
  DblNumVec _interpolation_nodes_y;
public:
    PatchSamples(const string& n, const string& p);
    PatchSamples(PatchSurf* surface);
    ~PatchSamples();  
    int setup(bool refined=false); 
  

    // A thrilling number of accessors
    PatchSurf*& bdry() {                 return _bdry; }
    double& spacing() {               return _spacing; }
    //int& surface_interpolation_num_samples() {
                                      //return _surface_interpolation_num_samples; }

    vector<int>& patch_partition() {  return _patch_partition; }
    vector<double>& step_size() {     return _step_size; }
    vector<int>& num_sample_points(){ return _num_sample_points; }
    vector<int>& num_sample_points_in_patch() { 
                                      return _num_sample_points_in_patch; }
    vector<int>& sample_point_starting_index() {        
                                      return _sample_point_starting_index; }
    vector<IntNumMat>& patch_sampling_index() {
                                      return _patch_sampling_index; }
    Vec sample_as_face_point() {      return _sample_as_face_point; }
    Vec sample_point_3d_position() {     return _sample_point_3d_position;  }
    Vec sample_point_normal() {       return _sample_point_normal;    }
    Vec sample_point_jacobian() {     return _sample_point_jacobian;  }
    Vec sample_point_blend_func_value() {  return _sample_point_blend_func_value; }
    Vec sample_point_quad_weight() {       return _sample_point_quad_weight;    }
    Vec sample_point_combined_weight() {  return _sample_point_combined_weight; }
    Vec sample_point_props() {        return _sample_point_props; }
    Vec sample_point_parametric_preimage(){ return _sample_point_parametric_preimage; }
    Vec sample_point_far_field() {        return _sample_point_far_field; }
    Vec sample_point_interpolant_spacing() {        return _sample_point_interpolant_spacing; }
    NumVec<OnSurfacePoint> sample_point_as_on_surface_point(){ return _sample_point_as_on_surface_point; }
    DblNumVec& barycentric_weights_x(){
        return _barycentric_weights_x;
    }
    DblNumVec& barycentric_weights_y(){
        return _barycentric_weights_y;
    }
    DblNumVec& interpolation_nodes_x(){
        return _interpolation_nodes_x;
    }
    DblNumVec& interpolation_nodes_y(){
        return _interpolation_nodes_y;
    }


    //! Collocation point on a patch P is a sample point of this patch, or a
    //! sample point of any other patch sharing a face with P
    //! Construct:
    //! colloc_point_global_indices    - the list of collocation points for the current partition
    //! num_colloc_point_in_patch - number of collocation points per patch 
    //! collocation_points_starting_indexvec - offset of the list of collocation points for the
    //! current patch in the global Vec of collocation points (colloc_point_as_face_point,
    //! colloc_point_3d_position stored in CollocationPointData, attached to Be3dOv derived classes 
    //! target_as_face_point is used to get the list of sample points for this partition 
    int refine_data(int dof, int refinement_factor, Vec dat, vector<DblNumMat>& refined_datvec); 
    Vec refine_density(int dof, int refinement_factor, Vec density, PatchSamples* refined_patch_samples); 

    //! Check whether (x,y) is a contained in the (pi)th chart. See
    //! BlendedPatch.is_xy_valid in bd3dbd.hpp for more.
    int is_interp_xy_valid(int pi, double* xy, bool& is_valid);

    //! Perform a change of coordinates from (x,y) \in D of Figure 1 in
    //! blendsurf paper to cd coordinates \in S of Figure 1, for an (x,y) \in D.
    int interpolated_position_and_derivatives(int pi, double* xy, int flag, double* res);
    
    int interpolate_data(int pi, double* xy, int refinement_factor, 
        int surface_interpolation_num_samples,
        int flag,
        vector<DblNumMat>& refined_datvec, 
        double* res); 
    
    //! Check whether the (i,j)th sample point is contained within the chart
    int is_sample_point_valid(int pi, int* ij, bool& is_valid);

    //! Get position, normal vector and jacobian of the (i,j)th sample point on 
    //! the (pi)th patch. See BlendedPatch.xy_to_patch_coords in bd3dbd.hpp for more.
    int get_sample_on_patch(int pi, int* ij, double* pos, double* nor, 
        double* jac);

    int get_sample_point(int pi, int* ij, int dof, Vec dat, double* res);

    int dim() { return 3; }
    int face_point_size_in_doubles() { return _bdry->face_point_size_in_doubles(); }  

    /**
     * The following methods all provide accessors into the cooresponding
     * arrays above, selecting out set of sample data for the pi-th patch 
     */
    DblNumMat sample_as_face_point(int pi);
    DblNumMat sample_point_3d_position(int pi);
    DblNumMat sample_point_normal(int pi);
    DblNumVec sample_point_jacobian(int pi);
    DblNumVec sample_point_blend_func_value(int pi);
    DblNumVec sample_point_quad_weight(int pi);
    DblNumVec sample_point_combined_weight(int pi);
    DblNumVec sample_point_props(int pi);
    DblNumMat sample_point_parametric_preimage(int pi);
    DblNumMat sample_point_data(int pi, int dof, Vec dat);
    DblNumVec sample_point_far_field(int pi);
    DblNumVec sample_point_interpolant_spacing(int pi);
    /**
     * Starting index of sample points for the (pi)th patch relative to starting 
     * index of partition 
     */
    int sample_point_local_start_index(int pi); 
      
    Vec generate_qbkix_points_from_sample_points(vector<int> qbkix_indices,
            vector<int> patches_to_sample);
  //simplify the following functions
  int  local_num_sample_points() {
      int64_t ret;
      VecGetLocalSize(_sample_point_3d_position,&ret);
      return ret/dim();
  }
  int  global_num_sample_points() {
      int64_t ret;
      VecGetSize(     _sample_point_3d_position,&ret);
      return ret/dim();
  }
  void local_sample_point_range(int& a, int& b) {
      int64_t tmp_a, tmp_b;
      VecGetOwnershipRange(_sample_point_3d_position, &tmp_a, &tmp_b); 
      a=tmp_a/dim(); 
      b=tmp_b/dim();
  }  
NumVec <OnSurfacePoint> closest_points_to_qbkix_points();
private:
    void initialize_sampling_vectors(const int num_patches);
    void initialize_parallel_vectors(const int local_num_sample_points);
    int initialize_sampling_indices();
    void sample_patches();
    void set_equal_sample_rate_param_space(const vector<Patch*>& patches);
    void set_equal_sample_rate_physical_space(const vector<Patch*>& patches);
};

Vec refine_function(int dof, Vec function,
        PatchSamples* coarse_samples,
        PatchSamples* fine_samples);

Vec interpolate_and_resample(int dof, Vec function,
        PatchSamples* interp_samples,
        PatchSamples* eval_samples);
void generate_samples_on_child_patches(MPI_Comm comm, 
        PatchSurfFaceMap* face_map,
        int num_samples_1d,
        Vec& uv_coordinates_single_patch, 
        Vec& sample_points_on_chlildren,
        vector<int64_t> patches_to_sample=vector<int64_t>());
END_EBI_NAMESPACE

#endif

