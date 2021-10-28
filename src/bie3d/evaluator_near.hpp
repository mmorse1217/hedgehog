#ifndef _BE3DOVNEA_HPP_
#define _BE3DOVNEA_HPP_

#include "evaluator.hpp"
#include "evaluator_on_surface.hpp"
#include "collocation_patch_samples.hpp"
BEGIN_EBI_NAMESPACE

enum ExpansionType {
    INTERPOLATE = 0,
    EXTRAPOLATE_ONE_SIDE = 1,
    EXTRAPOLATE_TWO_SIDE = 2,
    INTERPOLATE_ACROSS_SURFACE = 3,
    EXTRAPOLATE_ONE_SIDE_CHEBYSHEV = 4
};

class EvaluatorNear: public Evaluator
{

    private:
        bool _fmm_initialize;
    protected:
        // These are set at the initalization time by the caller
        // i.e. point marking happens 
        //!< 3d positions of targets  in Omega_2 (near targets)
        Vec _target_3d_position;
        //!< 3d positions of targets  in Omega_1 (intermediate targets)
        Vec _intermediate_target_3d_position;

        //!< in==1 or out==-1
        Vec _target_in_out;
        //!< 3d positions of closest sample points on the surface for each
        //!  target point 
        Vec _closest_sample_3d_position;
        //!< face points corresponding to closest sample points 
        Vec _closest_sample_as_face_point;

        //!< determines the size of Omega_1 zone (where values are computed by interpolation _alfcoef*spacing()
        double _alfcoef;
        //!< vector of steps for Lagrangian interpolation
        vector<double> _interpolation_nodes;
    
        int64_t _surface_interpolation_num_samples;
        int64_t _refinement_factor;
        //  Additional vectors, computed in the setup phase
        //!< distancs to the closest points on the boundary, length(_target_3d_position-_closest_sample_3d_position)  
        Vec _distance_closest_to_sample;

        //!< which of two regions a point belongs to (1 = Omega_2 + on the surface,  2 = Omega_1 (confusing)  paper p. 260  

        //!< target points in Omega_1 zone: original targets , or generated for 
        //   interpolation
        Vec _target_positions_intermediate; 

        //!< 3d positions of closest boundry points to targets in Omega_2 zone  
        Vec _closest_samples_3d_position_near;

        //!< face points of closest boundary points to targets in Omega_2 zone 
        Vec _closest_samples_as_face_point_near; 

        //KnlMat3d* _knlmatfar; //old
        // data structure for keeping track of
        CollocationPatchSamples* _refined_collocation_data;  
        CollocationPatchSamples* _collocation_data;  

        // a copy of on-the-surface evaluator (distinct
        // from the one used for the solver, initialized with positions/face ponts of closest samples to Omega_2 zone points
        EvaluatorOnSurface* _on_surface_evaluator;

        //TEMPORARY DATA
        //vector<DblNumMat> _refined_datvec;
        PatchSamples* _refined_patch_samples;
    public:
        FMM* fmm;
        EvaluatorNear(const string& n, const string& p):
            Evaluator(n,p), _target_3d_position(NULL), _target_in_out(NULL),
            _closest_sample_3d_position(NULL), _closest_sample_as_face_point(NULL),
            _alfcoef(1), _distance_closest_to_sample(NULL),
            _target_positions_intermediate(NULL), _closest_samples_3d_position_near(NULL),
            _closest_samples_as_face_point_near(NULL), _on_surface_evaluator(NULL), fmm(NULL)
    { 
        _refined_collocation_data = new CollocationPatchSamples();
        _collocation_data = new CollocationPatchSamples();
        _fmm_initialize = true;
    }
        ~EvaluatorNear() {
            if(_distance_closest_to_sample!=NULL) {
                VecDestroy(&_distance_closest_to_sample);
            }
            if(_target_positions_intermediate!=NULL) {
                VecDestroy(&_target_positions_intermediate);
            }
            if(_closest_samples_3d_position_near!=NULL) { 
                VecDestroy(&_closest_samples_3d_position_near);
            }
            if(_closest_samples_as_face_point_near!=NULL) {
                VecDestroy(&_closest_samples_as_face_point_near);
            }
            if(fmm!=NULL) {
                delete fmm;
            }
            if(_on_surface_evaluator!=NULL) {
                delete _on_surface_evaluator;
            }
            delete _refined_collocation_data;
            delete _collocation_data;
        }
        int setFromOptions();
        int setup();
        int eval(Vec den, Vec val);
        void set_refined_surface_discretization(PatchSamples* patch_samples){
            _refined_patch_samples = patch_samples;
        }
        //accessors:
        Vec& target_3d_position() {               return _target_3d_position; } 
        Vec& target_in_out() {                    return _target_in_out; }
        Vec& closest_sample_3d_position() {       return _closest_sample_3d_position; }
        Vec& closest_sample_as_face_point() {     return _closest_sample_as_face_point; }
        Vec& intermediate_target_3d_position(){   return _intermediate_target_3d_position; }
        double& alfcoef() {                       return _alfcoef; }
        vector<double>& interpolation_nodes() {
            return _interpolation_nodes;
        }
        void set_fmm_init(bool b){
            _fmm_initialize = b;
        }
    protected:

};
// returns the weights (in the last argument) of a polynomial interpolant 
// at u (3rd arg), of degree n-1 with data at p[i] (second arg)
int compute_interpolation_weights(int n, double*, double, double*); 
void jump_evaluation(Vec den, Vec val, vector<DblNumMat>& _refined_datvec, 
        Kernel3d kernel, PatchSamples* patch_samples, 
        CollocationPatchSamples* _collocation_data);

Vec interpolate_density_to_refined_grid(PatchSamples* patch_samples, 
        Vec& refined_sample_point_3d_position,
        Vec& refined_sample_point_combined_weight,
        vector<DblNumMat>& _refined_datvec,
        CollocationPatchSamples* _refined_collocation_data,
        int source_dof);

void interpolation_near(Vec& intermediate_targets, Vec& closest_surface_targets,
        Vec& near_targets, Vec& intermediate_potentials, 
        Vec& closest_surface_potential, Vec& jump_potential, 
        vector<DblNumMat>& _refined_datvec, Vec& near_target_in_out,
        Vec& val, int target_dof, ExpansionType expansion_type);

vector<double> get_interpolation_nodes(ExpansionType expansion_type);
vector<double> get_chebyshev_nodes();

Vec generate_auxiliary_interpolation_points(Vec near_targets, 
        Vec closest_surface_points_near, Kernel3d kernel,  
        ExpansionType expansion_type=INTERPOLATE,
        Vec interpolation_directions=NULL, 
        Vec sample_point_far_field=NULL,
        Vec sample_point_interpolant_spacing=NULL);

Vec generate_interior_qbkix_points(Vec near_targets, 
        Vec closest_surface_points_near, 
        vector<int> qbkix_point_index,
        Vec interpolation_directions=NULL, 
        Vec sample_point_far_field=NULL,
        Vec sample_point_interpolant_spacing=NULL);

END_EBI_NAMESPACE

#endif
