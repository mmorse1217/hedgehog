#include "patch_samples.hpp"
#include "common/vecmatop.hpp"
#include <iostream>
#include "patch_surf_blended.hpp"
#include "patch_surf_analytic.hpp"
#include "bie3d/evaluator_near.hpp"
#include "face_map_subpatch.hpp"
#include "common/interpolate.hpp"
#include "common/stats.hpp"
#include <unordered_map>
#include <sampling.hpp>
#include <profile.hpp>
BEGIN_EBI_NAMESPACE

using std::min;
using std::max;
using std::abs;
using std::cerr;
using Sampling::sample_1d;
using Sampling::sample_2d;
using Sampling::equispaced;
using Sampling::chebyshev1;
using Sampling::chebyshev2;


PatchSamples::PatchSamples(const string& n, const string& p):
    EbiObject(n,p), 
    _sample_as_face_point(NULL),
    _sample_point_3d_position(NULL),
    _sample_point_normal(NULL),
    _sample_point_jacobian(NULL),
    _sample_point_blend_func_value(NULL),
    _sample_point_quad_weight(NULL),
    _sample_point_props(NULL),
    _sample_point_far_field(NULL),
    _sample_point_interpolant_spacing(NULL),
    _sample_point_combined_weight(NULL),
    _sample_point_parametric_preimage(NULL)
{
    //----------------------------------------------------------------------  
    // Load options from Petsc
    //----------------------------------------------------------------------  
    _spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    _boundary_distance_ratio = Options::get_double_from_petsc_opts("-boundary_distance_ratio");
    _interpolation_spacing_ratio = Options::get_double_from_petsc_opts("-interpolation_spacing_ratio");
}

PatchSamples::PatchSamples(PatchSurf* surface): PatchSamples("", "") {
    _bdry = surface;
}

PatchSamples::~PatchSamples()
{
  if(_sample_as_face_point!=NULL) { VecDestroy(&_sample_as_face_point); }
  if(_sample_point_3d_position!=NULL) { VecDestroy(&_sample_point_3d_position); }
  if(_sample_point_normal!=NULL) { VecDestroy(&_sample_point_normal); }
  if(_sample_point_jacobian!=NULL) { VecDestroy(&_sample_point_jacobian); }
  if(_sample_point_blend_func_value!=NULL) { VecDestroy(&_sample_point_blend_func_value); }
  if(_sample_point_quad_weight!=NULL) { VecDestroy(&_sample_point_quad_weight); }
  if(_sample_point_combined_weight!=NULL) { VecDestroy(&_sample_point_combined_weight); }
  if(_sample_point_props!=NULL) { VecDestroy(&_sample_point_props); }
  if(_sample_point_far_field!=NULL) { VecDestroy(&_sample_point_far_field); }
  if(_sample_point_interpolant_spacing!=NULL) { VecDestroy(&_sample_point_interpolant_spacing); }
  if(_sample_point_parametric_preimage!=NULL) { VecDestroy(&_sample_point_parametric_preimage); }
}

void PatchSamples::initialize_sampling_vectors(const int num_patches){
    _step_size.resize( num_patches );
    _num_sample_points.resize( num_patches );
    _num_sample_points_in_patch.resize( num_patches, 0 );
    _sample_point_starting_index.resize( num_patches, 0 );
    _patch_sampling_index.resize( num_patches );  //_refined_datvec.resize( num_patches );
}

void PatchSamples::set_equal_sample_rate_param_space( const vector<Patch *> &patches) {
  for (int pi = 0; pi < patches.size(); pi++) {
    Patch *curpch = patches[pi];
    assert(curpch);

    _step_size[pi] = _spacing;
    _num_sample_points[pi] = floor(1. / _spacing) + 1;
  }
}

void PatchSamples::set_equal_sample_rate_physical_space( const vector<Patch *> &patches) {
  for (int pi = 0; pi < patches.size(); pi++) {
    Patch *curpch = patches[pi];
    assert(curpch);
    double bnd = curpch->bnd(); // Maximum value for x/y in parameter space

    double estjac;

    curpch->estimate_jacobian(&estjac);
    // ---------------------------------------------------------------------
    double tmp = _spacing / estjac; // max( ret[1].l2(), ret[2].l2() );
    _num_sample_points[pi] = (int)(ceil(bnd / tmp)) * 2; //

    // distance between regular samples of the bounding box enclosing the patch
    _step_size[pi] = 2.0 * bnd / double(_num_sample_points[pi]);
  }
}

int PatchSamples::initialize_sampling_indices(){
    int local_num_sample_points = 0;
    const int num_patches = _bdry->num_patches();
    DblNumMat parametric_samples;
    if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
        int num_samples = _num_sample_points[0];
        // tensor-product chebyshev for face-map
        parametric_samples = DblNumMat(2, num_samples*num_samples, true, 
                sample_2d<chebyshev1>(num_samples, Sampling::base_domain).data());
    }

    // over all patches assigned to this partition
    for(int pi=0; pi<num_patches; pi++) {
        if(_patch_partition[pi] == mpiRank()) {
            Patch* curpch = _bdry->patch(pi);
            double bnd = curpch->bnd();
            double init = -bnd;
            int num_samples = _num_sample_points[pi];		
            double step = _step_size[pi];

            IntNumMat& patch_sampling_index = _patch_sampling_index[pi];
            patch_sampling_index.resize(num_samples,num_samples);
            setvalue(patch_sampling_index, -1); //initial value is -1
            //reg ord := regular order of grid. index of sampled points in 
            //bounding box contained in the current patch
            int cnt = 0;
            for(int j=0; j<num_samples; j++) {
                for(int i=0; i<num_samples; i++) {
                    Point2 xy;
                    if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
                        xy = Point2(parametric_samples.clmdata(j*num_samples+i));
                    } else {
                        // MJM 1/2018 non-standard spacing for blendsurf points.
                        // leaving this alone.
                        //xy[0] = init + i*step;
                        //xy[1] = init + j*step;
                        xy = init*Point2(1) + Point2(i,j)*step;
                    }
                    bool is_valid;
                    // check if (i,j) sample is inside the star-shaped patch domain
                    curpch->is_xy_valid(xy.array(), is_valid);

                    if(is_valid) {
                        patch_sampling_index(i,j) = cnt++;
                    }
                }
            }
            _num_sample_points_in_patch[pi] = cnt;
            local_num_sample_points += cnt;
        }
    } 
    int cnt = 0;
    for(const auto num_points : _num_sample_points_in_patch)
        cnt += num_points;
    /*if(refined)
    stats.add_result("total refined samples", cnt);
    else
    stats.add_result("total coarse samples", cnt);*/
    int num_local_valid_samples;
    // compute global starting index for this partition
    MPI_Scan( &local_num_sample_points, 
            &num_local_valid_samples, 
            1,
            MPI_INT, 
            MPI_SUM, 
            mpiComm());
    num_local_valid_samples -= local_num_sample_points;

    // number of local valid samples on other processors
    int actual_num_local_valid_samples = num_local_valid_samples;

    if(dynamic_cast<PatchSurfFaceMap*>(_bdry) != NULL){
        // face-map
        for(int i =0; i < _sample_point_starting_index.size(); i++){
            _sample_point_starting_index[i] = actual_num_local_valid_samples + 
                     i*_num_sample_points[i]*_num_sample_points[i];
            _num_sample_points_in_patch[i] = _num_sample_points[i]*_num_sample_points[i];

        }
    } else{
        // blendsurf + analytic
        for(int pi=0; pi<num_patches; pi++) {
            if(_patch_partition[pi]==mpiRank()) {
                _sample_point_starting_index[pi] = num_local_valid_samples;
                num_local_valid_samples += _num_sample_points_in_patch[pi];
            }
        }
        vector<int> tmp_uniform_num_samples_vec(_num_sample_points_in_patch);
        // get number of samples and starting indices of other partitions
        MPI_Allreduce( &(tmp_uniform_num_samples_vec[0]), 
                &(_num_sample_points_in_patch[0]),
                num_patches, 
                MPI_INT, 
                MPI_SUM, 
                mpiComm() );

        vector<int> tmp_sample_point_starting_index(_sample_point_starting_index);

        MPI_Allreduce( &(tmp_sample_point_starting_index[0]),
                &(_sample_point_starting_index[0]), 
                num_patches, 
                MPI_INT,
                MPI_SUM, 
                mpiComm() );
    }
    cout << "making vecs." << endl;
    return local_num_sample_points;
}

void PatchSamples::initialize_parallel_vectors(const int local_num_sample_points){
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points*face_point_size_in_doubles(), _sample_as_face_point);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points*dim(), _sample_point_3d_position);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points*dim(),  _sample_point_normal);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points, _sample_point_jacobian);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points, _sample_point_blend_func_value);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points, _sample_point_quad_weight);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points, _sample_point_combined_weight);  
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points, _sample_point_props);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points, _sample_point_far_field);
    Petsc::create_mpi_vec(mpiComm(), local_num_sample_points, _sample_point_interpolant_spacing);
    Petsc::create_mpi_vec(mpiComm(), 2*local_num_sample_points, _sample_point_parametric_preimage);
}

DblNumMat integrate_tensor_product_lagrange_basis_funcs(const int num_samples) {
  DblNumMat lagrange_basis_integrals(num_samples, num_samples);
  // 1d interpolation nodes
  // assuming chebyshev in each dimension
  DblNumVec nodes(
      num_samples, true,
      sample_1d<chebyshev1>(num_samples, Sampling::base_domain).data());

  int64_t patch_order =
      Options::get_int_from_petsc_opts("-bd3d_facemap_patch_order");

  int quad_order = max(50, int(patch_order + 1));

  // TODO BUG THIS IS GENERATING SAMPLES Y FIRST THEN X AND IS
  // INCONSISTENT WITH REST OF THE CODE
  DblNumVec quadrature_weights_1d(num_samples);
  for (int i = 0; i < num_samples; i++) {
    quadrature_weights_1d(i) = Interpolate::integrate_ith_lagrange_basis_func(
        i, 0., 1., num_samples, nodes, quad_order, 1.);
  }
  for (int i = 0; i < num_samples; i++) {
    for (int j = 0; j < num_samples; j++) {

      lagrange_basis_integrals(i, j) =
          quadrature_weights_1d(i) * quadrature_weights_1d(j);
    }
  }
  return lagrange_basis_integrals;
}

void PatchSamples::sample_patches(){

    double sum = 0.;
    
    DblNumMat lagrange_basis_integrals;
    
    DblNumMat parametric_samples;
    if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
        const int num_samples = _num_sample_points[0]; // all patches have same num samples per patch
        lagrange_basis_integrals = integrate_tensor_product_lagrange_basis_funcs(num_samples);
        parametric_samples = DblNumMat(2, num_samples*num_samples, true, 
                sample_2d<chebyshev1>(num_samples, Sampling::base_domain).data());
    }

    string qbkix_convergence_type = Options::get_string_from_petsc_opts("-qbkix_convergence_type");
    bool qbkix_classical_conv = qbkix_convergence_type == "classic";
    bool qbkix_adaptive_conv = qbkix_convergence_type == "adaptive";
    assert(qbkix_classical_conv || qbkix_adaptive_conv);

#pragma omp parallel for
    for(int pi=0; pi<_bdry->num_patches(); pi++) {
        // for patches assigned to this partition
        if(_patch_partition[pi]==mpiRank()) {
            Patch* curpch = _bdry->patch(pi);
            //auto subpatch = dynamic_cast<FaceMapSubPatch*>(curpch);
            double bnd = curpch->bnd();
            // [-bnd,bnd]^2 for blendsurf
            double init = -bnd;
            
            //[0,1]^2 for face-map
            /*double patch_size;
            if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
                patch_size = curpch->characteristic_length();
                max_patch_size = max(max_patch_size, patch_size);
            }*/

            int num_samples = _num_sample_points[pi];
            double step = _step_size[pi];

            IntNumMat& patch_sampling_index = _patch_sampling_index[pi];

            DblNumMat face_point          = sample_as_face_point(pi);
            DblNumMat position            = sample_point_3d_position(pi);
            DblNumMat normal              = sample_point_normal(pi);
            DblNumMat parametric_preimage = sample_point_parametric_preimage(pi);

            DblNumVec jacobian            = sample_point_jacobian(pi);
            DblNumVec blend_func_value    = sample_point_blend_func_value(pi);
            DblNumVec quad_weight         = sample_point_quad_weight(pi);
            DblNumVec combined_weight     = sample_point_combined_weight(pi);
            DblNumVec properties          = sample_point_props(pi);
            DblNumVec far_field           = sample_point_far_field(pi);
            DblNumVec interpolant_spacing = sample_point_interpolant_spacing(pi);

            //FaceMapSubPatch* patch = dynamic_cast<FaceMapSubPatch*>(_bdry->patches()[pi]);
            //patch->compute_surface_area();
            //double max_sample_spacing;

            // Compute distance between further samples
            double patch_size = 0.;
            for(int j=0; j<num_samples; j++) {
                for(int i=0; i<num_samples; i++) {
                    int index = patch_sampling_index(i,j);
                    if(index!=-1) {
                        Point2 xy;
                        if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
                            xy = Point2(parametric_samples.clmdata(j*num_samples+i));
                        } else {
                            // MJM 1/2018 non-standard spacing for blendsurf points.
                            // leaving this alone.
                            //xy[0] = init + i*step;
                            //xy[1] = init + j*step;
                            xy = init*Point2(1) + Point2(i,j)*step;
                        }

                        FacePointOverlapping* face_point_arr = (FacePointOverlapping*)(face_point.clmdata(index));			 
                        curpch->xy_to_face_point(xy.array(), face_point_arr);

                        double alpha;  
                        Point3 pdd[3];
                        // DZ: bad name for xy_to_patch_coords -- produces 3d position?
                        // MJM: Suggestions for alternative? produces 3d position on
                        // surface (pdd[0]) and 1st partial derivatives in x/y
                        // direction 
                        curpch->xy_to_patch_coords(xy.array(), EVAL_VL|EVAL_FD, (double*)pdd);
                        // blending function value at xy
                        curpch->xy_to_patch_value(xy.array(), EVAL_VL, &alpha);
                        blend_func_value(index) = alpha;

                        // 3d position 
                        for(int d=0; d<dim(); d++) 
                            position(d,index) = pdd[0](d);

                        // normal
                        Point3 cn(cross(pdd[1],pdd[2]));
                        double len = cn.l2();

                        for(int d=0; d<dim(); d++)
                            normal(d,index) = cn[d]/len;
                        // jacobian
                        jacobian(index) = len;

                        // trapezoidal quadrature weight
                        quad_weight(index) = step*step;

                        if(dynamic_cast<PatchSurfFaceMap*>(_bdry) != NULL){
                            int on_surface_point_index = pi*num_samples*num_samples + index;
                            OnSurfacePoint on_surface_point(
                                    pi,         //parent patch
                                    1e-16,      // distance to target
                                    Point2(xy), // (u,v) coordinates
                                    NEAR,       // region marking
                                    on_surface_point_index); // corr. target id
                            on_surface_point.inside_domain = INSIDE; 
                            _sample_point_as_on_surface_point(on_surface_point_index) = 
                                on_surface_point;
                            
                            int level = FaceMapSubPatch::as_subpatch(curpch)->_level;
                            quad_weight(index) = lagrange_basis_integrals(i,j)/double(pow(4,level));
                            alpha = 1.;
                        }
                        patch_size += jacobian(index)*quad_weight(index);
                        /* For face-map only:
                         * approximate the discrete density samples of \phi on a patch by
                         * 
                         * \phi(x_k, y_k) = \phi_k = 
                         *     \sum_i=0^n \sum_j=0^n \phi_ij*L_i(x_k)*L_j(x_k)
                         * 
                         * where L_i(x) is the ith lagrange basis polynomial.
                         *
/bin/bash: :341: command not found
                         * weight needs an additional multplicative factor of 
                         * \int_P L_i(x)L_j(y) dP.
                         */

                        //far_field(index) = _spacing;
                        // assumes spacing and distance are hard coded
                        far_field(index) = _boundary_distance_ratio;
                        interpolant_spacing(index) = _interpolation_spacing_ratio;
                        
                        //combined integration weight
                        combined_weight(index) = len * alpha * quad_weight(index) ; //jac * alf * wgt

                        // combined integration weight
                        //combined_weight(index) = len * alpha * step*step; //jac * alf * wgt

                        parametric_preimage(0,index) = xy.x();
                        parametric_preimage(1,index) = xy.y();

                        // dominant flags
                        Tag* tagptr = (Tag*)&(properties(index));
                        tagptr->_gid = curpch->group_id();

                        bool dominant;
                        curpch->is_xy_dominant(xy.array(), dominant);

                        tagptr->_dmt = dominant;
                        sum += combined_weight(index);
                    }
                }
            }
            
            // computing patch size
            curpch->characteristic_length() = sqrt(patch_size);
            for(int j=0; j<num_samples; j++) {
                for(int i=0; i<num_samples; i++) {
                    //Integral of basis function L_i(x)*L_i(y)
                    int index = patch_sampling_index(i,j);
                    if(index!=-1) {
                        double L;
                        if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
                            L = curpch->characteristic_length(); // spacing is param space chebyshev spacing
                        } else {
                            L = 1.; // spacing is blendsurf direct 3 spacing
                        }

                        // For patch refinement of a large face-map patch,
                        // one must scale the contribution of a subpatch to
                        // be be proportional to it's size. Since the domain
                        // of a subpatch is still [0,1]^2, we need to scale 
                        // the quadrature weight by 1/(4^level).
                        if(!qbkix_classical_conv && qbkix_adaptive_conv){
                            // scale point spacing with L
                            far_field(index) = L*_boundary_distance_ratio;
                            interpolant_spacing(index) = L*_interpolation_spacing_ratio;
                        } else if(qbkix_classical_conv && !qbkix_adaptive_conv){
                            // scale point spacing with sqrt(L)
                            far_field(index) = sqrt(L)*_boundary_distance_ratio;
                            interpolant_spacing(index) = sqrt(L)*_interpolation_spacing_ratio;
                        } else{
                            assert(0);
                        }
                    }

                }
            }

        }

    }// end of patch loop
    cout <<"SUM "<<sum<<endl;
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::setup"
int PatchSamples::setup(bool refined)
{
    if(!refined){
        stats.start_timer("PatchSamples::setup()");
    }
    ebiFunctionBegin;
    
    // TODO MJM remove this pattern
    // translate "solver spacing" to "surface discretization spacing"
    // DZ what's the difference? Why is this needed?
    /*
    // MJM BUG this code is called in the executable, and passed to PatchSurf
    // surface rep; since _spacing was never used in PatchSurf, it's a wonder why
    // it was ever done.
    */

    // If we're using a refined surface (i.e. for near-evaluation) and blendsurf, reduce the
    // spacing appropriately.
    if (refined && dynamic_cast<PatchSurfBlended*>(_bdry)){
        if(_spacing <= 0.2 && _spacing >= 0.025){
            _spacing /= 2.; 
        } else if (_spacing <= 0.025){
            _spacing /= 3.; 
        } else {
            //assert(0);
        }
            _spacing = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
    }
    if (refined && !dynamic_cast<PatchSurfBlended*>(_bdry)){
            _spacing = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
    }

    cout << "SPACING = " << _spacing << endl;


    vector<Patch*>& patches = _bdry->patches(); //get patches

    int num_patches = patches.size();

    stats.add_result("Number of patches", num_patches);
    initialize_sampling_vectors(num_patches);

    //1. compute  step_size, num_sample_points
    if(dynamic_cast<PatchSurfBlended*>(_bdry)){
        set_equal_sample_rate_physical_space(patches);
    } else {
        set_equal_sample_rate_param_space(patches);
    }
    
    //1.1 decide the step and num_samples	 
    

    DblNumMat parametric_samples;
    if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
        int num_samples = _num_sample_points[0];
        // tensor-product chebyshev for face-map
        parametric_samples = DblNumMat(2, num_samples*num_samples, true, 
                sample_2d<chebyshev1>(num_samples, Sampling::base_domain).data());
    }

    for(int pi=0; pi<num_patches; pi++)	 
        cout<<_num_sample_points[pi]<<" ";  cout<<endl;

    //2. patch_sampling_index, num_sample_points_in_patch, sample_point_starting_index
    
    int local_num_sample_points = initialize_sampling_indices();

    //3. create all per sample point distrib storage
    
    initialize_parallel_vectors(local_num_sample_points);
    double sum = 0;
    //int64_t interpolant_num_samples = Options::get_int_from_petsc_opts("-LL");


    // Cache quadrature weights; we alwasys used num_samples x num_samples
    // points per patch, so we'll store the weights so we don't need to
    // recompute the integrals each time.
    int num_samples = 0; 
    if(dynamic_cast<PatchSurfFaceMap*>(_bdry) != NULL){
            // TODO BUG this break integral caching below if face-map patches 
            // have non-equal points per patch
            for(int pi=0; pi<num_patches; pi++) {
                for(int pj=0; pj<num_patches; pj++) {
                    ebiAssert(_num_sample_points[pi] == _num_sample_points[pj]);
                }
            }
            num_samples = _num_sample_points[0];
    } 
    _sample_point_as_on_surface_point.resize(num_patches*num_samples*num_samples);

    //DblNumMat lagrange_basis_integrals = integrate_tensor_product_lagrange_basis_funcs(num_samples);

    

    string qbkix_convergence_type = Options::get_string_from_petsc_opts("-qbkix_convergence_type");
    bool qbkix_classical_conv = qbkix_convergence_type == "classic";
    bool qbkix_adaptive_conv = qbkix_convergence_type == "adaptive";
    assert(qbkix_classical_conv || qbkix_adaptive_conv);

    // initialize per sample point quantities: face_point representations, 
    // 3d positions, normals, jacobians, blending func values, quad weights, 
    // combined weights, dominant tags
    
    // Precompute density barycentric interpolation weights 
    // Pre-compute interpolation weights
    if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
        _interpolation_nodes_x=DblNumVec(num_samples);
        _interpolation_nodes_y=DblNumVec(num_samples);
        for(int i =0; i < num_samples; i++){
            _interpolation_nodes_x(i) = parametric_samples(0,i);
            _interpolation_nodes_y(i) = parametric_samples(1,i*num_samples); 
        }
        cout << "PatchSamples init barycentric_weights 1d" << endl;

        _barycentric_weights_x = 
            Interpolate::compute_barycentric_weights_1d<double>(_interpolation_nodes_x);
        _barycentric_weights_y = 
            Interpolate::compute_barycentric_weights_1d<double>(_interpolation_nodes_y);
        cout << "PatchSamples inited barycentric_weights 1d" << endl;
    }
    
    
    double max_patch_size = -DBL_MAX;
    sample_patches();

    if(!refined){
        stats.stop_timer("PatchSamples::setup()");
    }
    if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
        for(auto p: _bdry->patches()){
            max_patch_size = max(max_patch_size, 
                    FaceMapSubPatch::as_subpatch(p)->characteristic_length());
        }
#pragma omp parallel for
        for(int i =0; i < _bdry->num_patches(); i++){ 
            auto p = _bdry->patch(i);
            FaceMapSubPatch::as_subpatch(p)->compute_near_zone_distance();
        }
    }
    stats.add_result("face-map max patch size", max_patch_size);

    ebiFunctionReturn(0);
}


// ---------------------------------------------------------------------- 
// refine data on the patches of the  surface in this partition, returning the result in refined_datavec 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::refine_data"
int PatchSamples::refine_data(int dof, int refinement_factor, Vec dat, 
        vector<DblNumMat>& refined_datvec)
{
  ebiFunctionBegin;
  
  vector<Patch*>& patches = _bdry->patches(); //get patches
  int num_patches = patches.size();
  refined_datvec.resize(num_patches);

  // for each patch in the current partition
  for(int pi=0; pi<num_patches; pi++) {
	 if(_patch_partition[pi]==mpiRank()) {

		int num_samples = _num_sample_points[pi];
		IntNumMat& patch_sampling_index = _patch_sampling_index[pi];
		DblNumMat sample_point_data(this->sample_point_data(pi, dof, dat));
		DblNumMat sample_point_parametric_preimage(this->sample_point_parametric_preimage(pi));

		// create a square array of size N x N where N is the patch resolution
		// fill it with data values at all valid points, zero at points outside the patch
		DblNumMat patch_samples(dof,num_samples*num_samples); clear(patch_samples);
		for(int j=0; j<num_samples; j++) {
		  for(int i=0; i<num_samples; i++) {

			 int index = patch_sampling_index(i,j);
			 if(index!=-1) {
				for(int d=0; d<dof; d++) 
                    patch_samples(d, i+j*num_samples) = sample_point_data(d,index);
			 }
		  }
		}
		int refined_num_samples = num_samples*refinement_factor;
        DblNumMat& refined_dat = refined_datvec[pi];
        refined_dat.resize(dof, refined_num_samples*refined_num_samples); clear(refined_dat);

        //refinement using fft: do fourier transform, increase the
        //number of samples by refinement factor padding with zeros, transform back 
        int mn[2];
        mn[0] = num_samples;
        mn[1] = num_samples;

        int rs[2];
        rs[0] = refinement_factor;
        rs[1] = refinement_factor;
        if(dynamic_cast<PatchSurfFaceMap*>(_bdry)){
                 ebiAssert(false);
        } else {
            fftrf2d(dof, patch_samples.data(), mn, rs, refined_dat.data());
        }
        //cout << "pi: " << pi << endl << refined_dat.m() << ", " << refined_dat.n() << " init: " << num_samples << endl;
     }
  }
  ebiFunctionReturn(0);
}

// Interpolate density from num_samples*num_samples tensor-product Chebyshev
// points to (num_samples*refinement_factor)*(num_samples*refinement_factor)
// tensor-product Chebyshev points
Vec PatchSamples::refine_density(int dof, int ref_factor, Vec density, 
        PatchSamples* refined_patch_samples)
{
  ebiFunctionBegin;
  
  vector<Patch*>& patches = _bdry->patches(); //get patches
  int num_patches = patches.size();
    double spacing, refined_spacing;
  PetscBool err;
  PetscOptionsGetReal(NULL, "",  "-bis3d_rfdspacing", &refined_spacing,  &err);
  assert(err);
  PetscOptionsGetReal(NULL, "",  "-bis3d_spacing", &spacing,  &err);
  assert(err);
  int refinement_factor = int(spacing/refined_spacing);




  int64_t num_local_points;
  VecGetSize(density, &num_local_points);
  int64_t num_refined_points = refined_patch_samples->local_num_sample_points();
  cout << num_local_points <<", "<< num_refined_points << ",  " << 
          num_local_points*refinement_factor*refinement_factor << endl;

  Vec refined_density;
  VecCreateMPI(this->mpiComm(),
          num_refined_points*dof,
          PETSC_DETERMINE,
          &refined_density);
    // for each patch in the current partition
  for(int pi=0; pi<num_patches; pi++) {
	 if(_patch_partition[pi]==mpiRank()) {

		int num_samples = _num_sample_points[pi];
		//IntNumMat& patch_sampling_index = _patch_sampling_index[pi];

        // Interpolation points and function values to interpolate
		DblNumMat sample_point_density(this->sample_point_data(pi, dof, density));
		DblNumMat sample_point_parametric_preimage(this->sample_point_parametric_preimage(pi));

        // samples to interpolate to and where to store result
		DblNumMat refined_sample_point_density(refined_patch_samples->sample_point_data(pi, dof, refined_density));
		DblNumMat refined_sample_point_parametric_preimage(refined_patch_samples->sample_point_parametric_preimage(pi));
        DblNumMat refined_density_samples(dof, pow(refinement_factor*num_samples-1,2));
        
        
        Interpolate::evaluate_barycentric_interpolant_2d(dof, 
                sample_point_parametric_preimage, 
                sample_point_density, 
                num_samples,
                refined_sample_point_parametric_preimage, 
                refined_sample_point_density);
     }
  }
  
    return refined_density;
}

Vec interpolate_and_resample(int dof, Vec function,
        PatchSamples* interp_samples,
        PatchSamples* eval_samples){

    cout << "Interpolate and resample density" << endl;
    auto surface = dynamic_cast<PatchSurfFaceMap*>(interp_samples->bdry());
    assert(surface);
    assert( dynamic_cast<PatchSurfFaceMap*>(eval_samples->bdry()));
    assert( interp_samples->bdry() == eval_samples->bdry());

    // for each patch
    //      get the coarse func_values and uv_coordinates...
    //      and the  fine func_values and uv_coordinates to evaluate at 

    // figure out which fine patches are contained in each coarse patch
    // TODO fix 
    // make num_interp_samples and num_eval_samples getting numpts from patchsamples
    /*double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    int num_samples_per_patch_1d = floor(1./spacing)+1; // TODO factor out bug prone
    int num_samples_per_patch = num_samples_per_patch_1d*num_samples_per_patch_1d;
*/
    int num_interp_samples_1d = interp_samples->num_sample_points()[0];
    int num_interp_samples_per_patch =num_interp_samples_1d*num_interp_samples_1d;
    int num_eval_samples_1d = eval_samples->num_sample_points()[0];
    int num_eval_samples_per_patch =num_eval_samples_1d*num_eval_samples_1d;
    int num_patches= surface->num_patches();
    // make the output vector of size
    // num_samples_per_patch*num_refined_patches*dof
    Vec refined_func;
    Petsc::create_mpi_vec(
            interp_samples->mpiComm(),
            num_eval_samples_per_patch*num_patches*dof,
            refined_func);
    cout << "num_eval_samples_per_patch: " << num_eval_samples_per_patch << endl;
    cout << "num_patches: " << num_patches << endl;

    // for each refined patch, interpolate values from parent patch values 
    //for(auto patch : surface->patches()){
    for(int pi= 0; pi < num_patches; pi++){

        //int parent_face_map_id = subpatch->_face_map_patch->V();
        //auto parent_patch = dynamic_cast<FaceMapSubPatch*>(coarse_surf->patches()[parent_id]);
        //auto parent_patch = surface->subpatch(parent_id);
        
        // Interpolation points and function values to interpolate
        DblNumMat interp_function_values(interp_samples->sample_point_data(pi, dof, function));
        DblNumMat interp_sample_uv_coords(interp_samples->sample_point_parametric_preimage(pi)); //2 x num_refined_patches*points_per_patch


        // samples to interpolate to and where to store result
        DblNumMat eval_function_values(eval_samples->sample_point_data(pi, dof, refined_func));
        DblNumMat eval_sample_uv_coords(eval_samples->sample_point_parametric_preimage(pi));
        
        /*// unscale samples to interpolate properly
        for(int si =0; si < fine_sample_uv_coords.n(); si++){
            // shift from [0,1]^2 w.r.t to the child patch to [a,b]^2 \subset
            // [0,1]^2 w.r.t the parent patch
            Point2 subpatch_sample(fine_sample_uv_coords.clmdata(si));
            Point2 patch_sample;
            subpatch->rescale(subpatch_sample.array(), patch_sample.array());
            


            for(int d =0; d < 2; d++){
                fine_sample_uv_coords_scaled(d,si) = patch_sample(d);
            }
        }*/

        // interpolate values 
        Interpolate::evaluate_barycentric_interpolant_2d(dof, 
                interp_sample_uv_coords, 
                interp_function_values, 
                num_interp_samples_1d,
                eval_sample_uv_coords, 
                eval_function_values);
        
    }

    return refined_func; 
}


Vec refine_function(int dof, Vec function,
        PatchSamples* coarse_samples,
        PatchSamples* fine_samples){

    auto coarse_surf = dynamic_cast<PatchSurfFaceMap*>(coarse_samples->bdry());
    auto fine_surf   = dynamic_cast<PatchSurfFaceMap*>(fine_samples->bdry());
    assert(coarse_surf && fine_surf );

    // figure out which fine patches are contained in each coarse patch
    // TODO use map below for vectorizing interpolation
    //PatchChildrenMap subpatch_map = coarse_surf->find_subpatches(fine_surf);
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    int num_samples_per_patch_1d = floor(1./spacing)+1; // TODO factor out bug prone
    int num_samples_per_patch = num_samples_per_patch_1d*num_samples_per_patch_1d;

    int num_refined_patches = fine_surf->patches().size();
    // make the output vector of size
    // num_samples_per_patch*num_refined_patches*dof
    Vec refined_func;
    Petsc::create_mpi_vec(
            coarse_samples->mpiComm(),
            num_samples_per_patch*num_refined_patches*dof,
            refined_func);


    // for each refined patch, interpolate values from parent patch values 
    // TODO vectorize somehow to do a single interpolation per coarse patch
#pragma omp parallel for
    for(int pi = 0; pi < fine_surf->num_patches(); pi++){
        auto subpatch = fine_surf->subpatch(pi);
        //auto subpatch = FaceMapSubPatch::as_subpatch(patch);
        int child_id = subpatch->V();
        int parent_id = subpatch->_coarse_parent_patch;

        //int parent_face_map_id = subpatch->_face_map_patch->V();
        //auto parent_patch = dynamic_cast<FaceMapSubPatch*>(coarse_surf->patches()[parent_id]);
        auto parent_patch = coarse_surf->subpatch(parent_id);
        
        // Interpolation points and function values to interpolate
        DblNumMat coarse_function_values(coarse_samples->sample_point_data(parent_id, dof, function));
        DblNumMat coarse_sample_uv_coords(coarse_samples->sample_point_parametric_preimage(parent_id)); //2 x num_refined_patches*points_per_patch
        DblNumMat coarse_sample_uv_coords_scaled(coarse_sample_uv_coords.m(), coarse_sample_uv_coords.n());

        for(int si =0; si < coarse_sample_uv_coords.n(); si++){
            // shift from [0,1]^2 w.r.t to the child patch to [a,b]^2 \subset
            // [0,1]^2 w.r.t the parent patch
            Point2 subpatch_sample(coarse_sample_uv_coords.clmdata(si));
            Point2 patch_sample;
            parent_patch->rescale(subpatch_sample.array(), patch_sample.array());
            

            for(int d =0; d < 2; d++){
                coarse_sample_uv_coords_scaled(d,si) = patch_sample(d);
            }
        }

        // samples to interpolate to and where to store result
        DblNumMat fine_function_values(fine_samples->sample_point_data(child_id, dof, refined_func));
        DblNumMat fine_sample_uv_coords(fine_samples->sample_point_parametric_preimage(child_id));
        DblNumMat fine_sample_uv_coords_scaled(fine_sample_uv_coords.m(), fine_sample_uv_coords.n());
        
        // unscale samples to interpolate properly
        for(int si =0; si < fine_sample_uv_coords.n(); si++){
            // shift from [0,1]^2 w.r.t to the child patch to [a,b]^2 \subset
            // [0,1]^2 w.r.t the parent patch
            Point2 subpatch_sample(fine_sample_uv_coords.clmdata(si));
            Point2 patch_sample;
            subpatch->rescale(subpatch_sample.array(), patch_sample.array());
            


            for(int d =0; d < 2; d++){
                fine_sample_uv_coords_scaled(d,si) = patch_sample(d);
            }
        }
        // interpolate values 
        Interpolate::evaluate_barycentric_interpolant_2d(dof, 
                coarse_sample_uv_coords_scaled, 
                coarse_function_values, 
                num_samples_per_patch_1d,
                /*coarse_samples->interpolation_nodes_x(),
                coarse_samples->interpolation_nodes_y(),
                coarse_samples->barycentric_weights_x(),
                coarse_samples->barycentric_weights_y(),*/

                fine_sample_uv_coords_scaled, 
                fine_function_values);
        /*
        void Interpolate::evaluate_barycentric_interpolant_2d(
                int dof,
                DblNumMat xy_coordinates, 
                DblNumMat function_values,
                int num_samples_1d,
                DblNumVec interpolation_nodes_x,
                DblNumVec interpolation_nodes_y,
                DblNumVec barycentric_weights_x,
                DblNumVec barycentric_weights_y,
                DblNumMat refined_xy_coordinates,
                DblNumMat& refined_function_values){
         */
    }

    return refined_func; 
}




// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::is_interp_xy_valid"
int PatchSamples::is_interp_xy_valid(int pi, double* xy, bool& is_valid)
{
  ebiFunctionBegin;
  vector<Patch*>& patches = _bdry->patches(); //get patches
  patches[pi]->is_xy_valid(xy, is_valid);
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::interpolated_position_and_derivatives"
int PatchSamples::interpolated_position_and_derivatives(int pi, double* xy, 
        int flag, double* res)
{
  ebiFunctionBegin;
  vector<Patch*>& patches = _bdry->patches(); //get patches
  //double start = omp_get_wtime();
  patches[pi]->xy_to_patch_coords(xy, flag, res);
  //stats.result_plus_equals("singular surface eval time", omp_get_wtime() - start);
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
// On patch pi, given values on a refined grid, interpolate at a  point
// xy,return result in res
// flag indicates whether to return value, value + first derivative, or
// value + second derivative, all packed into res sequentially 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::interpolate_data"
int PatchSamples::interpolate_data(int pi, double* xy, int refinement_factor, 
        int surface_interpolation_num_samples,
        int flag, vector<DblNumMat>& refined_datvec, double* res)
{
  ebiFunctionBegin;
  vector<Patch*>& patches = _bdry->patches(); //get patches
  Patch* curpch = patches[pi];

  double bnd = curpch->bnd();
  double refined_init = -bnd;
  int    refined_num_samples = _num_sample_points[pi]*refinement_factor;
  double refined_step = _step_size[pi]/double(refinement_factor);

  DblNumMat& refined_dat = refined_datvec[pi];
  int dof = refined_dat.m(); //degree of freedom
 surface_interpolation_num_samples = 4; 
  int mn[2];
  mn[0]=refined_num_samples;
  mn[1]=refined_num_samples;

  double ef[2];
  ef[0] = 2.0*bnd;
  ef[1] = 2.0*bnd;

  double ab[2];
  ab[0] = (xy[0]-refined_init)/refined_step;
  ab[1] = (xy[1]-refined_init)/refined_step;

  int    ij[2];
  ij[0] = (int)floor(ab[0]);
  ij[1] = (int)floor(ab[1]);

  if (ij[0] >= refined_num_samples || ij[1] >= refined_num_samples){
	 std::cout << ij[0] << " " << ij[1] << " refined_ " << refined_num_samples << endl;
  }

  if (ij[0] >= refined_num_samples) 
      ij[0] = refined_num_samples-1;

  if (ij[1] >= refined_num_samples) 
      ij[1] = refined_num_samples-1;

  ebiAssert(ij[0]>=0 && ij[0]<refined_num_samples);
 	
  ebiAssert(ij[1]>=0 && ij[1]<refined_num_samples);

  double uv[2];
  uv[0] = ab[0] - ij[0];
  uv[1] = ab[1] - ij[1];
  // important: _surface_interpolation_num_samples argument is semi-useless;
  // lagev2d function
  // is half-coverted from hard-coded LL=4 to arbitrary LL, but it does not
  // seem it would work properly for any  other LL
  
  //if(dynamic_cast<PatchSurfFaceMap*>(_bdry) != NULL){
      //lagev2d(flag, DMFLAG_CLOSED, dof, refined_dat.data(), mn, ef, ij, uv, res, surface_interpolation_num_samples);
  //} else {
  lagev2d(flag, DMFLAG_PERIOD, dof, refined_dat.data(), mn, ef, ij, uv, res, surface_interpolation_num_samples);
  //}
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::is_sample_point_valid"
int PatchSamples::is_sample_point_valid(int pi, int* ij, bool& is_valid)
{
  ebiFunctionBegin;

  IntNumMat& patch_sampling_index = _patch_sampling_index[pi];
  int index = patch_sampling_index(ij[0],ij[1]);
  is_valid = (index!=-1);

  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::get_sample_on_patch"
int PatchSamples::get_sample_on_patch(int pi, int* ij, double* pos, double* nor, double* jac)
{
  ebiFunctionBegin;
  
  IntNumMat& patch_sampling_index = _patch_sampling_index[pi];
  int index = patch_sampling_index(ij[0],ij[1]);
  
  DblNumMat tpos( sample_point_3d_position(pi) );
  DblNumMat tnor( sample_point_normal(pi) );
  DblNumVec tjac( sample_point_jacobian(pi) );

  for(int d=0; d<dim(); d++) 
      pos[d] = tpos(d,index);
  
  for(int d=0; d<dim(); d++) 
      nor[d] = tnor(d,index);

  jac[0] = tjac(index);
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSamples::get_sample_point"
int PatchSamples::get_sample_point(int pi, int* ij, int dof, Vec dat, double* res)
{
  ebiFunctionBegin;
  //int num_samples = _num_sample_points[pi];
  IntNumMat& patch_sampling_index = _patch_sampling_index[pi];
  int index = patch_sampling_index(ij[0],ij[1]);
  DblNumMat sample_point_data( this->sample_point_data(pi, dof, dat) );
  for(int d=0; d<dof; d++)
	 res[d] = sample_point_data(d,index);
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
DblNumMat PatchSamples::sample_as_face_point(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_as_face_point, &arr);
  double* buf=arr; VecRestoreArray(_sample_as_face_point, &arr);
  return DblNumMat(face_point_size_in_doubles(), _num_sample_points_in_patch[pi], false, buf + dif*face_point_size_in_doubles());
}

int PatchSamples::sample_point_local_start_index(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  return dif;
}
  
DblNumMat PatchSamples::sample_point_3d_position(int pi)
{ 
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_3d_position, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_3d_position, &arr);
  return DblNumMat(dim(),     _num_sample_points_in_patch[pi], false, buf + dif*dim());
}
DblNumMat PatchSamples::sample_point_normal(int pi)
{ 
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_normal, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_normal, &arr);
  return DblNumMat(dim(),     _num_sample_points_in_patch[pi], false, buf + dif*dim());
}
DblNumVec PatchSamples::sample_point_jacobian(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_jacobian, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_jacobian, &arr);
  return DblNumVec(_num_sample_points_in_patch[pi], false, buf + dif);
}
DblNumVec PatchSamples::sample_point_blend_func_value(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_blend_func_value, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_blend_func_value, &arr);
  return DblNumVec(_num_sample_points_in_patch[pi], false, buf + dif);
}
DblNumVec PatchSamples::sample_point_quad_weight(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_quad_weight, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_quad_weight, &arr);
  return DblNumVec(_num_sample_points_in_patch[pi], false, buf + dif);
}
DblNumVec PatchSamples::sample_point_props(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_props, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_props, &arr);
  return DblNumVec(_num_sample_points_in_patch[pi], false, buf + dif);
}

DblNumVec PatchSamples::sample_point_far_field(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_far_field, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_far_field, &arr);
  return DblNumVec(_num_sample_points_in_patch[pi], false, buf + dif);
}
DblNumVec PatchSamples::sample_point_interpolant_spacing(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_interpolant_spacing, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_interpolant_spacing, &arr);
  return DblNumVec(_num_sample_points_in_patch[pi], false, buf + dif);
}


DblNumMat PatchSamples::sample_point_parametric_preimage(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_parametric_preimage, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_parametric_preimage, &arr);
  return DblNumMat(2, _num_sample_points_in_patch[pi], false, buf + 2*dif);
}


DblNumMat PatchSamples::sample_point_data(int pi, int adof, Vec adat)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    adat, &arr);
  double* buf=arr; VecRestoreArray(adat, &arr);
  return DblNumMat(adof, _num_sample_points_in_patch[pi], false, buf + dif*adof);
}
DblNumVec PatchSamples::sample_point_combined_weight(int pi)
{
  int a,b; local_sample_point_range(a,b);
  int s = _sample_point_starting_index[pi];  assert(s>=a && s<b);
  int dif = s - a;
  double* arr;     VecGetArray(    _sample_point_combined_weight, &arr);
  double* buf=arr; VecRestoreArray(_sample_point_combined_weight, &arr);
  return DblNumVec(_num_sample_points_in_patch[pi], false, buf + dif);
}

Vec PatchSamples::generate_qbkix_points_from_sample_points(vector<int> qbkix_indices, 
        vector<int> patches_to_sample){
    
    Vec interpolation_directions;
    VecDuplicate(this->sample_point_normal(), &interpolation_directions);
    VecCopy(this->sample_point_normal(), interpolation_directions);
    VecScale(interpolation_directions, -1.);
    // generate qbkix targets from samples
    // qbkix point gen change
    // Create a Petsc IndexSet to extract the sample point data associated with
    // the patch id's listed in patches_to_sample. @BUG this most definitely
    // breaks in parallel :), and for non-equal number of samples per patch.
    for(int pi = 0; pi < _bdry->patches().size(); pi++){
        assert(_num_sample_points[pi] == _num_sample_points[0]);
    }

    int num_samples_per_patch = _num_sample_points[0]*_num_sample_points[0];

    cout << "qbkix points to generate:" << endl;
    for(int qi  =0; qi < qbkix_indices.size(); qi++)
        cout << qbkix_indices[qi] << ", ";
    cout << endl;

    
    // create an index set into the total set of patch sample info for the
    // subset of patches we need to generate qbkix points for...
    vector<int64_t> patches_to_sample_64t(patches_to_sample.begin(), patches_to_sample.end());
    IS patch_to_sample_point_index_set;
    IS patch_to_sample_value_index_set;
    ISCreateBlock(this->mpiComm(),
        num_samples_per_patch*DIM,
        patches_to_sample_64t.size(),
        patches_to_sample_64t.data(),
        PETSC_COPY_VALUES,
        &patch_to_sample_point_index_set);
    ISCreateBlock(this->mpiComm(),
        num_samples_per_patch,
        patches_to_sample_64t.size(),
        patches_to_sample_64t.data(),
        PETSC_COPY_VALUES,
        &patch_to_sample_value_index_set);
    
    Vec sample_point_3d_position_subset;
    Vec interpolation_directions_subset;
    Vec sample_point_far_field_subset;
    Vec sample_point_interpolant_spacing_subset;
    
    VecGetSubVector(this->sample_point_3d_position(),
            patch_to_sample_point_index_set,
            &sample_point_3d_position_subset);
    
    VecGetSubVector(interpolation_directions,
            patch_to_sample_point_index_set,
            &interpolation_directions_subset);

    VecGetSubVector(this->sample_point_far_field(),
            patch_to_sample_value_index_set,
            &sample_point_far_field_subset);
    
    VecGetSubVector(this->sample_point_interpolant_spacing(),
            patch_to_sample_value_index_set,
            &sample_point_interpolant_spacing_subset);

    Vec qbkix_points = 
        generate_interior_qbkix_points(
            sample_point_3d_position_subset,
            sample_point_3d_position_subset,
            qbkix_indices,
            interpolation_directions_subset,
            sample_point_far_field_subset,
            sample_point_interpolant_spacing_subset);
    
    VecRestoreSubVector(this->sample_point_3d_position(),
            patch_to_sample_point_index_set,
            &sample_point_3d_position_subset);

    VecRestoreSubVector(this->sample_point_far_field(),
            patch_to_sample_value_index_set,
            &sample_point_far_field_subset);
    
    VecRestoreSubVector(this->sample_point_interpolant_spacing(),
            patch_to_sample_value_index_set,
            &sample_point_interpolant_spacing_subset);
    
    VecRestoreSubVector(interpolation_directions,
            patch_to_sample_point_index_set,
            &interpolation_directions_subset);
    VecDestroy(&interpolation_directions);
    
    ISDestroy(&patch_to_sample_point_index_set);
    ISDestroy(&patch_to_sample_value_index_set);
    
    return qbkix_points;
    
}

NumVec <OnSurfacePoint> PatchSamples::closest_points_to_qbkix_points(){
    int qbkix_order = Options::get_double_from_petsc_opts("-near_interpolation_num_samples");
    int num_samples = global_num_sample_points();
    NumVec<OnSurfacePoint> closest_points(qbkix_order*num_samples);
    for(int i= 0; i < num_samples; i++){
        for(int qi = 0; qi < qbkix_order; qi++){
            closest_points(i*qbkix_order+qi) = _sample_point_as_on_surface_point(i);
        }
    }
    return closest_points;
}

void generate_samples_on_child_patches(MPI_Comm comm, 
        PatchSurfFaceMap* face_map,
        int num_samples_1d,
        Vec& uv_coordinates_single_patch, 
        Vec& sample_points_on_chlildren,
        vector<int64_t> patches_to_sample
        ){

    int num_patches_to_sample = patches_to_sample.size();//face_map->num_patches();
    const int num_children = 4;
    int num_samples_per_patch = num_samples_1d*num_samples_1d;
    
    // uv_coordinates to evaluate patches at child nodes.
    Petsc::create_mpi_vec(comm, 
            2*num_samples_per_patch*num_children, 
            uv_coordinates_single_patch);
    // sample position of child nodes in 3d
    Petsc::create_mpi_vec(comm, 
            DIM*num_samples_per_patch*num_children*num_patches_to_sample, 
            sample_points_on_chlildren);
    
    // note assumes chebyshev 
    DblNumMat uv_coords_local(2, uv_coordinates_single_patch);
    for(int i=0; i < 2; i++){
        for (int j = 0; j < 2; j++) {
            // iterate over child domains [0,.5]^2, [.5,1]x [0,.5], [0,.5] x
            // [.5,1], [.5,1]^2
            Rectangle child_domain;
            child_domain.first.first  = i == 0 ?  0.: .5;
            child_domain.first.second = i == 0 ? .5 : 1.;
            child_domain.second.first = j == 0 ?  0.: .5;
            child_domain.second.second= j == 0 ? .5 : 1.;

            DblNumMat parametric_samples = 
                DblNumMat(2, num_samples_per_patch, true, 
                    sample_2d<chebyshev1>(num_samples_1d, child_domain).data());
            //cout << i << ", " << j << endl;
            //cout << parametric_samples << endl;
            // copy
            for (int k = 0; k < parametric_samples.n(); k++) {
                for (int d = 0; d < parametric_samples.m(); d++) {
            //cout << k << ", " <<d  << "; " << (2*i+j)*num_samples_1d*num_samples_1d +k << endl;
                   uv_coords_local(d,(2*i+j)*num_samples_per_patch + k) = parametric_samples(d,k); 
                }
            }
        }
    }
    
    //evaluate samples patches
    DblNumMat sample_points_on_chlildren_local(DIM, sample_points_on_chlildren);
#pragma omp parallel for
    for (int ppi = 0; ppi < num_patches_to_sample; ppi++) {
        int pi = patches_to_sample[ppi];
        auto patch = face_map->patch(pi);
        for (int i = 0; i < uv_coords_local.n(); i++) {
            Point3 sample;
            Point2 uv(uv_coords_local.clmdata(i));
            patch->xy_to_patch_coords(uv.array(), PatchSamples::EVAL_VL, sample.array());
            for (int d = 0; d < DIM; d++) {
                sample_points_on_chlildren_local(d, ppi*num_samples_per_patch*num_children + i) = sample(d);
            }
        }
        
    }
}

END_EBI_NAMESPACE
