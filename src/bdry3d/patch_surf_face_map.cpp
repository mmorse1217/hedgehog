#include "patch_surf_face_map.hpp"
#include "face_map_subpatch.hpp"
#include <blended_evaluator.hpp>
#include <polynomial_evaluator.hpp>
#include <analytic_evaluator.hpp>
#include "common/interpolate.hpp"
#include "common/utils.hpp"
#include <sampling.hpp>

// TODO abstract library calls out
//#include "p4est_interface.hpp"
#include "p4est_refinement.hpp"
#include "bie3d/solver_utils.hpp"
extern "C" {
#include <p4est.h>
#include <p4est_vtk.h>
}

#include <evalfromw.hpp>
#include <bdsurf.hpp>

using std::istringstream;
using Sampling::sample_1d;
using Sampling::sample_2d;
using Sampling::base_domain;
using Sampling::chebyshev2;

BEGIN_EBI_NAMESPACE

PatchSurfFaceMap::PatchSurfFaceMap(const string& n, const string& p):
  PatchSurf(n,p), _surface_type(NOT_SET), _coarse(false) {;} 
PatchSurfFaceMap::~PatchSurfFaceMap() {
    delete _face_map;
}

// ---------------------------------------------------------------------- 
int PatchSurfFaceMap::setFromOptions()
{
    ebiAssert(_surface_type !=  NOT_SET);
    SurfaceEvaluator* evaluator;
    switch(_surface_type){
        case BLENDED:
            {
                string fileb = Options::get_string_from_petsc_opts("-bd3d_meshfile");
                ifstream tmpb(fileb.c_str());
                string filea = Options::get_string_from_petsc_opts("-bd3d_bdsurf_submatlibfile");
                ifstream tmpa(filea.c_str());

                CCSubMatLib* submatlib = new CCSubMatLib();
                cerr<<"submatlib setup"<<endl;
                submatlib->setup(tmpa);
                string filebd = Options::get_string_from_petsc_opts("-bd3d_bdsurf_bdulibfile");
                ifstream tmpbd(filebd.c_str());

                BdULib* bdulib = new BdULib();
                cerr<<"bdulib setup"<<endl;
                bdulib->setup(tmpbd);

                int cl = Options::get_int_from_petsc_opts("-bd3d_bdsurf_ctrllvl");
                assert(cl>=0 && cl<=2);

                int pc = Options::get_int_from_petsc_opts("-bd3d_bdsurf_pouctrl");
                assert(pc==0 || pc==1 || pc==2);

                int ct = Options::get_int_from_petsc_opts("-bd3d_bdsurf_chttyp");
                int bt = Options::get_int_from_petsc_opts("-bd3d_bdsurf_bsstyp");
                int pp = Options::get_int_from_petsc_opts("-bd3d_bdsurf_stpppg");
                int lb = Options::get_int_from_petsc_opts("-bd3d_bdsurf_lb");
                int ub = Options::get_int_from_petsc_opts("-bd3d_bdsurf_ub");
                int gl = Options::get_int_from_petsc_opts("-bd3d_bdsurf_indepboun");
                double fs = Options::get_double_from_petsc_opts("-bd3d_bdsurf_flats");
                int cot = Options::get_int_from_petsc_opts("-bd3d_bdsurf_concat");
                int sd = Options::get_int_from_petsc_opts("-bd3d_bdsurf_spdeg");
                
                

                int poudeg = Options::get_int_from_petsc_opts("-bd3d_bdsurf_poubsdeg");
                int loadW = Options::get_int_from_petsc_opts("-bd3d_bdsurf_loadw");
                int tmp = Options::get_int_from_petsc_opts("-bd3d_bdsurf_g");
                string matrix_dir = Options::get_string_from_petsc_opts("-bd3d_bdsurf_matdir");
                int activevert = Options::get_int_from_petsc_opts("-bd3d_bdsurf_activevert");
                MatricesDir = matrix_dir.c_str();
                NUM_PATCHES = tmp;
               
                
                _pouctrl = Options::get_int_from_petsc_opts("-bd3d_bdsurf_pouctrl");
                double spacing= Options::get_double_from_petsc_opts("-bdsurf_interpolant_spacing");
                double interpolate= Options::get_double_from_petsc_opts("-bdsurf_interpolate");
                //_pdeg= Options::get_int_from_petsc_opts("-bd3d_bdsurf_pdeg");
                _pdeg = 5; // i think might be worse than the petsc error above.

                //FaceMapSurf& bdsurf = new FaceMapSurf(  );
                BdSurf* bdsurf = new BdSurf();
                bdsurf->setParams(cl,pc,ct,bt,pp,lb,ub, gl, fs, cot, sd ,poudeg, spacing, interpolate,0,loadW );

                //bdsurf->setup(tmpb,(void*)submatlib, filebd); 
                bdsurf->setup(tmpb,0,Point3(1.0),(void*)submatlib,(void*)bdulib); // blendsurf v3 (doesn't work yet)
                evaluator =  new BlendedEvaluator(bdsurf);
            } 
            break;
        case POLYNOMIAL:
            {
                
                string meshfile= Options::get_string_from_petsc_opts("-bd3d_meshfile");
                string polynomials = Options::get_string_from_petsc_opts("-poly_coeffs_file");

                evaluator = new PolynomialEvaluator(polynomials, meshfile);
            }
            break;
        case ANALYTIC:
            {
                int analytic_func = Options::get_int_from_petsc_opts("-analytic_function");
                string meshfile = Options::get_string_from_petsc_opts("-bd3d_meshfile");
                evaluator = new AnalyticEvaluator(meshfile, AnalyticEvaluator::select_function(analytic_func));
            }
            break;
        default: 
            ebiAssert(0);
    }
                _filename= Options::get_string_from_petsc_opts("-bd3d_filename").c_str();

    PetscBool flg = PETSC_FALSE;

    int64_t refinement_factor = 0;  
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-facemap_refinement_factor",&refinement_factor,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    int64_t patch_order= 0;  
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-facemap_patch_order",&patch_order,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    int64_t adaptive= 0;  
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-facemap_adaptive",&adaptive,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    double fit_accuracy = Options::get_double_from_petsc_opts("-bd3d_facemap_fit_accuracy");

    cout << "adaptive: " << adaptive << endl;
    cout << "patch_order: " << patch_order << endl;
    cout << "refinement_factor: " << refinement_factor << endl;
    cout << "fit_accuracy: " << fit_accuracy << endl;

    _face_map = new FaceMapSurf(evaluator);
    _face_map->_legacy = true;
    _face_map->set_params(refinement_factor, patch_order, adaptive, fit_accuracy );
    _p4est = _face_map->_p4est;
    _on_surface_threshold = Options::get_double_from_petsc_opts("-markgrid_on_surface_threshold");

    ebiFunctionReturn(0);
}

void PatchSurfFaceMap::setup_from_existing_face_map(PatchSurfFaceMap* face_map){

    //p4est_destroy(_p4est);
    //delete _face_map;

    this->_surface_type = face_map->_surface_type;
    this->_quadrature_weights = face_map->_quadrature_weights;
    _boundary_component_center = face_map->boundary_component_center();
    _boundary_orientation  = face_map->boundary_orientation();
    this->_face_map = face_map->_face_map; 
    this->initialize_with_existing_p4est(face_map->_p4est);
    //initialize_p4est_leaves_with_subpatches(_p4est, this);
    _on_surface_threshold = Options::get_double_from_petsc_opts("-markgrid_on_surface_threshold");
}

void PatchSurfFaceMap::initialize_with_existing_p4est(p4est_t*& p4est){
    _p4est = object_safe_p4est_copy(p4est);
    _patches = p4est_to_face_map_subpatches(_p4est, this);
    //initialize_p4est_leaves_with_subpatches(_p4est, this);
   

}
void PatchSurfFaceMap::partial_teardown(){
    for(const auto& patch : _patches){
        //cout << "deleting Patch [" << patch  << "]" << endl;
        delete patch;
    }
    _patches.resize(0);
}
// ---------------------------------------------------------------------- 
int PatchSurfFaceMap::setup()
{
    ebiFunctionBegin;
    _face_map->setup();
    //1. read in and generate all scbs (??)
    string filedata;  iC( file2string(mpiComm(), _filename.c_str(), filedata) ); //parallel communication
    cerr << _filename.c_str() << endl;
    istringstream fin(filedata);

    int pcnt = 0;
    for(int V=0; V<_face_map->_patches.size(); V++) {
        _patches.push_back(  new FaceMapPatch(this, pcnt, V) );
        pcnt++;
    }

    //3. centers
    GpMesh& gpmesh = _face_map->gpmesh();
    _boundary_component_center = gpmesh.get_interior_points();
    _boundary_orientation  = gpmesh.get_intpt_orientation();

    for(auto b : _boundary_orientation){
        cout << b << endl;
    }


    int num_samples = 20;
    double spacing = 1./double(num_samples-1);

    NumMatrix samples = sample_2d<chebyshev2>(num_samples, base_domain);

    DblNumVec quad_weight(num_samples*num_samples);


    DblNumVec quadrature_nodes(num_samples, true,
            sample_1d<chebyshev2>(num_samples,base_domain).data());
    
    DblNumVec* quadrature_weights = new DblNumVec(num_samples);
    
    for(int si = 0; si < num_samples; si++){
        (*quadrature_weights)(si) = 
            Interpolate::integrate_ith_lagrange_basis_func(
                si,0., 1., num_samples, quadrature_nodes, 50, 1.);

    }
    _quadrature_weights = quadrature_weights;

    for(int si = 0; si < num_samples; si++){
        for(int sj = 0; sj < num_samples; sj++){
            int index = si*num_samples + sj;

            quad_weight(index) = 
                (*quadrature_weights)(si)*(*quadrature_weights)(sj);

        }
    }

    for(int pi = 0; pi < _patches.size(); pi++){
        Patch* patch = _patches[pi];

        DblNumMat normals(DIM,num_samples*num_samples);
        for(int i = 0; i < num_samples*num_samples; i++){
            Point2 sample_uv(samples.clmdata(i));

            vector<Point3> values(3, Point3());
            patch->xy_to_patch_coords(sample_uv.array(), 
                    PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, 
                    (double*)values.data());

            Point3 normal = cross(values[1], values[2]);
            for(int d = 0; d < DIM; d++)
                normals(d, i) = normal(d);
        }

        patch->characteristic_length(normals, quad_weight);

    }
    initialize_p4est_leaves_with_subpatches(_p4est, this);
    ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int PatchSurfFaceMap::face_point_size_in_doubles() {
    return (sizeof(FacePointFaceMap)-1)/sizeof(double)+1;
}

// ---------------------------------------------------------------------- 
int PatchSurfFaceMap::patches_containing_face_point(FacePointOverlapping* face_point, vector<int>& pivec)
{
    ebiFunctionBegin;
    pivec.clear();
    FacePointFaceMap* bb = (FacePointFaceMap*) face_point;
    //FaceMapSurf& bdsurf = *_face_map;

    int F = bb->F();
    //double* cd = bb->cd();
    FaceMapPatch* patch = (FaceMapPatch*) _patches[F];
    bool valid  = false;
    patch->is_xy_valid(bb->cd() , valid);
    if( valid){
        pivec.push_back(F);

    }

    ebiFunctionReturn(0);
}

void PatchSurfFaceMap::refine(Vec near_points){
    PatchSurfFaceMap* f = this; 
    refine_patches_for_qbkix_point_location(_p4est, f);
    
    if(_surface_type != POLYNOMIAL){
        //p4est_balance(p4est, P4EST_CONNECT_FULL, NULL);
    }
        //p4est_vtk_write_file(_p4est, NULL, "face_map_refinement");

    if(_coarse)
        set_coarse_patch_ids(_p4est);
    vector<Patch*> subpatches = p4est_to_face_map_subpatches(_p4est, this);
    //_patches = (vector<Patch*>)subpatches;
    _patches = subpatches;
}

void PatchSurfFaceMap::refine_test(Vec near_points){
    
    //p4est_connectivity_t* connectivity = build_connectivity_from_face_map(this);
    //p4est_t* p4est = p4est_new(this->mpiComm(), connectivity, sizeof(RefinementData<FaceMapSubPatch>), NULL, NULL);
    
    if(_surface_type != POLYNOMIAL){
        //p4est_balance(p4est, P4EST_CONNECT_FULL, NULL);
    }
        //p4est_vtk_write_file(_p4est, NULL, "face_map_refinement");

        //initialize_p4est_leaves_with_subpatches(_p4est, this);
    if(_coarse)
        set_coarse_patch_ids(_p4est);
    vector<Patch*> subpatches = p4est_to_face_map_subpatches(_p4est, this);
    _patches = subpatches;
    //_patches = (vector<Patch*>)subpatches;
    /*
    _patches = subpatches;
    _p4est = p4est;
    _p4est_connectivity = connectivity;*/
}


void PatchSurfFaceMap::refine_uniform(int level){
    
    PatchSurfFaceMap* f = this; 
    refine_patches_uniform(level, _p4est, f);
    if(_surface_type != POLYNOMIAL){
        //p4est_balance(p4est, P4EST_CONNECT_FULL, NULL);
    }
        //p4est_vtk_write_file(_p4est, NULL, "face_map_refinement");

        //initialize_p4est_leaves_with_subpatches(_p4est, this);
    if(_coarse)
        set_coarse_patch_ids(_p4est);
    vector<Patch*> subpatches = p4est_to_face_map_subpatches(_p4est, this);
    _patches = subpatches;
    //_patches = (vector<Patch*>)subpatches;
    /*
    _patches = subpatches;
    _p4est = p4est;
    _p4est_connectivity = connectivity;*/
}

void PatchSurfFaceMap::resolve_rhs(FunctionHandle function, int range_dim, double abs_err){
    
    PatchSurfFaceMap* f = this; 
    resolve_function(_p4est, f, function, range_dim, abs_err);
    //balance(_p4est);
    set_coarse_patch_ids(_p4est);
    if(_surface_type != POLYNOMIAL){
        //p4est_balance(p4est, P4EST_CONNECT_FULL, NULL);
    }
        //p4est_vtk_write_file(_p4est, NULL, "face_map_refinement");

        //initialize_p4est_leaves_with_subpatches(_p4est, this);
    if(_coarse)
        set_coarse_patch_ids(_p4est);
    vector<Patch*> subpatches = p4est_to_face_map_subpatches(_p4est, this);
    _patches = subpatches;
    //_patches = (vector<Patch*>)subpatches;
    /*
    _patches = subpatches;
    _p4est = p4est;
    _p4est_connectivity = connectivity;*/
}

void PatchSurfFaceMap::save(string filename){
    _p4est->connectivity = _p4est_connectivity;
    store_p4est(filename, _p4est);
}

PatchSurfFaceMap* PatchSurfFaceMap::load(string filename){

    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->setFromOptions();
    face_map->setup();
    load_p4est(filename, face_map->mpiComm(), face_map->_p4est, &face_map->_p4est_connectivity);
    face_map->patches() = 
        p4est_to_face_map_subpatches(face_map->_p4est, face_map);
    return face_map;

}

// Method below implement requirements from Bd3dov.hpp. Functionality should be
// similar to that of agsurf.cpp, except for blended surfaces instead of
// analytic ones.

FaceMapPatch::FaceMapPatch(PatchSurf* b, int pi, int V): Patch(b, pi), _V(V)
{
  ebiFunctionBegin;
    estimate_jacobian(&_approx_jacobian);
  //FaceMapSurf& face_map = ((PatchSurfFaceMap*)bdry())->face_map();
  

 /* 
  _bnd = .5;
  
  double patch_size;
  Point3 origin[1]; 
  Point3 corner[1];
  
  double xy0[2]; 
  xy0[0] = 0.;
  xy0[1] = 0.;

  double xy1[2];
  xy1[0] = 0.;
  xy1[1] = 1.;

  xy_to_patch_coords(xy0, PatchSamples::EVAL_VL, (double*)origin);
  xy_to_patch_coords(xy1, PatchSamples::EVAL_VL, (double*)corner);
  patch_size = (origin[0]-corner[0]).l2();

  xy1[0] = 1.;
  xy1[1] = 0.;
  xy_to_patch_coords(xy1, PatchSamples::EVAL_VL, (double*)corner);
  patch_size *= (origin[0]-corner[0]).l2();
  _characteristic_length = sqrt(patch_size);
  */
  ebiFunctionReturnVoid;
}

int FaceMapPatch::group_id()
{
  FaceMapSurf& face_map = ((PatchSurfFaceMap*)bdry())->face_map();
  return face_map.get_group_id(_V);
}

NumVec<OnSurfacePoint> FaceMapPatch::sample_patch(
        int num_samples,
        SamplingPattern sampling_pattern){

    NumVec<OnSurfacePoint> samples(num_samples*num_samples);

    double step = 1./double(num_samples+1);

    for(int j=0; j<num_samples; j++) {
        for(int i=0; i<num_samples; i++) {

            OnSurfacePoint sample;
            switch(sampling_pattern){
                case EQUISPACED:
                    sample.parametric_coordinates = Point2(i*step,j*step);
                    break;
                case CHEBYSHEV:
                    sample.parametric_coordinates = Point2(
                            (cos(i*M_PI/double(num_samples-1))+1.)/2.,
                            (cos(j*M_PI/double(num_samples-1))+1.)/2.);
                    break;
                default:
                    assert(0);
                    break;
            }
            sample.parent_patch = _V;
            sample.region = ON;

        }
    }
    return samples;
}

// ---------------------------------------------------------------------- 
int FaceMapPatch::is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid)
{
  ebiFunctionBegin;
  
  FacePointFaceMap* face_point_blended = (FacePointFaceMap*) face_point;
  double* cd = face_point_blended->cd();
  is_xy_valid(cd, is_valid);

  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 

int FaceMapPatch::face_point_to_xy(FacePointOverlapping* face_point, double* xy)
{
  ebiFunctionBegin;
  FacePointFaceMap* face_point_blended = (FacePointFaceMap*)face_point;  
  double* cd = face_point_blended->cd();
  xy[0] = cd[0];
  xy[1] = cd[1];
  
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
//
int FaceMapPatch::is_xy_valid(double* xy, bool& is_valid)
{
  ebiFunctionBegin;
  // Make sure xy is in [0,1]x[0,1]
  is_valid = xy[0] <= 1. &&
             0. <= xy[0]  &&
             xy[1] <= 1. &&
             0. <= xy[1];

  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
int FaceMapPatch::is_xy_dominant(double* xy, bool& dominant)
{
  ebiFunctionBegin;
  
  // All valid points are "dominant" for face-centered maps
  is_xy_valid(xy, dominant);
  //dominant = true;
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
int FaceMapPatch::xy_to_face_point(double* xy, FacePointOverlapping* face_point)
{
  ebiFunctionBegin;
  double cd[2];		
  cd[0] = xy[0];
  cd[1] = xy[1];

  FacePointFaceMap* tmp = (FacePointFaceMap*)face_point; 
  *tmp = FacePointFaceMap(_V, cd);
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
int FaceMapPatch::xy_to_patch_coords(double* xy, int flag, double* ret)
{
  ebiFunctionBegin;
  FaceMapSurf& face_map = ((PatchSurfFaceMap*)bdry())->face_map();
  
  bool is_valid;
  is_xy_valid(xy, is_valid);
  ebiAssert(is_valid);

  face_map.eval(flag, _V, xy, (Point3*)ret);

  // when evaluating first partial derivatives...
  if(flag & EVAL_FD){
      // note default behavior (-dom = 0, bc_orientation = 1) makes normals
      // pointing out of V
      
      // Let V be volume enclosed by the boundary \Gamma
      // exterior means solving the PDE on the exterior of the domain, i.e.
      // normals need to point into V for the PDE
      bool is_exterior_problem = Options::get_int_from_petsc_opts("-dom") == 1;

      // boundary component is oriented with normals pointing out of V,
      // regardless of PDE
      bool is_oriented_with_exterior_pointing_normals = orientation() == 1;

      bool normal_flip = false;
      if(is_oriented_with_exterior_pointing_normals && !is_exterior_problem){
      //normals are pointing out of V and we're solving an interior problem, all
      //good
      } else if(is_oriented_with_exterior_pointing_normals && is_exterior_problem){
      //normals are pointing out of V and we're solving an exterior problem,
      //need to flip
          normal_flip = true;
      } else if(!is_oriented_with_exterior_pointing_normals && !is_exterior_problem){
      //normals are pointing into V and we're solving an interior problem,
      // this is a multiply connected boundary interior problem. 
      // need to flip
          normal_flip = true;
      } else if(!is_oriented_with_exterior_pointing_normals && is_exterior_problem){
      //normals are pointing into V and we're solving an exterior problem,
      // this is a multiply connected boundary exterior problem. 
      // need to flip
          normal_flip = true;
       }
      if(normal_flip){
      // flip the normals by interchanging derivative order
     Point3* point_ret = (Point3*) ret;
     // flip the normals by interchanging derivative order
     Point3 temp = point_ret[1];
     point_ret[1] = point_ret[2];
     point_ret[2] = temp;
      }
 
  /*if(Options::get_int_from_petsc_opts("dom") == 1 &&// is exterior problem
          flag & EVAL_FD // we're evaluating 1st derivatives
          ){ 
      // flip the normals by interchanging derivative order
      Point3 temp = ret[1];
      ret[1] = ret[2];
      ret[2] = temp;
  }*/
  }
  // MJM NOTE 4/18 this breaks for second derivatives, I think

  ebiFunctionReturn(0);
}
Point3 FaceMapPatch::normal(double* xy) {
    vector<Point3> position_and_derivs(3,Point3(0.));
    xy_to_patch_coords(xy, EVAL_VL|EVAL_FD, (double*)position_and_derivs.data());

    Point3 normal = cross(position_and_derivs[1],position_and_derivs[2]);
    return normal.dir();

}

void FaceMapPatch::eval_unsafe(double* xy, int flag, double* ret) {
  FaceMapSurf& face_map = ((PatchSurfFaceMap*)bdry())->face_map();
  face_map.eval(flag, _V, xy, (Point3*)ret);
}
// ---------------------------------------------------------------------- 
int FaceMapPatch::estimate_jacobian(double* jac)
{
  ebiFunctionBegin;
  double xy[2];  xy[0]=.5;  xy[1]=.5;
  Point3 ret[3];
  FaceMapSurf& face_map = ((PatchSurfFaceMap*)bdry())->face_map();

  face_map.eval(BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, _V, xy, ret); //blendsurf v3

  jac[0] = max(ret[1].l2(), ret[2].l2());
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 

int FaceMapPatch::xy_to_patch_value(double* xy, int flag, double* ret)
{
  ebiFunctionBegin;

  ebiAssert(flag == BdSurf::EVAL_VALUE);

  *ret =1.;

  ebiFunctionReturn(0);
}

void FaceMapPatch::bounding_box(Point3& bounding_box_min, Point3& bounding_box_max){
    FaceMapSurf& face_map = ((PatchSurfFaceMap*)bdry())->face_map();
    face_map.bounding_box(_V, bounding_box_min, bounding_box_max);
    /*
    Vector<Point3> control_points = face_map.control_points(_V);

    bounding_box_min = Point3(DBL_MAX, DBL_MAX, DBL_MAX);
    bounding_box_max = -Point3(DBL_MAX, DBL_MAX, DBL_MAX);
    //bounding_box_max = Point3(DBL_MIN, DBL_MIN, DBL_MIN);

    for(int i = 0; i < control_points.length(); i++){
        Point3 control_point = control_points(i);
        bounding_box_min.x() = min(control_point.x(), bounding_box_min.x());
        bounding_box_min.y() = min(control_point.y(), bounding_box_min.y());
        bounding_box_min.z() = min(control_point.z(), bounding_box_min.z());
        
        bounding_box_max.x() = max(control_point.x(), bounding_box_max.x());
        bounding_box_max.y() = max(control_point.y(), bounding_box_max.y());
        bounding_box_max.z() = max(control_point.z(), bounding_box_max.z());
    }
    */
}

PatchChildrenMap PatchSurfFaceMap::find_subpatches(
        PatchSurfFaceMap* refined_face_map){
    PatchChildrenMap patch_children;
    //assert(0); // TODO BUG subpatch->_parent_id points to the global p4est tree id. 
               //Need yet another index that points to the coarse patch the
               //refined one is derived from.
    for(auto patch : refined_face_map->patches()){
        auto subpatch = dynamic_cast<FaceMapSubPatch*>(patch);
        patch_children[subpatch->_parent_id].push_back(subpatch->V());
    }
    return patch_children;
}
void FaceMapPatch::principal_curvatures(Point2 xy, double& k1, double& k2){
    double H = mean_curvature(xy);
    double K = gaussian_curvature(xy);
    k1 = H + sqrt(H*H-K);
    k2 = H - sqrt(H*H-K);
}

END_EBI_NAMESPACE
