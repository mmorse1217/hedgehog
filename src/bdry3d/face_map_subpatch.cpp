#include "face_map_subpatch.hpp"
#include <sampling.hpp>
#include <common/interpolate.hpp>
#include <common/vtk_writer.hpp>
#include "bie3d/error_estimate.hpp"
using Sampling::sample_2d;
using Sampling::chebyshev1;
using Sampling::chebyshev2;
BEGIN_EBI_NAMESPACE


FaceMapSubPatch::FaceMapSubPatch(FaceMapPatch* face_map_patch, 
        Interval x_interval, 
        Interval y_interval,
        int level,
        int id,
        int parent_id,
        DblNumVec* quad_weights):
    Patch(face_map_patch->bdry(), face_map_patch->pi()),
    _face_map_patch(face_map_patch),
    _x_interval(x_interval),
    _y_interval(y_interval),
    _id(id),
    _level(level),
    _parent_id(parent_id),
    _quadrature_weights(quad_weights)
{
    estimate_jacobian(&_approx_jacobian);

}
void FaceMapSubPatch::compute_near_zone_distance(){
    _near_zone_distance = ErrorEstimate::evaluate_near_zone_distance(this);    
}

Point3 FaceMapSubPatch::normal(double* xy) {
    double xy_scaled[2];
    rescale(xy, xy_scaled);
    return _face_map_patch->normal(xy_scaled);
}


// Pass virtual function calls to the FaceMapPatch containing the
// FaceMapSubPatch

double FaceMapSubPatch::bnd(){
    return _face_map_patch->bnd();
}

int FaceMapSubPatch::group_id(){
    return _face_map_patch->group_id();
}

int FaceMapSubPatch::is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid){
    return _face_map_patch->is_face_point_valid(face_point, is_valid);
}

int FaceMapSubPatch::face_point_to_xy(FacePointOverlapping* face_point, double* xy){
    //double xy_scaled[2];
    //rescale(xy, xy_scaled);
    return _face_map_patch->face_point_to_xy(face_point, xy);
}
// SubPatches are defined on [0,1]^2, but are actually contained in the parent
// face-map with domain [0,1]^2. 
// Scale the point in  [0,1]^2 for the subpatch to _x_interval x _y_interval for
// the parent face-map
void FaceMapSubPatch::rescale(double* xy, double* xy_scaled){
    xy_scaled[0] = (_x_interval.second - _x_interval.first)*xy[0] + _x_interval.first;
    xy_scaled[1] = (_y_interval.second - _y_interval.first)*xy[1] + _y_interval.first;
}
void FaceMapSubPatch::unscale(double* xy, double* xy_scaled){
    xy_scaled[0] = 1./(_x_interval.second - _x_interval.first)*(xy[0] - _x_interval.first);
    xy_scaled[1] = 1./(_y_interval.second - _y_interval.first)*(xy[1] - _y_interval.first);
}

// all points in patch domain [0,1]^2 are valid
int FaceMapSubPatch::is_xy_valid( double* xy, bool& is_valid){
    //double xy_scaled[2];
    //rescale(xy, xy_scaled);
    
    // pass to original implementation
    //return _face_map_patch->is_xy_valid(xy, is_valid);
    Point2 p(xy);
    is_valid = p.x() >= 0. &&
           p.y() >= 0. &&
           p.x() <= 1. &&
           p.y() <= 1.;
}

// all points in patch domain [0,1]^2 are dominant
int FaceMapSubPatch::is_xy_dominant( double* xy, bool& dominant){
    double xy_scaled[2];
    rescale(xy, xy_scaled);
    
    // pass to original implementation
    return _face_map_patch->is_xy_dominant(xy_scaled, dominant);
}


// TODO possible bug might need to change F and cd for the face point
int FaceMapSubPatch::xy_to_face_point(double* xy, FacePointOverlapping* face_point){
    double xy_scaled[2];
    rescale(xy, xy_scaled);
    return _face_map_patch->xy_to_face_point(xy_scaled, face_point);
}


int FaceMapSubPatch::xy_to_patch_coords(double* xy, int flag, double* ret){
    double xy_scaled[2];
    rescale(xy, xy_scaled);
    return _face_map_patch->xy_to_patch_coords(xy_scaled, flag, ret);
}

void FaceMapSubPatch::eval_unsafe(double* xy, int flag, double* ret){
    double xy_scaled[2];
    rescale(xy, xy_scaled);
    _face_map_patch->eval_unsafe(xy_scaled, flag, ret);
}

int FaceMapSubPatch::xy_to_patch_value(double* xy, int flag, double* ret){
    double xy_scaled[2];
    rescale(xy, xy_scaled);
    return _face_map_patch-> xy_to_patch_value(xy_scaled, flag, ret);
}

// TODO definitely wrong (check for correctness)
int FaceMapSubPatch::estimate_jacobian(double* jac){
    double xy_scaled[2];
    double xy[2];
    xy[0] = .5;
    xy[1] = .5;
    rescale(xy, xy_scaled);
  Point3 ret[3];
  FaceMapSurf& face_map = ((PatchSurfFaceMap*)_face_map_patch->bdry())->face_map();

  face_map.eval(BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, _face_map_patch->V(), xy_scaled, ret); 

  jac[0] = max(ret[1].l2(), ret[2].l2());
  return 0;
}

void FaceMapSubPatch::inflated_bounding_box(Point3& bounding_box_min, Point3& bounding_box_max){
    bounding_box(bounding_box_min, bounding_box_max);
    // Increase by the size of the near zone
    //double inflation_factor = error_estimate(n, target_accuracy);
    double target_accuracy = Options::get_double_from_petsc_opts("-target_accuracy");
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    int n = int(floor(1./spacing))+1;
    double mean_curvature = this->mean_curvature(Point2(.5, .5));
    double L = this->_characteristic_length;
    
    double inflation_factor = _near_zone_distance;
    
    bounding_box_max += Point3(inflation_factor);
    bounding_box_min -= Point3(inflation_factor);
}

void FaceMapSubPatch::bounding_box(Point3& bounding_box_min, Point3& bounding_box_max){
  FaceMapSurf& face_map = ((PatchSurfFaceMap*)_face_map_patch->bdry())->face_map();
  if(_level == 0){
      _face_map_patch->bounding_box(bounding_box_min, bounding_box_max);
  } else {
      face_map.bounding_box_on_subdomain(
              _face_map_patch->V(), 
              _x_interval,
              _y_interval, 
              bounding_box_min, 
              bounding_box_max);
  }
}

DblNumMat FaceMapSubPatch::control_points(){
  FaceMapSurf& face_map = ((PatchSurfFaceMap*)_face_map_patch->bdry())->face_map();

    vector<Point3> control_points =  
        face_map.control_points_on_subdomain(
          _face_map_patch->V(), 
          _x_interval,
          _y_interval);
    DblNumMat control_points_ret(DIM, control_points.size());
    for (int i = 0; i < control_points.size(); i++) {
        for (int d = 0; d < DIM; d++) {
           control_points_ret(d,i) = control_points[i](d); 
        }
    }
    return control_points_ret;
}


FaceMapSubPatch* PatchSurfFaceMap::subpatch(int pi){
    assert(pi >= 0);
    assert(pi < _patches.size());
    FaceMapSubPatch* p = dynamic_cast<FaceMapSubPatch*>(_patches[pi]);
    assert(p);
    return p;

}

void FaceMapSubPatch::single_bounding_box_triangle(DblNumMat& vertices, IntNumMat& triangles){
    // make a single triangle that has an edge along the diagonal between the
    // bounding box min and max; the last corner should be irrelevant.
    vertices.resize(3,3); 
    triangles.resize(3,1); 
    Point3 bbox_min;
    Point3 bbox_max;
    inflated_bounding_box(bbox_min, bbox_max);

    for (int d = 0; d < DIM; d++) {
        vertices(d,0) = bbox_min(d); 
        vertices(d,1) = bbox_max(d); 
    }

    vertices(0,2) = bbox_min(0); 
    vertices(1,2) = bbox_min(1); 
    vertices(2,2) = bbox_max(2); 
    
    triangles(0,0) = 0;
    triangles(1,0) = 1;
    triangles(2,0) = 2;

}

void FaceMapSubPatch::mesh_bounding_box(DblNumMat& vertices, IntNumMat& triangles){
    // find the bounding box
    Point3 bbox_min;
    Point3 bbox_max;
    bounding_box(bbox_min, bbox_max);
    
    // enumerate the 8 vertex corners of the box
    vertices.resize(DIM,8);
    for (int i = 0; i < 2; i++) { // 
        for (int j = 0; j < 2; j++) { // y
            for (int k = 0; k < 2; k++) { // z
                int index = 4*i+2*j+k;
                vertices(0,index) = i ? bbox_min(0) : bbox_max(0);
                vertices(1,index) = j ? bbox_min(1) : bbox_max(1);
                vertices(2,index) = k ? bbox_min(2) : bbox_max(2);
            }
        }
    }

    // TODO refactor 3 and 12 into const variables
    // 3 is num verts per triangle, 12 is num triangles per bounding box (6
    // faces with 2 triangles each)
    triangles.resize(3,12);

    // explicitly list indices that make up quad faces of bbox as enumerated
    // above
    vector<vector<int> > bbox_quads(6);
    bbox_quads[0]= {0,2,3,1};
    bbox_quads[1]= {1,3,7,5};
    bbox_quads[2]= {4,0,1,5};
    bbox_quads[3]= {4,6,2,0};
    bbox_quads[4]= {5,7,6,4};
    bbox_quads[5]= {2,6,7,3};
    for (int i = 0; i < 6; i++) {
        const auto& quad = bbox_quads[i];

        // index set for the two triangles that make up each quad
        vector<int> t1 = {0,1,2};
        vector<int> t2 = {0,2,3};
        for (int d = 0; d < 3; d++) {
            // copy them into the output array
            triangles(d,2*i) = quad[t1[d]];
            triangles(d,2*i+1) = quad[t2[d]];
        }
    }

}

bool FaceMapSubPatch::is_patch_valid(Vec function_values_at_parent, 
        Vec function_values_at_children, 
        Vec uv_coordinates_single_patch, int range_dim, double eps){

    DblNumMat function_values_at_children_local(range_dim, function_values_at_children);
    DblNumMat function_values_at_parent_local(range_dim, function_values_at_parent);
    DblNumMat uv_coordinates_single_patch_local(2,uv_coordinates_single_patch);

    int num_samples_per_patch_1d = int(sqrt(uv_coordinates_single_patch_local.n()/4));
    
    DblNumMat interp_nodes= 
        DblNumMat(2, num_samples_per_patch_1d*num_samples_per_patch_1d, true, 
                sample_2d<chebyshev1>(num_samples_per_patch_1d, Sampling::base_domain).data());
    DblNumMat interpolated_values(range_dim, uv_coordinates_single_patch_local.n());
    Interpolate::evaluate_barycentric_interpolant_2d(range_dim, 
            interp_nodes, 
            function_values_at_parent_local, 
            num_samples_per_patch_1d,
            uv_coordinates_single_patch_local, 
            interpolated_values);
    Vec error;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, interpolated_values.n()*interpolated_values.m(), error);
    DblNumMat error_local(interpolated_values.m(), error);
    bool is_valid = true;
    double error_abs_err = 0.;
    for (int i = 0; i < interpolated_values.n(); i++) {
        for (int d = 0; d < interpolated_values.m(); d++) {
            double err = fabs(interpolated_values(d,i) - 
                    function_values_at_children_local(d,i));
            bool abs_err_tolerance = err < eps;

            is_valid = is_valid && abs_err_tolerance;
            error_local(d,i) = err;
            error_abs_err += err*err;
        }
    }
    VecDestroy(&error);
    return is_valid;
}
double FaceMapSubPatch::gaussian_curvature(Point2 xy){
    Point2 xy_scaled;
    rescale(xy.array(), xy_scaled.array());
    return _face_map_patch->gaussian_curvature(xy_scaled);
}
double FaceMapSubPatch::mean_curvature(Point2 xy){
    Point2 xy_scaled;
    rescale(xy.array(), xy_scaled.array());
    return _face_map_patch->mean_curvature(xy_scaled);
}
void FaceMapSubPatch::principal_curvatures(Point2 xy, double& k1, double& k2){
    Point2 xy_scaled;
    rescale(xy.array(), xy_scaled.array());
    return _face_map_patch->principal_curvatures(xy_scaled,k1, k2);
}

END_EBI_NAMESPACE
