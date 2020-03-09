#ifndef _BD3DOV_HPP_
#define _BD3DOV_HPP_

#include "common/ebiobject.hpp"
#include "vec3t.hpp"
#include "vec2t.hpp"
#include "common/nummat.hpp"
#include "common/utils.hpp"
#include <sampling.hpp>

BEGIN_EBI_NAMESPACE

using std::vector;

class Scb3dOv;
class Patch; // Patch on the boundary
class FacePointOverlapping; // Boundary point?

enum SamplingPattern {
    EQUISPACED = 0,
    CHEBYSHEV = 1
};
//----------------------------------------------------------
// 3d Boundary representation
class PatchSurf: public EbiObject
{
protected:
  vector<Patch*> _patches;
  vector<Point3> _boundary_component_center; //center of components
  vector<int> _boundary_orientation; //orientation

  int64_t _pouctrl; 
  int64_t _pdeg;

public:
  PatchSurf(const string& n, const string& p): EbiObject(n,p) {;}
  virtual ~PatchSurf();

  //accessor methods
  int64_t& pouctrl() { return _pouctrl; }
  int64_t& pdeg() { return _pdeg; }

  // list of overlapping patches
  vector<Patch*>& patches() { return _patches; } 

  // list of points in each multiply connected domain (one point inside each domain)
  vector<Point3>& boundary_component_center() { return _boundary_component_center; } 
  
  // orientation of points in each domain
  vector<int>& boundary_orientation() { return _boundary_orientation; }

  //setup ..
  virtual int setFromOptions()=0;
  virtual int setup()=0;
  virtual Point3 ctr()=0;
  virtual void bounding_box(Point3& bbmin, Point3& bbmax)=0;

  //others
  virtual int patches_containing_face_point(FacePointOverlapping*, 
          vector<int>& pivec)=0; //from face_point to all the patches containing it 

  virtual int face_point_size_in_doubles()=0;
  int num_patches(){
    return _patches.size();
  }
  Patch* patch(int pi){
      assert(pi >= 0);
      assert(pi < num_patches());
      return _patches[pi];
  }
};

//----------------------------------------------------------
class Patch {
public:
  enum {	 
      EVAL_VL = 1,
      EVAL_FD = 2  
  };
protected:
  PatchSurf* _bdry;
  int _pi;
    double _characteristic_length;
public:
  double _approx_jacobian;
  Patch(PatchSurf* b, int pi): _bdry(b), _pi(pi) {;}
        double& characteristic_length(){
            return _characteristic_length;
        }

        double characteristic_length(DblNumMat normals, DblNumVec quadrature_weight);
        double compute_surface_area(DblNumMat normal, DblNumVec quadrature_weight);

  virtual ~Patch() {;}
  //...
  virtual double bnd()=0;
  virtual int group_id()=0;
  virtual int is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid)=0;
  virtual int face_point_to_xy(FacePointOverlapping* face_point, double* xy)=0;
  virtual int is_xy_valid( double* xy, bool& is_valid)=0;
  virtual int is_xy_dominant( double* xy, bool& dominant)=0;
  virtual int xy_to_face_point(double* xy, FacePointOverlapping* face_point) = 0;
  virtual int xy_to_patch_coords(double* xy, int flag, double*) = 0;
  virtual int xy_to_patch_value(double* xy, int flag, double*) = 0;
  virtual int estimate_jacobian(double*)=0; //estimate jacobian
  NumVec<Point2> sample_patch(SamplingPattern sampling_pattern);
  DblNumMat generate_qbkix_samples(NumVec<Point2> on_surface_patch_samples);
//bool is_patch_valid(void (*func)(Vec, int, Vec&), int range_dim, double eps);
  
  // return triangle mesh of samples of the patch.
  // The number of vertices in one dimension is n = (floor(1./spacing) + 1) (sorry future me)
  // func returns N = n*n vertices total and 2*(n-1)^2 triangle faces
  // equispaced sampling
  void mesh_patch(double spacing,  // sampling spacing
          Rectangle sample_domain,
          DblNumMat& vertices, //3 x num_vertices matrix of positions
          IntNumMat& faces //3 x num_faces matrix of positions
          );
  // same as above, vertices are provided, and assumed to be the correct amount
  void mesh_patch_no_sampling(  
          double spacing, //num vertices in 1d along a regular sampling
          IntNumMat& faces //3 x num_faces matrix of positions
          );
  // same as above, vertices are provided, and assumed to be the correct amount
  int mesh_num_vertices_1d(double spacing){
    return int(floor(1./spacing))+1;
  }
  int mesh_num_vertices(double spacing){
      int n = mesh_num_vertices_1d(spacing);
      return n*n;
  }
  int mesh_num_triangles(double spacing){
      int n = mesh_num_vertices_1d(spacing);
      return 2*(n-1)*(n-1);
  }


  //...
  PatchSurf* bdry() { return _bdry; }
  int pi() { return _pi; }
};
//----------------------------------------------------------
class FacePointOverlapping //unique identification of point on surface, non-overlapping
{
public:
  FacePointOverlapping() {;}
  virtual ~FacePointOverlapping() {;}  //virtual int cid() = 0;
};

END_EBI_NAMESPACE

#endif
