#ifndef _BD3DAG_HPP_
#define _BD3DAG_HPP_

#include <bdsurf.hpp>
#include "patch_surf.hpp"
#include "analytic_surf.hpp"
#include "vec2t.hpp"
#include "common/vecmatop.hpp"

BEGIN_EBI_NAMESPACE

//----------------------------------------------------------
class PatchSurfAnalytic: public PatchSurf
{
protected:
  //PARAMS
  //COMPONENTS
  AnalyticSurf _agsurf;
  BdSurf _bdsurf; // This probably shouldn't be here

public:
  string _filename;
  PatchSurfAnalytic(const string& n, const string& p):
	 PatchSurf(n,p), _agsurf(n+"AGSURF_",p+"agsurf_"), _bdsurf() {;}
  ~PatchSurfAnalytic() {;}


  //accessors 
  string& filename() { return _filename; }
  AnalyticSurf& agsurf() { return _agsurf; }
  BdSurf& bdsurf() { return _bdsurf; }
  //setup ..
  int setFromOptions();
  int setup();
  Point3 ctr() { return _agsurf.ctr(); } // Matt: center of the surface in what sense?
  void bounding_box(Point3& bbmin, Point3& bbmax) { _agsurf.bbx(bbmin, bbmax); }
  int patches_containing_face_point(FacePointOverlapping*, vector<int>& pivec); //from face_point to all the patches containing it 
  int face_point_size_in_doubles();
};

//----------------------------------------------------------
class AnalyticPatch: public Patch
{
protected:
  int _si;
  int _ci; //camera id
  int _pouctrl; 
  int _pdeg; 
public:
  AnalyticPatch(PatchSurf* b, int pi, int si, int ci, int po, int pd): Patch(b, pi), _si(si), _ci(ci), _pouctrl(po), _pdeg(pd) {;}
  int& si() { return _si; }
  int& ci() { return _ci; }
  
  double bnd();
  int group_id();
  int is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid); //whether face_point in this patch or not
  int face_point_to_xy(FacePointOverlapping* face_point, double* xy);
  int is_xy_valid( double* xy, bool& is_valid);
  int is_xy_dominant( double* xy, bool& dominant);
  int xy_to_face_point(double* xy, FacePointOverlapping* face_point);
  int xy_to_patch_coords(double* xy, int flag, double*);
  int xy_to_patch_value(double* xy, int flag, double*);
  int estimate_jacobian(double*);

  int& pouctrl() { return _pouctrl; }
  int& pdeg() { return _pdeg; }
protected:
  double POU(int flag,   double* tg);
  double SHAPE(int flag, double* uv);
  double EX()  { return PI/4.0; }
  double BUF() { return PI/4.0*0.5; }

};

//----------------------------------------------------------
class FacePointAnalytic: public FacePointOverlapping
{
protected:
  int _si; //which sphere 
  Point3 _unipos; //pt on the sphere
public:
  FacePointAnalytic(int si, Point3 unipos): _si(si), _unipos(unipos) {
  }
  int& si() { return _si; }
  Point3& unipos() { return _unipos; }
};

END_EBI_NAMESPACE

#endif
