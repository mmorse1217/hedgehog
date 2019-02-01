#ifndef _AGSURF_HPP_
#define _AGSURF_HPP_

#include "common/ebi.hpp"
#include "common/ebiobject.hpp"
#include "vec2t.hpp"
#include "vec3t.hpp"

BEGIN_EBI_NAMESPACE

using std::vector;

class AnalyticSurf: public EbiObject
{
public:
  enum { EVAL_VL=1, EVAL_FD=2, EVAL_SD=4 };
  class Group {
  protected:
	 Point3 _ctr;
	 Point3 _scl;
	 int _orient;
  public:
	 Point3& ctr() { return _ctr; }
	 Point3& scl() { return _scl; }
	 int& orient() { return _orient; }
  };
  
protected:
  vector<Group> _groups;
  
public:
  AnalyticSurf(const string& n, const string& p): EbiObject(n,p) {;}
  ~AnalyticSurf() {;}
  int setup(istream& in);  //int dump( ostream& ot);
  Point3 ctr(); //{ return _sft; }
  void bbx(Point3& bbmin, Point3& bbmax); //{ a = _sft-_scl; b = _sft+_scl; }
  
  int eval(int flags, int ci, double* tg, Point3* ret);
  
  vector<Group>& groups() { return _groups; }

  static void UVW2XYZ(int, double*, double*);
  static void XYZ2UVW(int, double*, double*);
};

//int rotate(const Quaternion& qr);  //int shift(const Point3& sh);

END_EBI_NAMESPACE

#endif
