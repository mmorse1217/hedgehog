#include "analytic_surf.hpp"
#include "common/vecmatop.hpp"

BEGIN_EBI_NAMESPACE

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticSurf::setup"
int AnalyticSurf::setup(istream& fin)
{
  ebiFunctionBegin;
  
  string tmp;
  int cnt;
  fin>>tmp; ebiAssert(tmp=="groups");
  fin>>cnt;
  
  _groups.resize(cnt);
  double x,y,z;
  int t;
  for(int i=0; i<cnt; i++) {
	 fin>>tmp; ebiAssert(tmp=="ctr");
	 fin>>x>>y>>z; _groups[i].ctr() = Point3(x,y,z);
	 fin>>tmp; ebiAssert(tmp=="scl");
	 fin>>x>>y>>z; _groups[i].scl() = Point3(x,y,z);
	 fin>>tmp; ebiAssert(tmp=="orient");
	 fin>>t; _groups[i].orient() = t;
  }
  
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
Point3 AnalyticSurf::ctr()
{
  Point3 a,b; bbx(a,b);
  return 0.5*(a+b);
}

void AnalyticSurf::bbx(Point3& bbmin, Point3& bbmax)
{
  bbmin = Point3( SCL_MAX);
  bbmax = Point3(-SCL_MAX);
  for(unsigned gi=0; gi<_groups.size(); gi++) {
	 Group& g = _groups[gi];
	 bbmin = min(bbmin, g.ctr()-g.scl());
	 bbmax = max(bbmax, g.ctr()+g.scl());
  }
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticSurf::eval"
int AnalyticSurf::eval(int flags, int idx, double* pt, Point3* ret)
{
  ebiFunctionBegin;
  //6 faces for each group
  int gi = idx/6;
  int ci = idx%6;
  Point3 ctr = _groups[gi].ctr();
  Point3 scl = _groups[gi].scl();
  int orient = _groups[gi].orient();
  scl = double(orient)*scl; //LEX VERY IMPORTANT
  
  ebiAssert( (flags==EVAL_VL) || (flags==(EVAL_VL|EVAL_FD)) );
  double theta = pt[0];
  double gamma = pt[1];
  double tt = tan(theta);		double tt2 = tt*tt;		//double tt3 = tt2*tt;
  double tg = tan(gamma);		double tg2 = tg*tg;		//double tg3 = tg2*tg;
  double r2 = 1.0 + tt2 + tg2;		double r = sqrt(r2);		double r3 = r2*r;  //double r4 = r2*r2;    
  Point3 u;
  u[0] = tt / r;
  u[1] = tg / r;
  u[2] = 1  / r;
  Point3 ut;
  ut[0] = (1+tt2)*(1+tg2)/r3;
  ut[1] = -tt*tg*(1+tt2)/r3;
  ut[2] = -tt*(1+tt2)/r3;
  Point3 ug;
  ug[0] = -tt*tg*(1+tg2)/r3;
  ug[1] = (1+tt2)*(1+tg2)/r3;
  ug[2] = -tg*(1+tg2)/r3;
  
  Point3 o;  UVW2XYZ(ci, u.array(),  o.array());  o  = ewmul(o, scl);
  Point3 ot; UVW2XYZ(ci, ut.array(), ot.array()); ot = ewmul(ot,scl);
  Point3 og; UVW2XYZ(ci, ug.array(), og.array()); og = ewmul(og,scl);
  //Point3 x;  DblNumVec ov( 3,false,o.array());   DblNumVec xv( 3,false,x.array());   iC( dgemv(1.0, _rot, ov,  0.0, xv) );
  //Point3 xt; DblNumVec otv(3,false,ot.array());  DblNumVec xtv(3,false,xt.array());  iC( dgemv(1.0, _rot, otv, 0.0, xtv) );
  //Point3 xg; DblNumVec ogv(3,false,og.array());  DblNumVec xgv(3,false,xg.array());  iC( dgemv(1.0, _rot, ogv, 0.0, xgv) );
  o = o + ctr;
  
  //Point3 xn = cross(xt,xg);
  //double len = xn.l2();
  //xn = xn/len;
  //if(orient==-1) og = -og; ///flip xg's sign if necessary
  if(flags & EVAL_VL) {	 ret[0] = o;  }
  if(flags & EVAL_FD) {	 ret[1] = ot;	 ret[2] = og;  }  
  ebiFunctionReturn(0);
}

void AnalyticSurf::UVW2XYZ(int ci, double* u, double* x) {
  switch(ci) {
  case 0: x[0] = u[0]; x[1] = u[1]; x[2] = u[2]; break;
  case 1: x[0] =-u[1]; x[1] =-u[0]; x[2] =-u[2]; break;
  case 2: x[0] = u[2]; x[1] = u[0]; x[2] = u[1]; break;
  case 3: x[0] =-u[2]; x[1] =-u[1]; x[2] =-u[0]; break;
  case 4: x[0] = u[1]; x[1] = u[2]; x[2] = u[0]; break;
  case 5: x[0] =-u[0]; x[1] =-u[2]; x[2] =-u[1]; break;
  }
  
}

void AnalyticSurf::XYZ2UVW(int ci, double* x, double* u) {
  switch(ci) {
  case 0: u[0] = x[0]; u[1] = x[1]; u[2] = x[2]; break;
  case 1: u[0] =-x[1]; u[1] =-x[0]; u[2] =-x[2]; break;
  case 2: u[0] = x[1]; u[1] = x[2]; u[2] = x[0]; break;
  case 3: u[0] =-x[2]; u[1] =-x[1]; u[2] =-x[0]; break;
  case 4: u[0] = x[2]; u[1] = x[0]; u[2] = x[1]; break;
  case 5: u[0] =-x[0]; u[1] =-x[2]; u[2] =-x[1]; break;
  }
}


END_EBI_NAMESPACE
