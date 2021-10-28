#include "patch_surf_analytic.hpp"

BEGIN_EBI_NAMESPACE

using std::istringstream;
using std::abs;
using std::cerr;

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSurfAnalytic::setup"
int PatchSurfAnalytic::setFromOptions()
{
  ebiFunctionBegin;
  PetscBool flg = PETSC_FALSE;
  //filename
  char file[100];  iC( PetscOptionsGetString(NULL, prefix().c_str(), "-filename", file, 100, &flg) ); ebiAssert(flg==PETSC_TRUE);
  _filename = file;
  // Pou := Partition of Unity
  // Default is pouctrl = 3, pdge = 5 => oprtimum blend, degree 5
  cerr << "start" << endl;
  iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-pouctrl",& _pouctrl,  &flg) );
 
  if (flg == false) { _pouctrl = 3; }
  cerr << _pouctrl << " " << _pdeg << endl;
  iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-pdeg",&_pdeg,  &flg) );

  if (flg == false) { _pdeg = 5; }
  cerr << _pouctrl << " " << _pdeg << endl;
  if (_pouctrl == 0 && _pdeg == 0){
	 _pouctrl = 3;
	 _pdeg = 5;
  }
  if (_pdeg == 0 && _pouctrl != 0){ cerr << "BAD POU" << endl; exit(-1); }
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSurfAnalytic::setup"
int PatchSurfAnalytic::setup()
{
  ebiFunctionBegin;

  //1. read in and generate all scbs (Matt: what is scbs?)
  string filedata;  iC( file2string(mpiComm(), _filename.c_str(), filedata) ); //parallel communication
  istringstream fin(filedata);
  iC( _agsurf.setup(fin) ); // pass input to BlendSurf

  //2. derive the patches from the scbs and make a list of patches
  int pcnt = 0;
  for(unsigned si=0; si<_agsurf.groups().size(); si++) {
	 for(int ci=0; ci<6; ci++) {
		_patches.push_back( new AnalyticPatch(this, pcnt, si, ci, _pouctrl, _pdeg) );
		pcnt++;
	 }
  }

  //3. Make lists of centers and orientations of partitions of unity
  for(unsigned si=0; si<_agsurf.groups().size(); si++) {
	 _boundary_component_center.push_back( _agsurf.groups()[si].ctr() );
	 _boundary_orientation.push_back( _agsurf.groups()[si].orient() );
  }
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int PatchSurfAnalytic::face_point_size_in_doubles() {
    return (sizeof(FacePointAnalytic)-1)/sizeof(double)+1;
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSurfAnalytic::patches_containing_face_point"
int PatchSurfAnalytic::patches_containing_face_point(FacePointOverlapping* face_point, vector<int>& pivec)
{
  // For a given boundary point, it returns a vector of patch indices
  // containing that boundary point
  ebiFunctionBegin;
  pivec.clear();
  FacePointAnalytic* bb = (FacePointAnalytic*)face_point;
  AnalyticSurf& agsurf = _agsurf;
  int si = bb->si();
  
  int off = si*6;
  Point3 uni(bb->unipos());
  for(int ci=0; ci<6; ci++) {
	 AnalyticPatch* curpch = (AnalyticPatch*)(_patches[off+ci]);
	 double FM = - curpch->bnd();
	 double TO =   curpch->bnd();
	 double u[3];	 agsurf.XYZ2UVW(ci, uni.array(), u);
	 double tg[2];		tg[0] = atan2(u[0], u[2]);		tg[1] = atan2(u[1], u[2]);
	 if(u[2]>0.0 &&
		 tg[0]>=FM && tg[0]<TO && tg[1]>=FM && tg[1]<TO ) {
		pivec.push_back(ci + off);
	 }
  }
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
// ---------------------------------------------------------------------- 
// ---------------------------------------------------------------------- 
// ---------------------------------------------------------------------- 
double AnalyticPatch::bnd()
{
  return PI/4.0*1.5;
}
int AnalyticPatch::group_id()
{
  return _si;
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::is_face_point_valid"
int AnalyticPatch::is_face_point_valid(FacePointOverlapping* b, bool& is_valid)
{
  // ??
  ebiFunctionBegin;
  double FM = -bnd();
  double TO =  bnd();
  FacePointAnalytic* face_point = (FacePointAnalytic*)b;
  if(face_point->si()!=_si) {
	 is_valid = false;
  } else {
	 Point3 uni( face_point->unipos() );
	 AnalyticSurf& agsurf = ((PatchSurfAnalytic*)bdry())->agsurf();
	 double u[3]; agsurf.XYZ2UVW(_ci, uni.array(), u);
	 double tg[2];  tg[0] = atan2(u[0], u[2]);  tg[1] = atan2(u[1], u[2]);
	 is_valid = ( u[2]>0.0 &&
				tg[0]>=FM && tg[0]<TO && tg[1]>=FM && tg[1]<TO );
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::face_point_to_xy"
int AnalyticPatch::face_point_to_xy(FacePointOverlapping* face_point, double* tg)
{
  // Returns .... (x,y)-coordinates of a boundary point ??
  ebiFunctionBegin;
  
  double FM = -bnd();
  double TO =  bnd();
  FacePointAnalytic* face_point_analytic = (FacePointAnalytic*) face_point;

  ebiAssert( face_point_analytic->si()==_si );
  Point3 uni( face_point_analytic->unipos() );

  AnalyticSurf& agsurf = ((PatchSurfAnalytic*)bdry())->agsurf();
  double u[3]; agsurf.XYZ2UVW(_ci, uni.array(), u);

  tg[0] = atan2(u[0], u[2]);  tg[1] = atan2(u[1], u[2]);
  //should always be good
  //
  ebiAssert( u[2]>0.0 &&
				 tg[0]>=FM && tg[0]<TO && tg[1]>=FM && tg[1]<TO );
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::is_xy_valid"
int AnalyticPatch::is_xy_valid(double* tg, bool& is_valid)
{
  ebiFunctionBegin;
  double FM = -bnd();
  double TO =  bnd();
  is_valid = ( tg[0]>=FM && tg[0]<TO && tg[1]>=FM && tg[1]<TO );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::is_xy_dominant"
int AnalyticPatch::is_xy_dominant(double* tg, bool& dominant)
{
  ebiFunctionBegin;
  double FM = -PI/4;
  double TO =  PI/4;
  dominant = ( tg[0]>=FM && tg[0]<TO && tg[1]>=FM && tg[1]<TO );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::xy_to_face_point"
int AnalyticPatch::xy_to_face_point(double* xy, FacePointOverlapping* face_point)
{
  // Maps (x,y) coordinates to the nearest boundary point??
  ebiFunctionBegin;
  
  double FM = -bnd();
  double TO =  bnd();
  ebiAssert( xy[0]>=FM && xy[0]<TO && xy[1]>=FM && xy[1]<TO );
  AnalyticSurf& agsurf = ((PatchSurfAnalytic*)bdry())->agsurf();
  double theta = xy[0];
  double gamma = xy[1];
  double tt = tan(theta); double tt2 = tt*tt;
  double tg = tan(gamma); double tg2 = tg*tg;
  double r = sqrt(1.0 + tt2 + tg2);
  Point3 u(DIM3);		u[0] = tt / r;		u[1] = tg / r;		u[2] = 1  / r;
  Point3 o(DIM3);		agsurf.UVW2XYZ(_ci, u.array(), o.array()); //the unit position
  FacePointAnalytic* tmp = (FacePointAnalytic*)face_point;
  *tmp = FacePointAnalytic(_si, o);
  
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::xy_to_patch_coords"
int AnalyticPatch::xy_to_patch_coords(double* xy, int flag, double* ret)
{
  ebiFunctionBegin;
  double FM = -bnd();
  double TO =  bnd();
  ebiAssert( xy[0]>=FM && xy[0]<TO && xy[1]>=FM && xy[1]<TO );
  AnalyticSurf& agsurf = ((PatchSurfAnalytic*)bdry())->agsurf();
  iC( agsurf.eval(flag, _si*6 + _ci, xy, (Point3*)ret) );  //AnalyticSurf::EVAL_VL|AnalyticSurf::EVAL_FD
  
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::xy_to_patch_value"
int AnalyticPatch::xy_to_patch_value(double* xy, int flag, double* ret)
{
  ebiFunctionBegin;
  ebiAssert(flag==EVAL_VL);
  ret[0] = POU(flag, xy);
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "AnalyticPatch::estimate_jacobian"
int AnalyticPatch::estimate_jacobian(double* jac)
{
  ebiFunctionBegin;
  double xy[2];  xy[0]=0;  xy[1]=0;
  Point3 ret[3];
  AnalyticSurf& agsurf = ((PatchSurfAnalytic*)bdry())->agsurf();
  iC( agsurf.eval(EVAL_VL|EVAL_FD, _ci, xy, ret) );
  jac[0] = max(ret[1].l2(), ret[2].l2()); //ESTIMATE JACOBIAN
  ebiFunctionReturn(0);
}
// ----------------------------------------------------------------------
// No longer used
double AnalyticPatch::SHAPE(int flag, double* uv) {
  exit(0);
  ebiAssert( flag == EVAL_VL);
  double lt = EX() - BUF();
  double rt = EX() + BUF();  //double lt = PI/4.0 * 0.75; //PI/2.0 - _TO  //double rt = PI/4.0 * 1.25; //_TO
  double t,a,b;  //double poucp[6] = {1,1,1,0,0,0};
  t = (abs(uv[0]) - lt) / (rt-lt);  //if(t<=0) a = 1.0; else if(t>=1) a = 0.0; else a = exp(2.0*exp(-1.0/(t))/(t-1.0));
  if(t<=0) { a = 1.0; }
  else if (t>=1) {a = 0.0; }
  else { a = exp(1-1.0/(1-exp(1-1.0/t))); }
  t = (abs(uv[1]) - lt) / (rt-lt);  //if(t<=0) b = 1.0; else if(t>=1) b = 0.0; else b = exp(2.0*exp(-1.0/(t))/(t-1.0));
  if(t<=0) { b = 1.0; }
  else if(t>=1) {b = 0.0; }
  else { b = exp(1-1.0/(1-exp(1-1.0/t))); }
  return a*b;
}

double AnalyticPatch::POU(int flag, double* thga) {
  double tt = tan(abs(thga[0]));
  double tg = tan(abs(thga[1]));
  double pa[2]; pa[0] = abs(thga[0]); pa[1] = abs(thga[1]);
  double pb[2]; pb[0] = atan2(tg,tt); pb[1] = atan2(1,tt);
  double pc[2]; pc[0] = atan2(tt,tg); pc[1] = atan2(1,tg);

  double resu[3];  
  double resv[3];
  double LB = EX() - BUF();
  double UB = EX() + BUF();

  //BdSurf& bdsurf = ((PatchSurfAnalytic*)bdry())->bdsurf();
  int pouctrl = _pouctrl;
  int pdeg = _pdeg;
  
  ( pou1d(flag, pa[0], LB, UB, resu, pouctrl, pdeg) );
  ( pou1d(flag, pa[1], LB, UB, resv, pouctrl, pdeg) );
  double sa = resu[0]*resv[0];

  ( pou1d(flag, pb[0], LB, UB, resu, pouctrl, pdeg) );
  ( pou1d(flag, pb[1], LB, UB, resv, pouctrl, pdeg) );
  double sb = resu[0]*resv[0];
  
  ( pou1d(flag, pc[0], LB, UB, resu, pouctrl, pdeg) );
  ( pou1d(flag, pc[1], LB, UB, resv, pouctrl, pdeg) );
  double sc = resu[0]*resv[0];
  double tmp = sa/(sa+sb+sc);
  return tmp;
}
/*
int AnalyticPatch::pou1d(int flags, double u, double LB, double UB, double *res){
  ebiFunctionBegin;
  ebiAssert(flags == EVAL_VL);

  double* val = res;
  double* ext = res+1;
  double* fnl = res+2;
  double t = (abs(u) - LB)/(UB-LB);  //if(t<=0) a = 1.0; else if(t>=1) a = 0.0; else a = exp(2.0*exp(-1.0/(t))/(t-1.0));

  //t = ((u)-(1.0/32.0))/(31.0/32.0 - 1.0/32.0);
  BdSurf& bdsurf = ((PatchSurfAnalytic*)bdry())->bdsurf();
  if(t<=0) { val[0] = 1.0; }
  else if (t>=1) {val[0] = 0.0; }
  else { val[0] = exp(1-1.0/(1-exp(1-1.0/t))); }
  
  //cerr << t << " " << val[0] << endl;
  double ERR = 1e-7;
  if(       t<=0+ERR) {
	 val[0] = 1;
  } else if(t>=1-ERR) {
	 val[0]=0;
  } else{
	 if(bdsurf.pouctrl()==0 ) {
		double s = 1-t;
		double t2 = t*t;		double t3 = t2*t;		double t4 = t2*t2;
		double s2 = s*s;		double s3 = s2*s;		double s4 = s2*s2;
		double a = 2*exp(-1/t)/(t-1);
		double b = 2*exp(-1/s)/(s-1);
		double da =  a*(1/(t*t) - 1/(t-1));
		double db = -b*(1/(s*s) - 1/(s-1));
		double dda = a*(-4*t3+7*t2-4*t+1+2*t4)/((t-1)*(t-1))/t4;
		double ddb = b*(-4*s3+7*s2-4*s+1+2*s4)/((s-1)*(s-1))/s4;
		double ea = exp(a);
		double eb = exp(b);
		double f = ea;
		double g = ea+eb;
		double df = ea*da;
		double dg = ea*da + eb*db;
		double ddf = ea*da*da + ea*dda;
		double ddg = ea*da*da + ea*dda + eb*db*db + eb*ddb;
		if(flags & EVAL_VL) {
		  val[0] = f/g;
		}
		else { ebiAssert(0); } //For now 
	 }
  }
  //cerr << t << " " << val[0] << endl;

  ebiFunctionReturn(0);
}
*/
END_EBI_NAMESPACE
		
