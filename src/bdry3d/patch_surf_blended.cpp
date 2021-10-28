#include "patch_surf_blended.hpp"
#include <evalfromw.hpp>
#include "common/stats.hpp"

using std::istringstream;
BEGIN_EBI_NAMESPACE

//double BlendedPatch::_UB = 0.75; 
// ---------------------------------------------------------------------- 
PatchSurfBlended::PatchSurfBlended(const string& n, const string& p):
  PatchSurf(n,p), _bdsurf(){;} //blendsurf v3
 // _bdsurf(n+"BDSURF_",p+"bdsurf_"){;}
 //blendsurf v1
// ---------------------------------------------------------------------- 
PatchSurfBlended::~PatchSurfBlended()
{;}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSurfBlended::setFromOptions"
int PatchSurfBlended::setFromOptions()
{
  // Load parameters from options file
  ebiFunctionBegin;
  PetscBool flg = PETSC_FALSE;
  //filename
  char file[100];  iC( PetscOptionsGetString(NULL, prefix().c_str(), "-filename", file, 100, &flg) ); ebiAssert(flg==PETSC_TRUE);
  _filename = file;
  
  // Pou := partition of unity
  // Default is pouctrl = 3, pdge = 5 => oprtimum blend, degree 5
  iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-pouctrl", &_pouctrl,  &flg) );
 
  if (flg == false) { _pouctrl = 3; }
  iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-pdeg", &_pdeg,  &flg) );


  if (flg == false) { _pdeg = 5; }

  cerr << _pouctrl << " " << _pdeg << endl;
    char fileb[100];
    iC( PetscOptionsGetString(NULL, prefix().c_str(), "-meshfile", fileb, 100, &flg) );
    ebiAssert(flg==PETSC_TRUE);
    ifstream tmpb(fileb);

    char filea[100];
    iC( PetscOptionsGetString(NULL, prefix().c_str(),
                             "-bdsurf_submatlibfile", filea, 100, &flg) );
    ebiAssert(flg==PETSC_TRUE);
    ifstream tmpa(filea);
    CCSubMatLib* submatlib = new CCSubMatLib();
    cerr<<"submatlib setup"<<endl;
    submatlib->setup(tmpa);

    char filebd[100];
    iC( PetscOptionsGetString(NULL,prefix().c_str(),
                             "-bdsurf_bdulibfile", filebd, 100, &flg) );
    ebiAssert(flg==PETSC_TRUE);
    ifstream tmpbd(filebd);
    BdULib* bdulib = new BdULib();
    cerr<<"bdulib setup"<<endl;
    bdulib->setup(tmpbd);

    int64_t cl;
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_ctrllvl",&cl,  &flg) );
    ebiAssert(flg==PETSC_TRUE);
    assert(cl>=0 && cl<=2);

    int64_t pc;
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_pouctrl",&pc,  &flg) );
    ebiAssert(flg==PETSC_TRUE);
    assert(pc==0 || pc==1 || pc==2);

    int64_t ct;     
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_chttyp",&ct,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    int64_t bt;     
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_bsstyp",&bt,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    int64_t pp;	    
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_stpppg",&pp,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    double lb;
    PetscOptionsGetReal(NULL,prefix().c_str(), "-bdsurf_lb", &lb, &flg);
    ebiAssert(flg==PETSC_TRUE);

    double ub;	
    PetscOptionsGetReal(NULL,prefix().c_str(), "-bdsurf_ub", &ub, &flg);
    ebiAssert(flg==PETSC_TRUE);
    
    int64_t gl;	    
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_indepboun",&gl,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    double fs;	
    PetscOptionsGetReal(NULL,prefix().c_str(), "-bdsurf_flats", &fs, &flg);
    ebiAssert(flg==PETSC_TRUE);

    int64_t cot;    
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_concat",&cot,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    int64_t sd;	   
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_spdeg",&sd,  &flg) );
    ebiAssert(flg==PETSC_TRUE);


    int64_t poudeg; 
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_poubsdeg",&poudeg,  &flg) );
    ebiAssert(flg==PETSC_TRUE);


    int64_t loadW;  
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_loadw", &loadW,  &flg) );
    ebiAssert(flg==PETSC_TRUE);
    //TODO: check why this is a sane default value
    //if(mi==opts.end())
    //    loadW = 0;
    //else {
    //    istringstream ss((*mi).second); ss>>loadW; 
    //}
    int64_t tmp;
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_g", &tmp,  &flg) );
    NUM_PATCHES = tmp;
    ebiAssert(flg==PETSC_TRUE);

    char matrix_dir[100];
    iC( PetscOptionsGetString(NULL,prefix().c_str(), "-bdsurf_matdir", matrix_dir, 100,&flg) );
    ebiAssert(flg==PETSC_TRUE);

    MatricesDir = matrix_dir;

    int64_t activevert = -1;  
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_activevert",&activevert,  &flg) );
    ebiAssert(flg==PETSC_TRUE);
    int64_t refinement_factor= 0;  
    iC( PetscOptionsGetInt(NULL, prefix().c_str(),  "-bdsurf_refinement_factor",&refinement_factor,  &flg) );
    ebiAssert(flg==PETSC_TRUE);

    double spacing = Options::get_double_from_petsc_opts("-bdsurf_interpolant_spacing");
    int interpolate = Options::get_int_from_petsc_opts("-bdsurf_interpolate");
    //BdSurf& bdsurf = new BdSurf(  );

    _bdsurf.setParams(cl,pc,ct,bt,
            pp,lb,ub, gl,
            fs, cot, sd ,poudeg,
            spacing, interpolate,0,0);//,0, loadW );

    _bdsurf.setup(tmpb,0,Point3(1.0),(void*)submatlib,(void*)bdulib); // blendsurf v3
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSurfBlended::setup"
int PatchSurfBlended::setup()
{
  ebiFunctionBegin;
  //1. read in and generate all scbs (??)
  string filedata;  iC( file2string(mpiComm(), _filename.c_str(), filedata) ); 
  cerr << _filename.c_str() << endl;
  istringstream fin(filedata);

    //2. computer patches from mesh
    int pcnt = 0;
    for(int V=0; V<_bdsurf.numVertices(); V++) {
       _patches.push_back(  new BlendedPatch(this, pcnt, V) );
       pcnt++;
    }
    //3. centers
    GpMesh& gpmesh = _bdsurf.gpmesh();
    
    _boundary_component_center = gpmesh.get_interior_points();
    _boundary_orientation  = gpmesh.get_intpt_orientation();

    ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int PatchSurfBlended::face_point_size_in_doubles() {
    return (sizeof(FacePointBlended)-1)/sizeof(double)+1;
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "PatchSurfBlended::patches_containing_face_point"
int PatchSurfBlended::patches_containing_face_point(FacePointOverlapping* face_point, vector<int>& pivec)
{
  ebiFunctionBegin;
  pivec.clear();
  FacePointBlended* bb = (FacePointBlended*) face_point;
  BdSurf& bdsurf = _bdsurf;
  for(int v=0; v<4; v++) {
	 int V;	 int f;
	 double cd[2];
     iC( bdsurf.Fcd2Vfcd(BdSurf::EVAL_VALUE, bb->F(), bb->cd(), v, V, f, cd) );
	 if(cd[0]<bdsurf.EVAL_UB() && cd[1]<bdsurf.EVAL_UB())
		pivec.push_back(V);
  }
  ebiFunctionReturn(0);
}

// Method below implement requirements from Bd3dov.hpp. Functionality should be
// similar to that of agsurf.cpp, except for blended surfaces instead of
// analytic ones.

#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::BlendedPatch"
BlendedPatch::BlendedPatch(PatchSurf* b, int pi, int V): Patch(b, pi), _V(V)
{
  ebiFunctionBegin;
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  iC( bdsurf.paramBound(_V, bdsurf.EVAL_UB(), _bnd) );   
  ebiFunctionReturnVoid;
}

int BlendedPatch::group_id()
{
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  GpMesh& gpmesh = bdsurf.gpmesh();
  return gpmesh.get_group_ids()[_V];
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::is_face_point_valid"
int BlendedPatch::is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid)
{
  ebiFunctionBegin;
  FacePointBlended* face_point_blended = (FacePointBlended*) face_point;
  int F = face_point_blended->F();
  double* cd = face_point_blended->cd();
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  bool found = false;
  int tV, tf;
  double st[6];
  for(int v=0; v<4; v++) {
	 iC( bdsurf.Fcd2Vfcd(BdSurf::EVAL_VALUE, F, cd, v, tV, tf, st) );
	 if(tV==_V) {
		found = true;
		break;
	 }
  }
  if (!found) { cerr << "is_face_point_valid error" << endl; exit(-1); }
  is_valid = (st[0]<bdsurf.EVAL_UB() && st[1]<bdsurf.EVAL_UB());  //is_valid = 1;  else	 is_valid = 0;
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::face_point_to_xy"
int BlendedPatch::face_point_to_xy(FacePointOverlapping* face_point, double* xy)
{
  ebiFunctionBegin;
  FacePointBlended* face_point_blended = (FacePointBlended*)face_point;  
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  int F = face_point_blended->F();
  double* cd = face_point_blended->cd();
  bool found = false;
  int tV, tf;
  double st[6];
  for(int v=0; v<4; v++) {
	 iC( bdsurf.Fcd2Vfcd(BdSurf::EVAL_VALUE, F, cd, v, tV, tf, st) );
	 if(tV==_V) {
		found = true;
		break;
	 }
  }
  if (!found) { cerr << "face_point_xy error" << endl; exit(-1); }
 
  ebiAssert(st[0]<bdsurf.EVAL_UB() && st[1]<bdsurf.EVAL_UB());

  iC( bdsurf.Vfcd2Vxy(BdSurf::EVAL_VALUE, tV, tf, st, xy) ); 
  
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::is_xy_valid"
int BlendedPatch::is_xy_valid(double* xy, bool& is_valid)
{
  ebiFunctionBegin;
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  int f;
  double cd[2];		
  
  iC( bdsurf.Vxy2Vfcd(BdSurf::EVAL_VALUE, _V, xy, f, cd) ); 
  
  is_valid = (cd[0]<bdsurf.EVAL_UB() && cd[1]<bdsurf.EVAL_UB());
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::is_xy_dominant"
int BlendedPatch::is_xy_dominant(double* xy, bool& dominant)
{
  ebiFunctionBegin;
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  int f;
  double cd[2];		
  iC( bdsurf.Vxy2Vfcd(BdSurf::EVAL_VALUE, _V, xy, f, cd) );
  dominant = (cd[0]<=0.5 && cd[1]<=0.5);
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::xy_to_face_point"
int BlendedPatch::xy_to_face_point(double* xy, FacePointOverlapping* face_point)
{
  ebiFunctionBegin;
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  int f;
  double cd[2];		
  iC( bdsurf.Vxy2Vfcd(BdSurf::EVAL_VALUE, _V, xy, f, cd) ); 

  ebiAssert(cd[0]<bdsurf.EVAL_UB() && cd[1]<bdsurf.EVAL_UB()); //assure it is good
  int F;
  double st[2];		  
  iC( bdsurf.Vfcd2Fcd(BdSurf::EVAL_VALUE, _V, f, cd, F, st) );

  FacePointBlended* tmp = (FacePointBlended*)face_point; 
  *tmp = FacePointBlended(F, st);
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::xy_to_patch_coords"
int BlendedPatch::xy_to_patch_coords(double* xy, int flag, double* ret)
{
  ebiFunctionBegin;
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();

  // TODO exterior problems on multiply connected domains is possibly broken
  iC( bdsurf.eval(flag, _V, xy, (Point3*)ret) ); 
  
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::estimate_jacobian"
int BlendedPatch::estimate_jacobian(double* jac)
{
  ebiFunctionBegin;
  double xy[2];  xy[0]=0;  xy[1]=0;
  Point3 ret[3];
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();

  iC( bdsurf.eval(BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, _V, xy, ret) ); 

  jac[0] = max(ret[1].l2(), ret[2].l2());
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BlendedPatch::xy_to_patch_value"
int BlendedPatch::xy_to_patch_value(double* xy, int flag, double* ret)
{
  ebiFunctionBegin;
  BdSurf& bdsurf = ((PatchSurfBlended*)bdry())->bdsurf();
  int f;  double cd[2];		
 
  iC( bdsurf.Vxy2Vfcd(BdSurf::EVAL_VALUE, _V, xy, f, cd) );
  ebiAssert(flag == BdSurf::EVAL_VALUE);
  iC( bdsurf.blendFuncEval(flag, cd, 1.0-bdsurf.EVAL_UB(), bdsurf.EVAL_UB(), ret) );
  

  ebiFunctionReturn(0);
}

END_EBI_NAMESPACE
