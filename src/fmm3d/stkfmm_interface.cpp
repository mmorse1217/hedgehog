#include "stkfmm_interface.hpp"
#include "common/utils.hpp"
#include "common/kernel3d.hpp"
#ifdef HAS_STKFMM
#include <StokesLayerKernel.hpp>
#include <STKFMM.hpp>
#define COUT(str) (std::cout<<str<<std::endl)
#define CERR(str,action) (std::cerr<<"[ERROR]["<< __FUNCTION__ <<"] "<<str<<std::endl, action)
#define ASSERT(expr, msg) ( (expr) ? assert(true) : CERR(msg,abort()))
#define COUTDEBUG(str) (std::cout<<"[DEBUG] "<<str<<std::endl)
const stkfmm::KERNEL get_stkfmm_kernel(int _kernel_type) {
  //const pvfmm::Kernel<double>* kernel;
  stkfmm::KERNEL kernel;
  switch (_kernel_type) {
  case hedgehog::KNL_STK_S_U:
    //kernel = &pvfmm::StokesLayerKernel<double>::Vel();
    kernel = stkfmm::KERNEL::Stokes;
    break;
  case hedgehog::KNL_STK_S_P:
    //kernel = &pvfmm::StokesLayerKernel<double>::PVel();
    kernel = stkfmm::KERNEL::Stokes;
    break;
  case hedgehog::KNL_STK_D_U:
    //kernel = &ker_stokes_dl;
    //kernel = &pvfmm::StokesLayerKernel<double>::Vel();
    kernel = stkfmm::KERNEL::Stokes;
    break;
  case hedgehog::KNL_STK_D_P:
    // kernel = &pvfmm::StokesKernel<double>::pressure();
    //kernel = &ker_stokes_pressure_dl;
    //kernel = &pvfmm::StokesLayerKernel<double>::PVel();
    kernel = stkfmm::KERNEL::PVel;
    break;
  case hedgehog::KNL_LAP_D_U:
    // kernel = &pvfmm::LaplaceKernel<double>::potential();
    //kernel = &ker_laplace_dl;
    //kernel = &pvfmm::LaplaceKernel<double>::potential();
    kernel = stkfmm::KERNEL::LapPGrad;
    break;
  case hedgehog::KNL_LAP_S_U:
    // kernel = &pvfmm::LaplaceKernel<double>::potential();
    //kernel = &ker_laplace_sl;
    //kernel = &pvfmm::LaplaceKernel<double>::potential();
    kernel = stkfmm::KERNEL::LapPGrad;
    break;
  default:
    ASSERT(false, "KernelNotImplementedError");
  }
  return kernel;
}

class hedgehog::STKFMM::STKFMMImpl{
  //typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  //typedef pvfmm::FMM_Tree<Mat_t> Tree_t;
  int mult_order;
  int max_pts;
  int max_depth;
  int sl;
  int dl;
  int periodic;
  int source_dof;
  int target_dof;
  int dense_eval;
  stkfmm::PAXIS periodicity_type;
  MPI_Comm comm;

  //pvfmm::BoundaryType boundary;
  //const pvfmm::Kernel<real_t>* ker;

  std::unique_ptr<stkfmm::STKFMM> fmm;
  public:

  STKFMMImpl(Kernel3d legacy_kernel)  {

    max_depth = Options::get_int_from_petsc_opts("-bis3d_maxlevel");
    mult_order = Options::get_int_from_petsc_opts("-bis3d_np");
    max_pts = Options::get_int_from_petsc_opts("-bis3d_ptsmax");
    periodic = false;
    string periodicity_type_str = Options::get_string_from_petsc_opts("-stkfmm_periodicity_type");
    dense_eval = Options::get_int_from_petsc_opts("-direct_eval");
    
    if (periodicity_type_str == "none") {
      // free space
      periodicity_type = stkfmm::PAXIS::NONE;
    } else if (periodicity_type_str == "x") {
      // x-axis periodicity
      periodicity_type = stkfmm::PAXIS::PX;
    } else if (periodicity_type_str == "xy") {
      // x- and y-axis periodicity
      periodicity_type = stkfmm::PAXIS::PXY;
    } else if (periodicity_type_str == "xyz") {
      // x-, y- and z-axis periodicity
      periodicity_type = stkfmm::PAXIS::PXYZ;
    } else {
      ASSERT(false, "NoSuchPeriodicityType");
    }

    source_dof = legacy_kernel.srcDOF();
    target_dof = legacy_kernel.trgDOF();


    // TODO factor this out when user specified communicators are added
    comm = MPI_COMM_WORLD;
    using namespace stkfmm;
    Equation_type equation_type;
    Kernel_type kernel_type;
    Kernel_variable kernel_variable;
    Kernel3d::parse_kernel_enum((int) legacy_kernel.kernelType(), equation_type, 
          kernel_type,  kernel_variable);

    stkfmm::KERNEL kernel = get_stkfmm_kernel(kernel_type);
    const int maxPoints = max_pts;
    fmm = std::make_unique<stkfmm::Stk3DFMM>(mult_order, max_pts, periodicity_type, static_cast<unsigned int>(kernel));
    fmm->showActiveKernels();

  };

  ~STKFMMImpl(){
  };
};

#else
class hedgehog::STKFMM::STKFMMImpl{};
#endif
hedgehog::STKFMM::~STKFMM(){
}
hedgehog::STKFMM::STKFMM(Vec source_positions, Vec source_normals,
        Vec target_positions, Kernel3d kernel){
    /*_store_problem_data(source_positions, source_normals, target_positions, kernel);
  
    _setup_context();
  
    // Scale sources and targets to [0,1] x [0,1] x [0,1]
    //Note: resulting scaled positions are restored in _src_positions and
    //_trg_positions in order to provide access to points in wrapper lib
    _scale_factor = scale_sources_and_targets(_src_positions, _trg_positions);
    
    _kernel_opts.scale_factor = _scale_factor;*/
}

int hedgehog::STKFMM::setFromOptions(){}
int hedgehog::STKFMM::setup(){}

/* Required evaluation function of the FMM library.
 * Call from inside GMRES matvec.
 */
// Note: function signature may need to change.
void hedgehog::STKFMM::evaluate( const Vec& srcDen, Vec& trgVal){}
void hedgehog::STKFMM::evaluate(const DblNumVec& srcDen, DblNumVec& trgVal){}
void hedgehog::STKFMM::interaction_matrix(DblNumMat source_positions, 
		DblNumMat source_normals, DblNumMat target_positions, 
		DblNumMat& interaction_matrix){}

void hedgehog::STKFMM::set_kernel(Kernel3d& kernel){}

// TODO MJM verify that thse &'s are needed
void hedgehog::STKFMM::set_src_positions(DblNumMat*& src_positions){}
void hedgehog::STKFMM::set_src_normals(DblNumMat*& src_normals){}
void hedgehog::STKFMM::set_trg_positions(DblNumMat*& trg_positions){}
void hedgehog::STKFMM::set_rebuild_tree(bool rebuild){}

	/*
	   if (config.direct) {
	   auto srcSLCoord = point.srcLocalSL; // a copy
		auto srcSLValue = value.srcLocalSL; // a copy
		auto srcDLCoord = point.srcLocalDL; // a copy
		auto srcDLValue = value.srcLocalDL; // a copy
		auto trgLocalCoord = point.trgLocal;
		PointDistribution::collectPtsAll(srcSLCoord); // src from all ranks
		PointDistribution::collectPtsAll(srcSLValue); // src from all ranks
		PointDistribution::collectPtsAll(srcDLCoord); // src from all ranks
		PointDistribution::collectPtsAll(srcDLValue); // src from all ranks

		const int nSL = srcSLCoord.size() / 3;
		const int nDL = srcDLCoord.size() / 3;
		const int nTrg = trgLocalCoord.size() / 3;
		trgLocalValue.clear();
		trgLocalValue.resize(kdimTrg * nTrg, 0);

		timer.tick();
		if (nSL)
			fmmPtr->evaluateKernel(kernel, 0, PPKERNEL::SLS2T,                //
								   nSL, srcSLCoord.data(), srcSLValue.data(), //
								   nTrg, trgLocalCoord.data(), trgLocalValue.data());
		if (nDL)
			fmmPtr->evaluateKernel(kernel, 0, PPKERNEL::DLS2T,                //
								   nDL, srcDLCoord.data(), srcDLValue.data(), //
								   nTrg, trgLocalCoord.data(), trgLocalValue.data());
		timer.tock("evaluateKernel");

		const auto &time = timer.getTime();
		treeTime = 0;
		runTime = time[0];
	} else {
		double origin[3] = {config.origin[0], config.origin[1], config.origin[2]};
		const double box = config.box;
		const int nSL = point.srcLocalSL.size() / 3;
		const int nDL = point.srcLocalDL.size() / 3;
		const int nTrg = point.trgLocal.size() / 3;
		trgLocalValue.clear();
		trgLocalValue.resize(kdimTrg * nTrg, 0);

		fmmPtr->clearFMM(kernel);
		fmmPtr->setBox(origin, box);
		fmmPtr->setPoints(nSL, point.srcLocalSL.data(), nTrg, point.trgLocal.data(), nDL, point.srcLocalDL.data());

		timer.tick();
		fmmPtr->setupTree(kernel);
		timer.tock("setupTree");

		timer.tick();
		fmmPtr->evaluateFMM(kernel, nSL, value.srcLocalSL.data(), //
							nTrg, trgLocalValue.data(),           //
							nDL, value.srcLocalDL.data());
		timer.tock("evaluateFMM");
		const auto &time = timer.getTime();
		treeTime = time[0];
		runTime = time[1];
 */
