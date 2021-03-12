#include "stkfmm_interface.hpp"


class hedgehog::STKFMM::STKFMMImpl{
};
hedgehog::STKFMM::~STKFMM(){
}

void hedgehog::STKFMM::collect_fmm_data(Vec source_positions, Vec source_normals,
		Vec target_positions, Kernel3d kernel){}
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
