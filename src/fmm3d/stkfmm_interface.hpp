#pragma once

#include "fmm_interface.hpp"
#include "common/kernel3d.hpp"
#include "common/ebi_petsc.hpp"
namespace hedgehog {
	class STKFMM : public FMM {
		private:
			class STKFMMImpl;
			unique_ptr<STKFMMImpl> _impl;
		public:
			STKFMM();
			~STKFMM();

			void collect_fmm_data(Vec source_positions, Vec source_normals,
					Vec target_positions, Kernel3d kernel);
			int setFromOptions();
			int setup();

			/* Required evaluation function of the FMM library.
			 * Call from inside GMRES matvec.
			 */
			// Note: function signature may need to change.
			void evaluate( const Vec& srcDen, Vec& trgVal);
			void evaluate(const DblNumVec& srcDen, DblNumVec& trgVal);
			void interaction_matrix(DblNumMat source_positions, 
					DblNumMat source_normals, DblNumMat target_positions, 
					DblNumMat& interaction_matrix);

			void set_kernel(Kernel3d& kernel);

			// TODO MJM verify that thse &'s are needed
			void set_src_positions(DblNumMat*& src_positions);
			void set_src_normals(DblNumMat*& src_normals);
			void set_trg_positions(DblNumMat*& trg_positions);
			void set_rebuild_tree(bool rebuild);
	};
}
