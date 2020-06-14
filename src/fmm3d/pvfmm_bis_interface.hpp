#ifndef __PVFMM_INTERFACE_HPP__
#define __PVFMM_INTERFACE_HPP__

#include "fmm_interface.hpp"
#include <memory>
#include <vector>
#include "common/kernel3d.hpp"
#include "fmm3d/tree.hpp"
#include "common/ebi_petsc.hpp"

BEGIN_EBI_NAMESPACE

typedef std::vector<double> vec;

        struct KernelOptions{
            bool initialized;
            double scale_factor;
            double navier_mu;
            double navier_nu;
            double helmholtz_frequency;
        };
class PvFMM: public FMM {
    private:
    void _interleave_density(DblNumVec srcDen);
        void _interface_evaluate(vec source_density, vec& potv); 
        void _setup_tree(vec source_density); 
        void _store_problem_data(Vec source_positions, Vec source_normals,
                Vec target_positions, Kernel3d kernel);
        void _setup_context();
        typedef  double real_t;
        void _copy_potential(int ntrg, int target_dof, std::vector<real_t> potv,
                            real_t *pot, bool rescale=false); 
    protected:
        vec _src_positions;
        vec _src_normals;
        DblNumMat _src_normals_temp;
        DblNumMat* _markgrid_targets;
        vec _trg_positions;
        int _kernel_type;

        Equation_type _equation_type;
        Kernel_type _kernel_layer;
        Kernel_variable _kernel_variable;

        void *_context; 
        class PVFMMImpl;
        unique_ptr<PVFMMImpl> _impl;
        double _scale_factor;
        double _frequency;
        int _direct_eval;
        int _source_dof;
        int _target_dof;
        int _num_sources;
        int _num_targets;
        bool _rebuild_tree;
        KernelOptions _kernel_opts;
        //PvFMMTree* _tree;

    public:

        PvFMM();
        /*PvFMM(Vec source_positions,
                Vec target_positions, Kernel3d kernel);*/
        PvFMM(Vec source_positions, Vec source_normals,
                Vec target_positions, Kernel3d kernel);
        ~PvFMM();
        
        //PvFMMTree* tree(){ return _tree; }
        int setup();
        int setFromOptions();

        void evaluate(const Vec& srcDen, Vec& trgVal);
        void evaluate(const DblNumVec& srcDen, DblNumVec& trgVal);
        
        // O(N^2) kernel summaiton 
        void evaluate_direct(const DblNumVec& srcDen, DblNumVec& trgVal);
        void evaluate_direct(const Vec& srcDen, Vec& trgVal);
        void interaction_matrix(DblNumMat source_positions, DblNumMat source_normals, 
                DblNumMat target_position, DblNumMat& interaction_matrix);
        double scale_sources_and_targets(vec sources, vec targets);
        int compute_pvfmm_to_kifmm_scaling(const DblNumVec& srcDen, DblNumVec& trgVal);
        //void initialize_tree(DblNumVec srcDen);
        void clear_tree();
        void transpose_data();
        void set_src_positions(DblNumMat*& src_positions); 
        void set_src_normals(DblNumMat*& src_normals);
        void set_trg_positions(DblNumMat*& trg_positions);
        int& set_num_points(int64_t& num_points);
        void set_kernel(Kernel3d& kernel);

        void set_rebuild_tree(bool rebuild){
            _rebuild_tree = rebuild;
        }
        DblNumMat get_markgrid_targets(){
            return *_markgrid_targets;
        }

void laplace_sldl_pvfmm(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *pot,
                    const int &rebuild_tree, void **context);
void pvfmm_direct_summation(
        vec source_positions, vec source_densities, vec target_positions,
        vec& target_potentials);
/*void pvfmm_direct_summation(size_t num_source_points, real_t *source_positions,
                            real_t *source_densities, size_t ntrg, real_t *trg,
                            real_t *pot);*/
void pvfmm_interaction_matrix(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *pot,
                    const int &rebuild_tree, void **context);
};

END_EBI_NAMESPACE
#endif
