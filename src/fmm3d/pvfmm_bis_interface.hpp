#ifndef __PVFMM_INTERFACE_HPP__
#define __PVFMM_INTERFACE_HPP__

#include "fmm_interface.hpp"
#include <pvfmm_interface.h>
#include <vector>
#include "common/kernel3d.hpp"
#include "fmm3d/tree.hpp"
#include "common/ebi_petsc.hpp"

BEGIN_EBI_NAMESPACE

typedef std::vector<double> vec;

class PvFMM: public FMM {
        
    protected:
        vec _src_positions;
        vec _src_normals;
        DblNumMat _src_normals_temp;
        DblNumMat* _markgrid_targets;
        vec _trg_positions;
        int _num_points;
        Kernel _kernel_type;

        Equation_type _equation_type;
        Kernel_type _kernel_layer;
        Kernel_variable _kernel_variable;

        void *_context; 
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
        ~PvFMM();
        
        //PvFMMTree* tree(){ return _tree; }
        int setup();
        int setFromOptions();

        int evaluate(const Vec& srcDen, Vec& trgVal);
        int evaluate(const DblNumVec& srcDen, DblNumVec& trgVal);
        
        // O(N^2) kernel summaiton 
        int evaluate_direct(const DblNumVec& srcDen, DblNumVec& trgVal);
        int evaluate_direct(const Vec& srcDen, Vec& trgVal);
        int interaction_matrix(DblNumMat source_positions, DblNumMat source_normals, 
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

};

END_EBI_NAMESPACE
#endif
