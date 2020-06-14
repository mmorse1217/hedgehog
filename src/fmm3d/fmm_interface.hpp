#ifndef __FMM_INTERFACE_HPP__
#define __FMM_INTERFACE_HPP__
#include "common/ebi.hpp"
#include "common/numtns.hpp"
#include "common/vecmatop.hpp"
//#include "common/syms.hpp"
#include "common/numvec.hpp"
#include <vec3t.hpp>
#include "common/kernel3d.hpp"
#include <stack>

BEGIN_EBI_NAMESPACE
/*
 * A template class for a standalone kernel-independent fast multipole
 * method implementation. Implementing the virtual functions here will run the
 * boundary integral solver with the desired FMM library.
 */


class FMM {
    protected:
        DblNumMat* _srcPos;
        DblNumMat* _srcNor;
        DblNumMat* _trgPos;
        Kernel3d _kernel;
        // FMM properties. Loaded from options files
        int64_t _multipole_order;
        int64_t _max_level;
        int64_t _max_points_per_box;


    public:
        FMM();
        //virtual ~FMM(){;}
        virtual ~FMM();
        
        void collect_fmm_data(Vec source_positions, Vec source_normals,
                Vec target_positions, Kernel3d kernel);
        virtual int setFromOptions() = 0;
        virtual int setup() = 0;

        /* Required evaluation function of the FMM library.
         * Call from inside GMRES matvec.
         */
        // Note: function signature may need to change.
        virtual void evaluate( const Vec& srcDen, Vec& trgVal) = 0;
        virtual void evaluate(const DblNumVec& srcDen, DblNumVec& trgVal) = 0;
        virtual void interaction_matrix(DblNumMat source_positions, 
                DblNumMat source_normals, DblNumMat target_positions, 
                DblNumMat& interaction_matrix) = 0;
        
        virtual void set_kernel(Kernel3d& kernel) = 0;

        // TODO MJM verify that thse &'s are needed
        virtual void set_src_positions(DblNumMat*& src_positions) = 0;
        virtual void set_src_normals(DblNumMat*& src_normals) = 0;
        virtual void set_trg_positions(DblNumMat*& trg_positions) = 0;
        virtual void set_rebuild_tree(bool rebuild) = 0;
        DblNumMat*& srcPos() { return _srcPos; } 
        DblNumMat*& srcNor() { return _srcNor; } 
        DblNumMat*& trgPos() { return _trgPos; } 
        

        int64_t multipole_order(){
            return _multipole_order;
        }

        int64_t max_level(){
            return _max_level;
        }

        int64_t max_points_per_box(){
            return _max_points_per_box;
        }
        /* Minimal required data required to pass along to a general FMM library
         * A particular implementation may require more details, but additional
         * inputs can be passed as system options and parsed in setup().
         */
        // FMM set up 
        void initialize_fmm(Vec source_positions, Vec source_normals, 
                Vec target_positions, Kernel3d kernel, bool setup=true);



};



END_EBI_NAMESPACE
#endif
