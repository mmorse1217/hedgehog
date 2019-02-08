#include "../catch.hpp"
#include "common/vecmatop.hpp"
#include <sampling.hpp>
#include "common/vtk_writer.hpp"

using namespace Ebi;
using namespace Sampling;
using Sampling::equispaced;
using Sampling::sample_2d;
using Sampling::base_domain;

TEST_CASE("Test patch refinement", "[refinment][fft][patch][no-fmm]"){
    //cout.precision(16);
    SECTION("Scalar test periodic function"){
        /*
         * fftrf2d(int dof, double* data, int* mm, int* rr, double* res)
         * dof - context: source degree of freedom
         *       useage:  stride size?
         * data- context: scaled density values at regular sample points on patch
         *       useage:  2d periodic function to perform FFT on
         * mm  - context: number of samples in 1d along square domain to FFT
         *       useage:  original num samples in 1d
         * rr  - context: number of refined samples in 1d along domain to FFT
         *       useage:  final num samples in 1d
         * res - context: refined resulting denisty 
         *       useage:  output refined data
         */
        cout << "fft test" << endl;

        int num_samples = 10;
        DblNumMat uv(2, num_samples*num_samples, true,
                sample_2d<equispaced>(num_samples, base_domain).data());

        int refinement_factor = 10;
        int ref_num_samples = refinement_factor*num_samples;
        DblNumMat refined_uv(2, ref_num_samples*ref_num_samples, true,
                sample_2d<equispaced>(ref_num_samples, base_domain).data());

        DblNumMat f(1, num_samples*num_samples);
        for (int i = 0; i < uv.n(); i++) {
            double u = uv(0,i);
            double v = uv(1,i);

            f(0,i) = cos(2*M_PI*u)*cos(2*M_PI*v);
            //f(0,i) = 1.;

        }

        int refined_num_samples = num_samples*refinement_factor;
        DblNumMat refined_f(1, ref_num_samples*ref_num_samples); clear(refined_f);

        //refinement using fft: do fourier transform, increase the
        //number of samples by refinement factor padding with zeros, transform back 
        int mn[2];
        mn[0] = num_samples;
        mn[1] = num_samples;

        int rs[2];
        rs[0] = refinement_factor;
        rs[1] = refinement_factor;
        fftrf2d(1, f.data(), mn, rs, refined_f.data());
        //DblNumMat error(refined_f.m(), refined_f.n());
        for (int i = 0; i < refined_uv.n(); i++) {
            double computed_f_value = refined_f(0,i); 
            double u = refined_uv(0,i);
            double v = refined_uv(1,i);
            double true_f_value = cos(2*M_PI*u)*cos(2*M_PI*v);
            //double true_f_value = 1.;
            cout << computed_f_value <<", " << true_f_value << endl;
            double rel_error = fabs(computed_f_value - true_f_value)/fabs(true_f_value);
            //error(0,i) = rel_error;
            CHECK( rel_error<= 1e-4);
        }
        
        
        
        
        DblNumMat out_ref(3,ref_num_samples*ref_num_samples);
        DblNumMat out(3,num_samples*num_samples);
NumVec<OnSurfacePoint> temp(refined_uv.n());
        for (int i = 0; i < refined_uv.n(); i++) {
            temp(i).parent_patch = -1;
            out_ref(0,i) = refined_uv(0,i); 
            out_ref(1,i) = refined_uv(1,i); 
            out_ref(2,i) = refined_f(0,i); 
        }
        for (int i = 0; i < uv.n(); i++) {
            out(0,i) = uv(0,i); 
            out(1,i) = uv(1,i); 
            out(2,i) = f(0,i); 
        }
        write_qbkix_points_to_vtk(out, 
               temp, 
                 0) ;
        write_qbkix_points_to_vtk(out_ref, 
               temp, 
                 0, "data/ref_") ;
    }
}
