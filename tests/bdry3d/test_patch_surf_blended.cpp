
#include "../catch.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "common/utils.hpp"
#include "common/stats.hpp"
#include "../utils/evaluation_utils.hpp"
#include "../utils/regression_utils.hpp"
#include "common/interpolate.hpp"
#include <bdsurf.hpp>
#include <sampling.hpp>
#include <utils.hpp>
using namespace hedgehog;
using Sampling::sample_2d;
using Sampling::equispaced;
using Sampling::base_domain;

Index2 idx(Point2 xy, double _init, double _step, int _num_samples){
    double xya[] = {xy.x(), xy.y()};
    double hashed_xy_arr[] = {(xya[0] - _init)/_step,
                            (xya[1] - _init)/_step};
    int index_xy_arr[] = {int(floor(hashed_xy_arr[0])), int(floor(hashed_xy_arr[1]))};
    Index2 array_based(
            min(max(index_xy_arr[0],1),_num_samples-3),
            min(max(index_xy_arr[1],1),_num_samples-3)
            );
    return array_based;
}
TEST_CASE("Test blendsurf interpolation code", "[blendsurf][interpolate]"){
    SECTION("interpolate a vanilla function"){
        int num_samples = 10;
        double step = 1./double(num_samples-1);
        int num_eval= 11;

         /*
        NumMatrix interpolation_nodes=//(2,num_samples*num_samples,true, 
            sample_2d<equispaced>(num_samples, base_domain);
          */
        NumMatrix interpolation_nodes(2,num_samples*num_samples);
        for (int i = 0; i < num_samples; i++) {
            for (int j = 0; j < num_samples; j++) {
                int index = i*num_samples+j;
                interpolation_nodes(0,index) = double(i)*step;
                interpolation_nodes(1,index) = double(j)*step;
            }
        } 
        NumMatrix evaluation_points = //(2,num_eval*num_eval,true, 
            sample_2d<equispaced>(num_eval, 
                    Rectangle(Interval(.3, .7),Interval(.3, .7)));

        NumMatrix function_values(3,num_samples*num_samples);
        NumMatrix interpolated_values(3,num_eval*num_eval);
        NumMatrix true_values(3,num_eval*num_eval);

        // evaluate function
        for (int i = 0; i < num_samples; i++) {
            for (int j = 0; j < num_samples; j++) {
                int index = i*num_samples+j;
                Point2 p(interpolation_nodes.clmdata(index));
                for (int d = 0; d < 3; d++) {
                    function_values(d,index) = p.x()*p.y();//cos(2*M_PI*p.y())*cos(2*M_PI*p.x());

                }

            }
        }

        //compute true values at nodes
        for (int i = 0; i < num_eval; i++) {
            for (int j = 0; j < num_eval; j++) {
                int index = i*num_eval+j;
                Point2 p(evaluation_points.clmdata(index));
                for (int d = 0; d < 3; d++) {
                    true_values(d,index) = p.x()*p.y();//cos(2*M_PI*p.y())*cos(2*M_PI*p.x());

                }
            }
        }

        NumVector temp_mat(4);
        for (int i = 0; i < 4; i++) {
            temp_mat(i) = interpolation_nodes(1,i); 
        }
        NumVector bary_weights =  bary_weights_1d(temp_mat);
        for (int i = 0; i < num_eval; i++) {
            for (int j = 0; j < num_eval; j++) {
                int index = i*num_eval+j;
                Point2 p(evaluation_points.clmdata(index));
                Point2 scaled_p(p.x()/step, p.y()/step);
                Point3 f(interpolated_values.clmdata(index));
                Index2 index_xy = idx(p, 0., step, num_samples);
                bary2d_2(index_xy, 3, interpolation_nodes, function_values, 
                         bary_weights.data(), p, f);

               /*cout << index_xy << endl;
               cout << p << endl;
               cout << f << endl;*/
                for (int d = 0; d < 3; d++) {
                interpolated_values(d,index) = f(d);
                }
            }
        }
                /*cout << temp_mat<< endl;
                cout << bary_weights<< endl;
                cout << interpolation_nodes << endl;
                cout << function_values<< endl;*/
        for (int i = 0; i < num_eval; i++) {
            for (int j = 0; j < num_eval; j++) {
                int index = i*num_eval+j;
                Point3 p_f(interpolated_values.clmdata(index));
                Point3 f(true_values.clmdata(index));
                //cout << p_f << ", " << f << p_f.l2()/f.l2() << endl;
                CHECK((f-p_f).length() <= 1e-3);
            }
        }

    }
}
