
#include "../catch.hpp"
#include <bdry3d/patch_surf_nanospline.hpp>
#include <memory>
#include <random>
#include <sampling.hpp>

using Sampling::sample_2d;
using Sampling::base_domain;
using Sampling::equispaced;

TEST_CASE("Test nanospline interface for flat patch", "[nanospline]"){
    unique_ptr<hedgehog::PatchSurfNanospline> surface(new hedgehog::PatchSurfNanospline("", ""));
    auto patch = surface->patch(0);
    const int num_samples = 20;
    hedgehog::DblNumMat parametric_samples = DblNumMat(2, num_samples*num_samples, true, 
                sample_2d<equispaced>(num_samples, base_domain).data());
    CHECK(patch->characteristic_length() == Approx(2.));
    for (int i=0; i<num_samples; i++) {
      for (int j = 0; j < num_samples; j++) {
          int index = i*num_samples +j;
          Point2 uv(parametric_samples.clmdata(index));
          Point3 values[6];
          patch->xy_to_patch_coords(uv.array(), patch->EVAL_VL|patch->EVAL_FD|patch->EVAL_SD, (double *)values);
          // TODO add test for positions (values[0]) one parametrization is correct
          
          CHECK(values[1].length() ==Approx(2.));
          CHECK(values[2].length() ==Approx(2.));
          
          CHECK(values[1].x() ==Approx(2.));
          CHECK(values[2].y() ==Approx(2.));
          for (int i=0; i<3; i++) {
              CHECK(values[3+i].length() <=1e-12);
          }
      }
    }
}
TEST_CASE("Test nanospline deform for flat patch", "[nanospline][deform]"){
    unique_ptr<hedgehog::PatchSurfNanospline> surface(new hedgehog::PatchSurfNanospline("", ""));
    auto patch = dynamic_cast<hedgehog::NanosplinePatch*>(surface->patch(0));
    std::uniform_real_distribution<double> unif(0., 1.);
    std::default_random_engine re(0);

    //int num_constraints = 5000;
    int num_constraints_x = 10;
    int num_constraints_y = 10;
    int num_constraints = num_constraints_x*num_constraints_y;
    DblNumMat changes_in_position(3, num_constraints);
    DblNumMat parameter_values(2, num_constraints);
    for (int i =0; i < num_constraints_x; i++) {
    for (int j =0; j < num_constraints_y; j++) {
        int index = i*num_constraints_x + j;
        changes_in_position(0,index) = 0;
        changes_in_position(1,index) = 0;
        changes_in_position(2,index) = .1;
        parameter_values(0,index) =  double(i)/double(num_constraints_x);
        parameter_values(1,index) = double(j)/double(num_constraints_y);
    }
    }
    Point3 p;
    Point2 uv(0.);
    patch->xy_to_patch_coords(uv.array(), patch->EVAL_VL, p.array());
    cout << "0,0: " << p << endl;
    uv.x() = 1.;
    patch->xy_to_patch_coords(uv.array(), patch->EVAL_VL, p.array());
    cout << "1,0: " << p << endl;
    uv.x() = 0.;
    uv.y() = 1.;
    patch->xy_to_patch_coords(uv.array(), patch->EVAL_VL, p.array());
    cout << "0,1: " << p << endl;
    uv.x() = 1.;
    uv.y() = 1.;
    patch->xy_to_patch_coords(uv.array(), patch->EVAL_VL, p.array());
    cout << "1,1: " << p << endl;
    patch->deform_periodic(parameter_values, changes_in_position);


}
