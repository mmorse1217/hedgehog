
#include "../catch.hpp"
#include <bdry3d/patch_surf_nanospline.hpp>
#include <memory>
#include <sampling.hpp>

using Sampling::sample_2d;
using Sampling::base_domain;
using Sampling::equispaced;

TEST_CASE("Test nanospline interface for flat patch", "[nanospline]"){
    unique_ptr<Ebi::PatchSurfNanospline> surface(new Ebi::PatchSurfNanospline("", ""));
    auto patch = surface->patch(0);
    const int num_samples = 20;
    Ebi::DblNumMat parametric_samples = DblNumMat(2, num_samples*num_samples, true, 
                sample_2d<equispaced>(num_samples, base_domain).data());
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
