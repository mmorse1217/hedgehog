
#include "../catch.hpp"
#include <bdry3d/patch_surf_nanospline.hpp>
#include <memory>
TEST_CASE("Test nanospline interface", "[nanospline]"){
    unique_ptr<Ebi::PatchSurfNanospline> surface(new Ebi::PatchSurfNanospline("", ""));
}
