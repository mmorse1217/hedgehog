#include "common/kernel3d.hpp"
#include "common/nummat.hpp"
#include <pvfmm_interface.h>
#include "fmm3d/pvfmm_bis_interface.hpp"
using namespace std;
using namespace hedgehog;
int main(int argc, char** argv){

    Kernel3d solver_kernel(121, vector<double>(2,1.));
    cout << "hello world" << endl;

}
