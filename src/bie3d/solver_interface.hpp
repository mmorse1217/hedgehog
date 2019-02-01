#ifndef _BIS3DOV_HPP_
#define _BIS3DOV_HPP_

#include "common/kernel3d.hpp"
#include "bdry3d/patch_samples.hpp"
#include "common/numtns.hpp"
#include "fmm3d/pvfmm_bis_interface.hpp"
#include "common/nummat.hpp"

BEGIN_EBI_NAMESPACE

class SolutionDensity: public EbiObject {
    private:
        Vec _solution;
        Vec _density;
        Vec _pole_coeffs;

        int _local_num_points;
        int _global_num_points;
        int _number_of_poles;
        int _domain_boundedness;
        Kernel3d _kernel;
        bool _read_only;
        
        // Pointer to local array of _solution
        double* _vec_ptr;

    public:
    SolutionDensity(): EbiObject("", "") {;}
    
    SolutionDensity(Vec &v, int local_num_points, int global_num_points,
            int number_of_poles, int domain_bounded, Kernel3d kernel,
            bool read_only=false):
        EbiObject("","") {

        // initialize member variables
        _solution = v;
        _local_num_points = local_num_points;
        _global_num_points = global_num_points;
        _number_of_poles = number_of_poles;
        _kernel = kernel;
        _domain_boundedness = domain_bounded;
        _read_only = read_only;

        set_petsc_vec(_solution);
    }
    SolutionDensity(const SolutionDensity& s): EbiObject("", ""){
        _solution = s._solution;
        _local_num_points = s._local_num_points;
        _global_num_points = s._global_num_points;
        _number_of_poles = s._number_of_poles;
        _kernel = s._kernel;
        _domain_boundedness = s._domain_boundedness;
    }

    ~SolutionDensity(){
        
        // Destroy new sub-Vec's

        // Restore the subarrays _density and _pole_coeffs to the original 
        // input Vec _solution.
        VecRestoreArray(_solution, &_vec_ptr);
        // MJM MEMORY LEAK add these destroy's in (test them with GMRES)
        //VecDestroy(&_density);
        //VecDestroy(&_pole_coeffs);
       
    }

    // yank out density and pole coeffs from input vector:
    // in = [density  pole_coeffs] where density is of size "ga" and 
    // pole_coeffs is of size "gb"
    void set_petsc_vec(Vec v){
        _solution = v;
        // make each separete or return a struct
        int la, lb, lm, ga, gb, gm;
        localSize(la, lb, lm);
        globalSize(ga, gb, gm);

        // get a pointer to the local array of the Petsc Vec
        VecGetArray(_solution, &_vec_ptr);

        // Allocate the actual sub-Vec's for the density and the pole_coefficients
        VecCreateMPIWithArray(this->mpiComm(), 1, la, ga, _vec_ptr, &_density);
        VecCreateMPIWithArray(this->mpiComm(), 1, lb, gb, _vec_ptr+la, &_pole_coeffs);
    }

    void localSize( int&, int&, int&);
    void globalSize(int&, int&, int&);
    
    Vec density(){
        return _density;
    }
    Vec pole_coeffs(){
        return _pole_coeffs;
    }

    Vec solution(){
        return _solution;
    }
};

enum {
  // Domain type
  DOM_BND = 0, // Bounded
  DOM_UNBND = 1, // Unbounded
  DOM_ERR = -1
};

//in nea eval, IN==1, OUT=-1
//in bdry, normalout==1, normalin==-1

//-----------------------------------------------
class SolverInterface: public EbiObject
{
protected:
  //PARAMS
  PatchSurf* _bdry;
  PatchSurf* _refined_surface;
  SolutionDensity* _reference_solution;

  // Indicates which parallel partition the ith patch is a part of
  vector<int> _patch_partition; 

  // Whether the domain is bounded (=0) or unbounded (=1) 
  int _dom;


  vector<double> _eqcoefs;

    
  // Discretization of patch representation
  PatchSamples* _patch_samples;

  // Refeined discretization of patch representation
  PatchSamples* _refined_patch_samples;

  

  // FMM 
  PvFMM* _pvfmm;
  
public:
  SolverInterface(const string& n, const string& p);
  ~SolverInterface();  
  //virtual void setFromOptions();

  virtual int setup() = 0;
  /*
   * Calls Petsc's GMRES to solve the discretized integral equations (6), (7),
   * or (9). The fast numerical quadrature outlined in Section 3 is used to
   * evaluate (D \phi)(x) inside of GMRES for x \in \Gamma.
   *
   * @param Vec         b0      right-hand side of the integral equation of
   *                            interest, denoted \mathbf{b}(x) in the paper.
   * @param Vec         x       the solution density \phi
   *
   */
  virtual int solve(Vec b0, Vec x) = 0;

  /*
   * Evaluates the solution density in the well-seperated region (\Omega_0 in 
   * Section 4). Uses the trapezoidal rule and the kernel independent FMM to
   * evaluate targets of distance >= \sqrt(h) from the boundary in O(N^{3/2})
   * time.
   *
   * @param Vec        tp       target positions
   * @param int        qt       flag indicating what to be evaluated (velocity
   *                                or pressure for stokes)
   * @param Vec        den      solution density found via solve(...)
   * @param Vec        val      PDE solution evaluted at target points in the
   *                            volume
   */
  virtual int fareval(Vec tp, int qt, Vec den, Vec val) = 0;


  /*
   * Evaluates the solution density in the intermediate region (\Omega_1 in 
   * Section 4). Uses the trapezoidal rule and the kernel independent FMM to
   * evaluate targets of h < distance <= \sqrt(h) from the boundary in O(N^{3/2})
   * time.
   *
   * @param Vec        tp       target positions
   * @param Vec        ti       \in {1, -1}, indicates if a target point is
   *                            inside or outside the domain \Omega (?)
   * @param Vec        bp       
   * @param int        qt       flag indicating what to be evaluated (velocity
   *                                or pressure for stokes)
   * @param Vec        den      solution density found via solve(...)
   * @param Vec        val      PDE solution evaluted at target points in the
   *                            volume
   * @param Vec        reg      indicates if a sample point is in the
   *                            intermediate region (reg[i]==2) or the 
   *                            near-region (reg[i] == 1) 
   */
  int neaeval(Vec tp, Vec ti, Vec bp, Vec bb, int qt, Vec den,
      Vec val, Vec reg);

 //virtual int localSize( int&, int&, int&)=0;
 //virtual int globalSize(int&, int&, int&)=0;
  
  
  //accessor methods
  PatchSurf*& bdry() { return _bdry; }
  PatchSurf*& refined_surface() { return _refined_surface; }
  PatchSamples* patch_samples() { return _patch_samples; }
  PatchSamples* refined_patch_samples() { return _refined_patch_samples; }
  PvFMM*& pvfmm() { return _pvfmm; } 

  vector<int>& patch_partition() { return _patch_partition; }
  vector<double>& eqcoefs() { return _eqcoefs; }
 
  void localSize( int& la, int&lb , int& lm ){ 
      _reference_solution->localSize(la, lb, lm);
  }
  void globalSize(int& ga, int& gb, int& gm ){ 
      _reference_solution->globalSize(ga, gb, gm);
  }
  int local_sample_dof(){
      int la; int lb; int lm;
      localSize(la, lb, lm);
      return la;
  }
  int local_pole_dof(){
      int la; int lb; int lm;
      localSize(la, lb, lm);
      return lb;
  }
  int local_total_dof(){
      int la; int lb; int lm;
      localSize(la, lb, lm);
      return lm;
  }
  int global_sample_dof(){
      int la; int lb; int lm;
      globalSize(la, lb, lm);
      return la;
  }
  int global_pole_dof(){
      int la; int lb; int lm;
      globalSize(la, lb, lm);
      return lb;
  }
  int global_total_dof(){
      int la; int lb; int lm;
      globalSize(la, lb, lm);
      return lm;
  } 
  int& dom() { return _dom; } 


};

END_EBI_NAMESPACE

#endif
