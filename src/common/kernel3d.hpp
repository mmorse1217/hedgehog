/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#ifndef _KERNEL3D_HPP_
#define _KERNEL3D_HPP_

#include "nummat.hpp"

BEGIN_EBI_NAMESPACE

using std::vector;

//eqt: 1 2 3 4 5 6
//lyr: s d r p
//qnt: u p ...
enum Equation_type {
    LAPLACE         = 100,
    MOD_HELMHOLTZ   = 200,
    STOKES          = 300,
    UNSTEADY_STOKES = 400,
    NAVIER          = 500,
    OTHER           = 900 
};

enum Kernel_type {
    SINGLE_LAYER    = 10,
    DOUBLE_LAYER    = 20,
    ROTLET          = 30,
    IDENTITY        = 90
};

enum Kernel_variable {
    VAR_U           = 1,
    VAR_P           = 2
};

enum {
  /*! Laplace kernel - Single Layer */
  KNL_LAP_S_U = 111,
  KNL_LAP_E_U = 112,
  /*! Laplace kernel - Double Layer */
  KNL_LAP_D_U = 121,
  /*! Laplace kernel - Identity Tensor */
  KNL_LAP_I   = 191,
  /*! Modified La Single Layer */
  KNL_MODHEL_S_U = 211,
  /*! Modified Lap Double Layer */
  KNL_MODHEL_D_U = 221,
  /*! Stokes kernel - F Velocity */
  ///kernel used by FMM3d algorithm for stokes equation
  KNL_STK_F_U = 301,
  /*! Stokes kernel - Single Layer Velocity */
  KNL_STK_S_U = 311,
  /*! Stokes kernel - Single Layer Pressure */
  KNL_STK_S_P = 312,
  /*! Stokes kernel - Double Layer Velocity */
  KNL_STK_D_U = 321,
  /*! Stokes kernel - Double Layer Pressure */
  KNL_STK_D_P = 322,
  /*! Stokes kernel - R Velocity */
  KNL_STK_R_U = 331,
  /*! Stokes kernel - R Pressure */
  KNL_STK_R_P = 332,
  /*! Stokes kernel - Identity Tensor */
  KNL_STK_I   = 391,
  /*! Stokes kernel - Levi-Civita Tensor */
  KNL_STK_E   = 392,

  /*! Unsteady Stokes */
  KNL_UNSTK_F_U = 401,
  KNL_UNSTK_S_U = 411,
  /*! Unsteady Stokes kernel - Single Layer Pressure */
  KNL_UNSTK_S_P = 412,
  /*! Unsteady Stokes kernel - Double Layer Velocity */
  KNL_UNSTK_D_U = 421,
  /*! Unsteady Stokes kernel - Double Layer Pressure */
  KNL_UNSTK_D_P = 422,

  //navier kernels  //KNL_NAV_F_U = 501, //used for fmm
  /*! Stokes kernel - Levi-Civita Tensor */
  KNL_NAV_S_U = 511, //single displacement
  KNL_NAV_D_U = 521, //double displacement
  KNL_NAV_R_U = 531,
  KNL_NAV_I   = 591, //identity tensor
  KNL_NAV_E   = 592, //levi-civita tensor
  //other kernels
  KNL_SQRTLAP = 901,
  KNL_EXP     = 902,
  //error
  KNL_ERR = -1
};

//#define _mindif 1e-12
//-------------------------------
class Kernel3d
{
protected:
  int _kernelType;
  vector<double> _coefs;  
  double _mindif; //minimal difference
  
  // Enum component decomposition
  
  // Equation type:
  // 100 -- Laplace
  // 200 -- Helmholtz
  // 300 -- Stokes
  // >=400 -- Things explode for now
  Equation_type _equation_type;

  // Kernel type:
  // 10 -- Single Layer
  // 20 -- Double Layer
  // 30 -- Rotlet Kernel
  // 40 -- Grad Single Layer Kernel
  // 90 -- Identity and Levi-Cevita (whatever that is)
  Kernel_type _kernel_type;

  // Kernel Variable:
  // 1 -- analoguous to velocity in Stokes
  // 2 -- analogous to pressure in Stokes
  Kernel_variable _kernel_var;
  
  // Target degree of freedom
  int _tdof;

  // Source degree of freedom
  int _sdof;

  // Total degrees of freedom
  int _pdof;

  // Rotational degrees of freedom
  int _rdof;

  // Degrees of freedom introduced by multiply connected boundaries
  // MJM: this maybe shouldn't be here
  int _gdof;

public:
  Kernel3d(): _kernelType(KNL_ERR) {;}

  Kernel3d(int kernelType, const vector<double>& coefs, double eps=1e-12): 
      _kernelType(kernelType), _coefs(coefs), _mindif(eps) {
    parse_kernel_enum(_kernelType, _equation_type, _kernel_type, _kernel_var);
    initialize_degrees_of_freedom();
  }

  Kernel3d(const Kernel3d& c): _kernelType(c._kernelType), _coefs(c._coefs),_mindif(c._mindif) {
    parse_kernel_enum(_kernelType, _equation_type, _kernel_type, _kernel_var);
    initialize_degrees_of_freedom();
  }

  Kernel3d& operator=(const Kernel3d& c) {
      _kernelType = c._kernelType;
      _coefs = c._coefs;
      _mindif= c._mindif;

      parse_kernel_enum(_kernelType, _equation_type, _kernel_type, _kernel_var);
      initialize_degrees_of_freedom();

      return *this;
  }
  int& kernelType()                          { return _kernelType; }
  const int& kernelType() const              { return _kernelType; }
  vector<double>& coefs()             { return _coefs; }
  const vector<double>& coefs() const { return _coefs; }
  void set_eps(double eps)  { _mindif = eps; }
  double coefs(int i) { return _coefs[i]; }
  int dim() { return 3; }
  uint srcDOF() const;
  uint trgDOF() const;
  /*!
   * Utility function to strip out Equation_type, Kernel_type, and
   * Kernel_variable enums above from the anonymous enum type fully describing
   * the kernel:
   * @param int     kernel_enum         Full anonymous enum describing kernel
   * @param int&    equation_type       Return variable storing the value of
   *                                    the Equation_type (100, 200, 300, etc.)
   * @param int&    kernel_type         Return variable storing the value of
   *                                    the Kernel_type (10, 20, etc.)
   * @param int&    kernel_variable     Return variable storing the value of
   *                                    the Kernel_variable (10, 20)
   */
  static void parse_kernel_enum(int kernel_enum, Equation_type& equation_type, 
          Kernel_type& kernel_type, Kernel_variable& kernel_variable){

      // pull out the 1's digit
      kernel_variable = (Kernel_variable) (kernel_enum % 10); 
      // pull out the 10's digit
      kernel_type =  (Kernel_type) ((kernel_enum % 100) - kernel_variable);
      // pull  out the 100's digits
      int eq = kernel_enum - kernel_type - kernel_variable;  
      equation_type = (Equation_type) eq; 
  }


  /**
   * Return an n x n identity matrix as a DblNumMat.
   * @param int     n                   number of rows & cols of the identity
   *                                    matrix
   */
  static DblNumMat get_identity_matrix(int n);

  /**
   * Return the 3x3 cross product matrix corresponding to the skew-symmetric
   * matrix expression of the usual cross product:
   *
   * a x b = C*b =  [0   -a_3   a_2] [b_1]
   *                [a_3    0  -a_1] [b_2]
   *                [-a_2 a_1     0] [b_3]
   *
   * where a = target - source
   *
   * @param DblNumVec   source           
   * @param DblNumVec   target            
   */
  static DblNumMat get_cross_product_matrix(DblNumVec source, DblNumVec target);


  /**
   * Indicates whether a particular equation type associated with a given
   * kernel requires modifications in the solver due to an incompressibility
   * condition.
   */
  bool is_incompressible(){
      switch(_equation_type){
          case STOKES:
              return true;
          default:
              return false;
      }
  }

  void initialize_degrees_of_freedom(){
    if (_equation_type == LAPLACE || _equation_type == MOD_HELMHOLTZ){
        _tdof = 1;
        _sdof = 1;

        // For Laplace and Helmholtz, there are no additional degrees of 
        // freedom for rotation and multiply connected boundaries; 0 them
        // out so things break if they're inappropriately used.
        _rdof = 0;
    
    } else if (_equation_type == STOKES){
        if(_kernel_var == VAR_P){
            _tdof = 1; // TODO MJM this is bug. Pressure eval with blended surfaces 
                        // doesn't initialize this properly inside core_evaluation
                        // However, trgDOF does work, weirdly.
        } else{
            _tdof = 3;
        }
        _sdof = 3;
        _rdof = 3;
    } else if (_equation_type == NAVIER){
        _tdof = 3;
        _sdof = 3;
        _rdof = 3;
    } else{
      cerr << "KernelNotYetImplementedError: " << (int) kernelType() <<", Equation: " 
          << (int) _equation_type << endl;
      ebiAssert(false);
    }
    _gdof = _sdof;
    _pdof = _rdof + _gdof;
  }


  void dirichlet_bc_from_singularities(
          DblNumMat source_positions,
          DblNumMat source_strengths,
          DblNumMat target_positions, 
          DblNumMat& boundary_data);
  
  void neumann_bc_from_singularities(
          DblNumMat source_positions,
          DblNumMat source_strengths,
          DblNumMat target_positions, 
          DblNumMat target_normals, 
          DblNumMat& boundary_data);
  void laplace_dirichlet(DblNumMat source_positions,
          DblNumMat source_strengths,
          DblNumMat target_positions, 
          DblNumMat& boundary_data);
  
  void laplace_neumann(DblNumMat source_positions,
          DblNumMat source_strengths,
          DblNumMat target_positions, 
          DblNumMat target_normals, 
          DblNumMat& boundary_data);
  void navier_dirichlet(DblNumMat source_positions,
          DblNumMat source_strengths,
          DblNumMat target_positions, 
          DblNumMat& boundary_data);
  
  void navier_neumann(DblNumMat source_positions,
          DblNumMat source_strengths,
          DblNumMat target_positions, 
          DblNumMat target_normals, 
          DblNumMat& boundary_data);

  double density_unscaling_value(double position_scale);
  //---------------------------------------------------------------------------
  // Accessors
  int get_tdof() {
    return _tdof;
  }

  int get_sdof() {
    return _sdof;
  }

  int get_rdof() {
    return _rdof;
  }

  int get_gdof() {
    return _gdof;
  }

  int get_pdof() {
    return _pdof;
  }

  Equation_type equation_type(){
    return _equation_type;
  }

  Kernel_type  kernel_type(){
    return _kernel_type;
  }

  Kernel_variable kernel_var(){
    return _kernel_var;
  }


  ///homogeneous or not
  bool homogeneous() const;

  //---------------------------------------------------------------------------
  
  ///homogeneous degree, vector size == sourceDegreeOfFreedom
  void homogeneousDeg(vector<double>&) const; 

  //each kernelType handles coef in its own way, can be empty
  int kernel(const DblNumMat& srcPos, const DblNumMat& srcNor, const DblNumMat& trgPos, DblNumMat& inter);
};

inline bool operator==(const Kernel3d& a, const Kernel3d& b) {
  return (a.kernelType()==b.kernelType() && a.coefs()==b.coefs());
}

END_EBI_NAMESPACE

#endif
