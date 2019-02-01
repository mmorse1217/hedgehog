#ifndef _MOBO_BLAS_H_
#define _MOBO_BLAS_H_

#include "ebi_namespace.hpp"
#include "petsc.h"
extern "C"{
} 

BEGIN_EBI_NAMESPACE

/*! Allows for linking to external fortran blas libraries for DGEMM */
#define DGEMM dgemm_
/*! Allows for linking to external fortran blas libraries for DGEMV */
#define DGEMV dgemv_
/*! Allows for linking to external fortran blas libraries for DAXPY */
#define DAXPY daxpy_
/*! Allows for linking to external fortran blas libraries for DGER */
#define DGER  dger_
/*! Allows for linking to external fortran blas libraries for DSCAL */
#define DSCAL dscal_
/*! Allows for linking to external fortran blas libraries for DSCAL */
#define DCOPY dcopy_

extern "C"
{
  /*! DAXPY compute y := alpha * x + y where alpha is a scalar and x and y are n-vectors.
	*  See http://www.netlib.org/blas/daxpy.f for more information.
	*/
  void DAXPY(PetscBLASInt* N, double* ALPHA, double* X, PetscBLASInt* INCX, double* Y, PetscBLASInt* INCY);
  /*!  DGEMM  performs one of the matrix-matrix operations
	*
	*     C := alpha*op( A )*op( B ) + beta*C,
	*
	*  where  op( X ) is one of
	*
	*     op( X ) = X   or   op( X ) = X',
	*
	*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
	*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
	*  See http://www.netlib.org/blas/dgemm.f for more information.
	*/
  void DGEMM(char* TRANSA, char* TRANSB, PetscBLASInt* M, PetscBLASInt* N, PetscBLASInt* K, double* ALPHA, double* A, PetscBLASInt* LDA, double* B, PetscBLASInt* LDB, double* BETA, double* C, PetscBLASInt* LDC);
  /*!  DGEMV  performs one of the matrix-vector operations
	*
	*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
	*
	*  where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
	*  See http://www.netlib.org/blas/dgemv.f for more information
	*/
  void DGEMV(char* TRANS, PetscBLASInt* M, PetscBLASInt* N, double* ALPHA, double* A, PetscBLASInt* LDA, double* X, PetscBLASInt* INCX, double* BETA, double* Y, PetscBLASInt* INCY);

  /*!  DGER   performs the rank 1 operation
	*
	*     A := alpha*x*y' + A,
	*
	*  where alpha is a scalar, x is an m element vector, y is an n element
	*  vector and A is an m by n matrix.
	*  See http://www.netlib.org/blas/dger.f for more information
	*/
  void DGER (PetscBLASInt* M, PetscBLASInt * N, double* ALPHA, double* X, PetscBLASInt* INCX, double* Y, PetscBLASInt* INCY, double* A, PetscBLASInt* LDA);
  /*! DSCAL computes y := alpha * y where alpha is a scalar and y is an n-vector.
	*  See http://www.netlib.org/blas/dscal.f for more information
	*/
  void DSCAL(PetscBLASInt* N, double* ALPHA, double* X, PetscBLASInt* INCX);

  /*! DCOPY copies x to y where x and y are n-vectors.
	*/
  void DCOPY(PetscBLASInt* N, double* X, PetscBLASInt* INCX, double* Y, PetscBLASInt* INCY);
  
}

END_EBI_NAMESPACE

#endif

