#ifndef _MOBO_LAPACK_H_
#define _MOBO_LAPACK_H_

#include "ebi_namespace.hpp"
#include "petsc.h"
extern "C"{
}

//blas and lapack workspace for the code

BEGIN_EBI_NAMESPACE

//blas and lapack workspace for the code
/*! Allows for linking to external fortran lapack libraries for DGESVD */
#define DGESVD dgesvd_
/*! Allows for linking to external fortran lapack libraries for DGESDD */
#define DGESDD dgesdd_
/*! Allows for linking to external fortran lapack libraries for DGETRF*/
#define DGETRF dgetrf_
/*! Allows for linking to external fortran lapack libraries for DGETRI */
#define DGETRI dgetri_

#define DCUHRE dcuhre_

//#define QSHEP3 qshep3_
//#define QS3VAL qs3val_

//EXTERN_C_BEGIN
extern "C"
{
    /*!	 DGESVD computes the singular value decomposition (SVD) of a real
	*  M-by-N matrix A, optionally computing the left and/or right singular
	*  vectors. The SVD is written
	*
	*       A = U * SIGMA * transpose(V)
	*
	*  where SIGMA is an M-by-N matrix which is zero except for its
	*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
	*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
	*  are the singular values of A; they are real and non-negative, and
	*  are returned in descending order.  The first min(m,n) columns of
	*  U and V are the left and right singular vectors of A.
	*
	* See http://www.netlib.org/lapack/double/dgesvd.f for more information 
	*/
  //extern void DGESVD(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *INFO);
  extern void DGESVD(char *jobu, char *jobvt, PetscBLASInt *m, PetscBLASInt *n, double *a, PetscBLASInt *lda, double *s, double *u, PetscBLASInt *ldu, double *vt,  PetscBLASInt *ldvt, double *work, PetscBLASInt *lwork, PetscBLASInt *info);
  /*! DGESDD computes the singular value decomposition (SVD) of a real
	*  M-by-N matrix A, optionally computing the left and right singular
	*  vectors.  If singular vectors are desired, it uses a
	* divide-and-conquer algorithm.
	*
	*  The SVD is written
	*
	*       A = U * SIGMA * transpose(V)
	*
	*  where SIGMA is an M-by-N matrix which is zero except for its
	*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
	*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
`	*  are the singular values of A; they are real and non-negative, and
	*  are returned in descending order.  The first min(m,n) columns of
	*  U and V are the left and right singular vectors of A.
	*
	*  See http://www.netlib.org/lapack/double/dgesdd.f for more information 
	*/
  extern void DGESDD(char *jobz, PetscBLASInt* m, PetscBLASInt* n, double* a, PetscBLASInt* lda, double* s, double* u, PetscBLASInt* ldu, double* vt, PetscBLASInt* ldvt, double* work, PetscBLASInt* lwork, PetscBLASInt* iwork, PetscBLASInt* info);
  /*!  DGETRF computes an LU factorization of a general M-by-N matrix A
	*  using partial pivoting with row interchanges.
	*
	*  The factorization has the form
	*
	*	     A = P * L * U
	*
	*  where P is a permutation matrix, L is lower triangular with unit
	*  diagonal elements (lower trapezoidal if m > n), and U is upper
	*  triangular (upper trapezoidal if m < n).
	*
	*  See http://www.netlib.org/lapack/double/dgetrf.f for more information
	*/
  extern void DGETRF(PetscBLASInt *M, PetscBLASInt *N, double *A, PetscBLASInt *LDA, PetscBLASInt *IPIV, PetscBLASInt *INFO);
  /*!  DGETRI computes the inverse of a matrix using the LU factorization
	*  computed by DGETRF.
	*
	*  This method inverts U and then computes inv(A) by solving the system
	*  inv(A)*L = inv(U) for inv(A).
	*
	*  See http://www.netlib.org/lapack/double/dgetri.f for more information
	*/
  extern void DGETRI(PetscBLASInt *N, double *A, PetscBLASInt *LDA, PetscBLASInt *IPIV, double *WORK, PetscBLASInt *LWORK, PetscBLASInt *INFO);


  extern	void DCUHRE(long *ndim,long *numfun,double *a,double *b, 
							long *minpts,long *maxpts,void *funsub,
							double *epsabs,double *epsrel,long *key, 
							long *nw,long *restar,double *result,
							double *abserr,long *neval,long *ifail,
							double *work);

  //  extern void QSHEP3(long int *N, float* X, float *Y, float *Z, float *F, long int *NQ, long int *NW, long int *NR,
  //					long int *LCELL, long int *LNEXT, float *XYZMIN, float *XYZDEL, float *RMAX, float *RSQ,
  //					float *A, long int *IER);

  //  extern double QS3VAL(float *PX, float *PY, float *PZ, long int *N, float *X, float *Y, float *Z, float *F, long int *NR,
  //					 long int *LCELL, long int *LNEXT, float *XYZMIN, float *XYZDEL, float *RMAX, float *RSQ,
  //					 float *A);

}
//EXTERN_C_END
  /*
void DBDSQR(char *UPLO, int *N, int *NCVT, int *NRU, int *NCC, double *D, double *E, double *VT, int *LDVT, double *U, int *LDU, double *C, int *LDC, double *WORK, int *INFO);
void DDISNA(char *JOB, int *M, int *N, double *D, double *SEP, int *INFO);
void DGBBRD(char *VECT, int *M, int *N, int *NCC, int *KL, int *KU, double *AB, int *LDAB, double *D, double *E, double *Q, int *LDQ, double *PT, int *LDPT, double *C, int *LDC, double *WORK, int *INFO);
void DGBCON(char *NORM, int *N, int *KL, int *KU, double *AB, int *LDAB, int *IPIV, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DGBEQU(int *M, int *N, int *KL, int *KU, double *AB, int *LDAB, double *R, double *C, double *ROWCND, double *COLCND, double *AMAX, int *INFO);
void DGBRFS(char *TRANS, int *N, int *KL, int *KU, int *NRHS, double *AB, int *LDAB, double *AFB, int *LDAFB, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DGBSV(int *N, int *KL, int *KU, int *NRHS, double *AB, int *LDAB, int *IPIV, double *B, int *LDB, int *INFO);
void DGBSVX(char *FACT, char *TRANS, int *N, int *KL, int *KU, int *NRHS, double *AB, int *LDAB, double *AFB, int *LDAFB, int *IPIV, char *EQUED, double *R, double *C, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DGBTF2(int *M, int *N, int *KL, int *KU, double *AB, int *LDAB, int *IPIV, int *INFO);
void DGBTRF(int *M, int *N, int *KL, int *KU, double *AB, int *LDAB, int *IPIV, int *INFO);
void DGBTRS(char *TRANS, int *N, int *KL, int *KU, int *NRHS, double *AB, int *LDAB, int *IPIV, double *B, int *LDB, int *INFO);
void DGEBAK(char *JOB, char *SIDE, int *N, int *ILO, int *IHI, double *SCALE, int *M, double *V, int *LDV, int *INFO);
void DGEBAL(char *JOB, int *N, double *A, int *LDA, int *ILO, int *IHI, double *SCALE, int *INFO);
void DGEBD2(int *M, int *N, double *A, int *LDA, double *D, double *E, double *TAUQ, double *TAUP, double *WORK, int *INFO);
void DGEBRD(int *M, int *N, double *A, int *LDA, double *D, double *E, double *TAUQ, double *TAUP, double *WORK, int *LWORK, int *INFO);
void DGECON(char *NORM, int *N, double *A, int *LDA, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DGEEQU(int *M, int *N, double *A, int *LDA, double *R, double *C, double *ROWCND, double *COLCND, double *AMAX, int *INFO);
void DGEES(char *JOBVS, char *SORT, int *SELECT, int *N, double *A, int *LDA, int *SDIM, double *WR, double *WI, double *VS, int *LDVS, double *WORK, int *LWORK, int *BWORK, int *INFO);
void DGEESX(char *JOBVS, char *SORT, int *SELECT, char *SENSE, int *N, double *A, int *LDA, int *SDIM, double *WR, double *WI, double *VS, int *LDVS, double *RCONDE, double *RCONDV, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *BWORK, int *INFO);
void DGEEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI, double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO);
void DGEEVX(char *BALANC, char *JOBVL, char *JOBVR, char *SENSE, int *N, double *A, int *LDA, double *WR, double *WI, double *VL, int *LDVL, double *VR, int *LDVR, int *ILO, int *IHI, double *SCALE, double *ABNRM, double *RCONDE, double *RCONDV, double *WORK, int *LWORK, int *IWORK, int *INFO);
void DGEGS(char *JOBVSL, char *JOBVSR, int *N, double *A, int *LDA, double *B, int *LDB, double *ALPHAR, double *ALPHAI, double *BETA, double *VSL, int *LDVSL, double *VSR, int *LDVSR, double *WORK, int *LWORK, int *INFO);
void DGEGV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *B, int *LDB, double *ALPHAR, double *ALPHAI, double *BETA, double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO);
void DGEHD2(int *N, int *ILO, int *IHI, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DGEHRD(int *N, int *ILO, int *IHI, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DGELQ2(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DGELQF(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DGELS(char *TRANS, int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
void DGELSS(int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *S, double *RCOND, int *RANK, double *WORK, int *LWORK, int *INFO);
void DGELSX(int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *JPVT, double *RCOND, int *RANK, double *WORK, int *INFO);
void DGEQL2(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DGEQLF(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DGEQPF(int *M, int *N, double *A, int *LDA, int *JPVT, double *TAU, double *WORK, int *INFO);
void DGEQR2(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DGEQRF(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DGERFS(char *TRANS, int *N, int *NRHS, double *A, int *LDA, double *AF, int *LDAF, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DGERQ2(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DGERQF(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DGESV(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
void DGESVD(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *INFO);
void DGESVX(char *FACT, char *TRANS, int *N, int *NRHS, double *A, int *LDA, double *AF, int *LDAF, int *IPIV, char *EQUED, double *R, double *C, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DGETF2(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void DGETRF(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void DGETRI(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
void DGETRS(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
void DGGBAK(char *JOB, char *SIDE, int *N, int *ILO, int *IHI, double *LSCALE, double *RSCALE, int *M, double *V, int *LDV, int *INFO);
void DGGBAL(char *JOB, int *N, double *A, int *LDA, double *B, int *LDB, int *ILO, int *IHI, double *LSCALE, double *RSCALE, double *WORK, int *INFO);
void DGGGLM(int *N, int *M, int *P, double *A, int *LDA, double *B, int *LDB, double *D, double *X, double *Y, double *WORK, int *LWORK, int *INFO);
void DGGHRD(char *COMPQ, char *COMPZ, int *N, int *ILO, int *IHI, double *A, int *LDA, double *B, int *LDB, double *Q, int *LDQ, double *Z, int *LDZ, int *INFO);
void DGGLSE(int *M, int *N, int *P, double *A, int *LDA, double *B, int *LDB, double *C, double *D, double *X, double *WORK, int *LWORK, int *INFO);
void DGGQRF(int *N, int *M, int *P, double *A, int *LDA, double *TAUA, double *B, int *LDB, double *TAUB, double *WORK, int *LWORK, int *INFO);
void DGGRQF(int *M, int *P, int *N, double *A, int *LDA, double *TAUA, double *B, int *LDB, double *TAUB, double *WORK, int *LWORK, int *INFO);
void DGGSVD(char *JOBU, char *JOBV, char *JOBQ, int *M, int *N, int *P, int *K, int *L, double *A, int *LDA, double *B, int *LDB, double *ALPHA, double *BETA, double *U, int *LDU, double *V, int *LDV, double *Q, int *LDQ, double *WORK, int *IWORK, int *INFO);
void DGGSVP(char *JOBU, char *JOBV, char *JOBQ, int *M, int *P, int *N, double *A, int *LDA, double *B, int *LDB, double *TOLA, double *TOLB, int *K, int *L, double *U, int *LDU, double *V, int *LDV, double *Q, int *LDQ, int *IWORK, double *TAU, double *WORK, int *INFO);
void DGTCON(char *NORM, int *N, double *DL, double *D, double *DU, double *DU2, int *IPIV, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DGTRFS(char *TRANS, int *N, int *NRHS, double *DL, double *D, double *DU, double *DLF, double *DF, double *DUF, double *DU2, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DGTSV(int *N, int *NRHS, double *DL, double *D, double *DU, double *B, int *LDB, int *INFO);
void DGTSVX(char *FACT, char *TRANS, int *N, int *NRHS, double *DL, double *D, double *DU, double *DLF, double *DF, double *DUF, double *DU2, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DGTTRF(int *N, double *DL, double *D, double *DU, double *DU2, int *IPIV, int *INFO);
void DGTTRS(char *TRANS, int *N, int *NRHS, double *DL, double *D, double *DU, double *DU2, int *IPIV, double *B, int *LDB, int *INFO);
void DHGEQZ(char *JOB, char *COMPQ, char *COMPZ, int *N, int *ILO, int *IHI, double *A, int *LDA, double *B, int *LDB, double *ALPHAR, double *ALPHAI, double *BETA, double *Q, int *LDQ, double *Z, int *LDZ, double *WORK, int *LWORK, int *INFO);
void DHSEIN(char *SIDE, char *EIGSRC, char *INITV, int *SELECT, int *N, double *H, int *LDH, double *WR, double *WI, double *VL, int *LDVL, double *VR, int *LDVR, int *MM, int *M, double *WORK, int *IFAILL, int *IFAILR, int *INFO);
void DHSEQR(char *JOB, char *COMPZ, int *N, int *ILO, int *IHI, double *H, int *LDH, double *WR, double *WI, double *Z, int *LDZ, double *WORK, int *LWORK, int *INFO);
void DLAED0(int *ICOMPQ, int *QSIZ, int *N, double *D, double *E, double *Q, int *LDQ, double *QSTORE, int *LDQS, double *WORK, int *IWORK, int *INFO);
void DLAED1(int *N, double *D, double *Q, int *LDQ, int *INDXQ, double *RHO, int *CUTPNT, double *WORK, int *IWORK, int *INFO);
void DLAED2(int *K, int *N, double *D, double *Q, int *LDQ, int *INDXQ, double *RHO, int *CUTPNT, double *Z, double *DLAMDA, double *Q2, int *LDQ2, int *INDXC, double *W, int *INDXP, int *INDX, int *COLTYP, int *INFO);
void DLAED3(int *K, int *KSTART, int *KSTOP, int *N, double *D, double *Q, int *LDQ, double *RHO, int *CUTPNT, double *DLAMDA, double *Q2, int *LDQ2, int *INDXC, int *CTOT, double *W, double *S, int *LDS, int *INFO);
void DLAED4(int *N, int *I, double *D, double *Z, double *DELTA, double *RHO, double *DLAM, int *INFO);
void DLAED5(int *I, double *D, double *Z, double *DELTA, double *RHO, double *DLAM);
void DLAED6(int *KNITER, int *ORGATI, double *RHO, double *D, double *Z, double *FINIT, double *TAU, int *INFO);
void DLAED7(int *ICOMPQ, int *N, int *QSIZ, int *TLVLS, int *CURLVL, int *CURPBM, double *D, double *Q, int *LDQ, int *INDXQ, double *RHO, int *CUTPNT, double *QSTORE, int *QPTR, int *PRMPTR, int *PERM, int *GIVPTR, int *GIVCOL, double *GIVNUM, double *WORK, int *IWORK, int *INFO);
void DLAED8(int *ICOMPQ, int *K, int *N, int *QSIZ, double *D, double *Q, int *LDQ, int *INDXQ, double *RHO, int *CUTPNT, double *Z, double *DLAMDA, double *Q2, int *LDQ2, double *W, int *PERM, int *GIVPTR, int *GIVCOL, double *GIVNUM, int *INDXP, int *INDX, int *INFO);
void DLAED9(int *K, int *KSTART, int *KSTOP, int *N, double *D, double *Q, int *LDQ, double *RHO, double *DLAMDA, double *W, double *S, int *LDS, int *INFO);
void DLAEDA(int *N, int *TLVLS, int *CURLVL, int *CURPBM, int *PRMPTR, int *PERM, int *GIVPTR, int *GIVCOL, double *GIVNUM, double *Q, int *QPTR, double *Z, double *ZTEMP, int *INFO);
void DLAGTF(int *N, double *A, double *LAMBDA, double *B, double *C, double *TOL, double *D, int *IN, int *INFO);
void DLAMRG(int *N1, int *N2, double *A, int *DTRD1, int *DTRD2, int *INDEX);
void DLASQ1(int *N, double *D, double *E, double *WORK, int *INFO);
void DLASQ2(int *M, double *Q, double *E, double *QQ, double *EE, double *EPS, double *TOL2, double *SMALL2, double *SUP, int *KEND, int *INFO);
void DLASQ3(int *N, double *Q, double *E, double *QQ, double *EE, double *SUP, double *SIGMA, int *KEND, int *OFF, int *IPHASE, int *ICONV, double *EPS, double *TOL2, double *SMALL2);
void DLASQ4(int *N, double *Q, double *E, double *TAU, double *SUP);
void DLASRT(char *ID, int *N, double *D, int *INFO);
void DLATZM(char *SIDE, int *M, int *N, double *V, int *INCV, double *TAU, double *C1, double *C2, int *LDC, double *WORK);
void DOPGTR(char *UPLO, int *N, double *AP, double *TAU, double *Q, int *LDQ, double *WORK, int *INFO);
void DOPMTR(char *SIDE, char *UPLO, char *TRANS, int *M, int *N, double *AP, double *TAU, double *C, int *LDC, double *WORK, int *INFO);
void DORG2L(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DORG2R(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DORGBR(char *VECT, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DORGHR(int *N, int *ILO, int *IHI, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DORGL2(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DORGLQ(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DORGQL(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DORGQR(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DORGR2(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *INFO);
void DORGRQ(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DORGTR(char *UPLO, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void DORM2L(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *INFO);
void DORM2R(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *INFO);
void DORMBR(char *VECT, char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO);
void DORMHR(char *SIDE, char *TRANS, int *M, int *N, int *ILO, int *IHI, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO);
void DORML2(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *INFO);
void DORMLQ(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO);
void DORMQL(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO);
void DORMQR(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO);
void DORMR2(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *INFO);
void DORMRQ(char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO);
void DORMTR(char *SIDE, char *UPLO, char *TRANS, int *M, int *N, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO);
void DPBCON(char *UPLO, int *N, int *KD, double *AB, int *LDAB, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DPBEQU(char *UPLO, int *N, int *KD, double *AB, int *LDAB, double *S, double *SCOND, double *AMAX, int *INFO);
void DPBRFS(char *UPLO, int *N, int *KD, int *NRHS, double *AB, int *LDAB, double *AFB, int *LDAFB, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DPBSTF(char *UPLO, int *N, int *KD, double *AB, int *LDAB, int *INFO);
void DPBSV(char *UPLO, int *N, int *KD, int *NRHS, double *AB, int *LDAB, double *B, int *LDB, int *INFO);
void DPBSVX(char *FACT, char *UPLO, int *N, int *KD, int *NRHS, double *AB, int *LDAB, double *AFB, int *LDAFB, char *EQUED, double *S, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DPBTF2(char *UPLO, int *N, int *KD, double *AB, int *LDAB, int *INFO);
void DPBTRF(char *UPLO, int *N, int *KD, double *AB, int *LDAB, int *INFO);
void DPBTRS(char *UPLO, int *N, int *KD, int *NRHS, double *AB, int *LDAB, double *B, int *LDB, int *INFO);
void DPOCON(char *UPLO, int *N, double *A, int *LDA, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DPOEQU(int *N, double *A, int *LDA, double *S, double *SCOND, double *AMAX, int *INFO);
void DPORFS(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *AF, int *LDAF, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DPOSV(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);
void DPOSVX(char *FACT, char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *AF, int *LDAF, char *EQUED, double *S, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DPOTF2(char *UPLO, int *N, double *A, int *LDA, int *INFO);
void DPOTRF(char *UPLO, int *N, double *A, int *LDA, int *INFO);
void DPOTRI(char *UPLO, int *N, double *A, int *LDA, int *INFO);
void DPOTRS(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);
void DPPCON(char *UPLO, int *N, double *AP, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DPPEQU(char *UPLO, int *N, double *AP, double *S, double *SCOND, double *AMAX, int *INFO);
void DPPRFS(char *UPLO, int *N, int *NRHS, double *AP, double *AFP, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DPPSV(char *UPLO, int *N, int *NRHS, double *AP, double *B, int *LDB, int *INFO);
void DPPSVX(char *FACT, char *UPLO, int *N, int *NRHS, double *AP, double *AFP, char *EQUED, double *S, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DPPTRF(char *UPLO, int *N, double *AP, int *INFO);
void DPPTRI(char *UPLO, int *N, double *AP, int *INFO);
void DPPTRS(char *UPLO, int *N, int *NRHS, double *AP, double *B, int *LDB, int *INFO);
void DPTCON(int *N, double *D, double *E, double *ANORM, double *RCOND, double *WORK, int *INFO);
void DPTEQR(char *COMPZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *INFO);
void DPTRFS(int *N, int *NRHS, double *D, double *E, double *DF, double *EF, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *INFO);
void DPTSV(int *N, int *NRHS, double *D, double *E, double *B, int *LDB, int *INFO);
void DPTSVX(char *FACT, int *N, int *NRHS, double *D, double *E, double *DF, double *EF, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *INFO);
void DPTTRF(int *N, double *D, double *E, int *INFO);
void DPTTRS(int *N, int *NRHS, double *D, double *E, double *B, int *LDB, int *INFO);
void DSBEV(char *JOBZ, char *UPLO, int *N, int *KD, double *AB, int *LDAB, double *W, double *Z, int *LDZ, double *WORK, int *INFO);
void DSBEVD(char *JOBZ, char *UPLO, int *N, int *KD, double *AB, int *LDAB, double *W, double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);
void DSBEVX(char *JOBZ, char *RANGE, char *UPLO, int *N, int *KD, double *AB, int *LDAB, double *Q, int *LDQ, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, double *WORK, int *IWORK, int *IFAIL, int *INFO);
void DSBGST(char *VECT, char *UPLO, int *N, int *KA, int *KB, double *AB, int *LDAB, double *BB, int *LDBB, double *X, int *LDX, double *WORK, int *INFO);
void DSBGV(char *JOBZ, char *UPLO, int *N, int *KA, int *KB, double *AB, int *LDAB, double *BB, int *LDBB, double *W, double *Z, int *LDZ, double *WORK, int *INFO);
void DSBTRD(char *VECT, char *UPLO, int *N, int *KD, double *AB, int *LDAB, double *D, double *E, double *Q, int *LDQ, double *WORK, int *INFO);
void DSPCON(char *UPLO, int *N, double *AP, int *IPIV, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DSPEV(char *JOBZ, char *UPLO, int *N, double *AP, double *W, double *Z, int *LDZ, double *WORK, int *INFO);
void DSPEVD(char *JOBZ, char *UPLO, int *N, double *AP, double *W, double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);
void DSPEVX(char *JOBZ, char *RANGE, char *UPLO, int *N, double *AP, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, double *WORK, int *IWORK, int *IFAIL, int *INFO);
void DSPGST(int *ITYPE, char *UPLO, int *N, double *AP, double *BP, int *INFO);
void DSPGV(int *ITYPE, char *JOBZ, char *UPLO, int *N, double *AP, double *BP, double *W, double *Z, int *LDZ, double *WORK, int *INFO);
void DSPRFS(char *UPLO, int *N, int *NRHS, double *AP, double *AFP, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DSPSV(char *UPLO, int *N, int *NRHS, double *AP, int *IPIV, double *B, int *LDB, int *INFO);
void DSPSVX(char *FACT, char *UPLO, int *N, int *NRHS, double *AP, double *AFP, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DSPTRD(char *UPLO, int *N, double *AP, double *D, double *E, double *TAU, int *INFO);
void DSPTRF(char *UPLO, int *N, double *AP, int *IPIV, int *INFO);
void DSPTRI(char *UPLO, int *N, double *AP, int *IPIV, double *WORK, int *INFO);
void DSPTRS(char *UPLO, int *N, int *NRHS, double *AP, int *IPIV, double *B, int *LDB, int *INFO);
void DSTEBZ(char *RANGE, char *ORDER, int *N, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, double *D, double *E, int *M, int *NSPLIT, double *W, int *IBLOCK, int *ISPLIT, double *WORK, int *IWORK, int *INFO);
void DSTEDC(char *COMPZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);
void DSTEIN(int *N, double *D, double *E, int *M, double *W, int *IBLOCK, int *ISPLIT, double *Z, int *LDZ, double *WORK, int *IWORK, int *IFAIL, int *INFO);
void DSTEQR(char *COMPZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *INFO);
void DSTERF(int *N, double *D, double *E, int *INFO);
void DSTEV(char *JOBZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *INFO);
void DSTEVD(char *JOBZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);
void DSTEVX(char *JOBZ, char *RANGE, int *N, double *D, double *E, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, double *WORK, int *IWORK, int *IFAIL, int *INFO);
void DSYCON(char *UPLO, int *N, double *A, int *LDA, int *IPIV, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DSYEV(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);
void DSYEVD(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);
void DSYEVX(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *IFAIL, int *INFO);
void DSYGS2(int *ITYPE, char *UPLO, int *N, double *A, int *LDA, double *B, int *LDB, int *INFO);
void DSYGST(int *ITYPE, char *UPLO, int *N, double *A, int *LDA, double *B, int *LDB, int *INFO);
void DSYGV(int *ITYPE, char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *B, int *LDB, double *W, double *WORK, int *LWORK, int *INFO);
void DSYRFS(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *AF, int *LDAF, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DSYSV(char *UPLO, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
void DSYSVX(char *FACT, char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *AF, int *LDAF, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, int *LWORK, int *IWORK, int *INFO);
void DSYTD2(char *UPLO, int *N, double *A, int *LDA, double *D, double *E, double *TAU, int *INFO);
void DSYTF2(char *UPLO, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void DSYTRD(char *UPLO, int *N, double *A, int *LDA, double *D, double *E, double *TAU, double *WORK, int *LWORK, int *INFO);
void DSYTRF(char *UPLO, int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
void DSYTRI(char *UPLO, int *N, double *A, int *LDA, int *IPIV, double *WORK, int *INFO);
void DSYTRS(char *UPLO, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
void DTBCON(char *NORM, char *UPLO, char *DIAG, int *N, int *KD, double *AB, int *LDAB, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DTBRFS(char *UPLO, char *TRANS, char *DIAG, int *N, int *KD, int *NRHS, double *AB, int *LDAB, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DTBTRS(char *UPLO, char *TRANS, char *DIAG, int *N, int *KD, int *NRHS, double *AB, int *LDAB, double *B, int *LDB, int *INFO);
void DTGEVC(char *SIDE, char *HOWMNY, int *SELECT, int *N, double *A, int *LDA, double *B, int *LDB, double *VL, int *LDVL, double *VR, int *LDVR, int *MM, int *M, double *WORK, int *INFO);
void DTGSJA(char *JOBU, char *JOBV, char *JOBQ, int *M, int *P, int *N, int *K, int *L, double *A, int *LDA, double *B, int *LDB, double *TOLA, double *TOLB, double *ALPHA, double *BETA, double *U, int *LDU, double *V, int *LDV, double *Q, int *LDQ, double *WORK, int *NCYCLE, int *INFO);
void DTPCON(char *NORM, char *UPLO, char *DIAG, int *N, double *AP, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DTPRFS(char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, double *AP, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DTPTRI(char *UPLO, char *DIAG, int *N, double *AP, int *INFO);
void DTPTRS(char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, double *AP, double *B, int *LDB, int *INFO);
void DTRCON(char *NORM, char *UPLO, char *DIAG, int *N, double *A, int *LDA, double *RCOND, double *WORK, int *IWORK, int *INFO);
void DTREVC(char *SIDE, char *HOWMNY, int *SELECT, int *N, double *T, int *LDT, double *VL, int *LDVL, double *VR, int *LDVR, int *MM, int *M, double *WORK, int *INFO);
void DTREXC(char *COMPQ, int *N, double *T, int *LDT, double *Q, int *LDQ, int *IFST, int *ILST, double *WORK, int *INFO);
void DTRRFS(char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO);
void DTRSEN(char *JOB, char *COMPQ, int *SELECT, int *N, double *T, int *LDT, double *Q, int *LDQ, double *WR, double *WI, int *M, double *S, double *SEP, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);
void DTRSNA(char *JOB, char *HOWMNY, int *SELECT, int *N, double *T, int *LDT, double *VL, int *LDVL, double *VR, int *LDVR, double *S, double *SEP, int *MM, int *M, double *WORK, int *LDWORK, int *IWORK, int *INFO);
void DTRSYL(char *TRANA, char *TRANB, int *ISGN, int *M, int *N, double *A, int *LDA, double *B, int *LDB, double *C, int *LDC, double *SCALE, int *INFO);
void DTRTI2(char *UPLO, char *DIAG, int *N, double *A, int *LDA, int *INFO);
void DTRTRI(char *UPLO, char *DIAG, int *N, double *A, int *LDA, int *INFO);
void DTRTRS(char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);
void DTZRQF(int *M, int *N, double *A, int *LDA, double *TAU, int *INFO);
  */

END_EBI_NAMESPACE

#endif
