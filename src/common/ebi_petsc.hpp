/*! \file */
#ifndef _EBI_HPP_EBI_
#define _EBI_HPP_EBI_

#include "petscsnes.h"
extern "C" {
}
#include "ebi_namespace.hpp"

BEGIN_EBI_NAMESPACE

/* ********************************************************************** */

#define STRING_LENGTH 1024
#define EbiNOERROR 0
#define __MODULENAME__ "all"
static const int DIM=3;
typedef void (*FunctionHandle)(Vec, int, Vec&);
/// All error code returing functions should use ebiEC as their return type.
typedef int ebiEC;

/* ********************************************************************** */
struct ebiGlobals{
  // used with the PETSc Vec operations
  static const double zero;
  static const double one;
  static const double mone;  
  
  int  ierr;
  int  line; // line where error first occured
  char dbg_functionname[STRING_LENGTH];
  char dbg_filename[STRING_LENGTH];
  char dbg_modulename[STRING_LENGTH];
  PetscBool  debug;
};
extern ebiGlobals ebi_globals;


#define DBG_FUNCT (!strcmp(ebi_globals.dbg_functionname,__FUNCT__))
#define DBG_FILE  (!strcmp(ebi_globals.dbg_filename,__FILE__))
#define DBG_MODULE (!strcmp(ebi_globals.dbg_modulename,__MODULENAME__))
#define DBG_MSG PetscPrintf( PETSC_COMM_SELF, "Code in function \"%s()\", file \"%s\", line %d\n", __FUNCT__, __FILE__,__LINE__)
#define DBG_CHK ( DBG_FUNCT || DBG_FILE || DBG_MODULE || ebi_globals.debug == PETSC_TRUE)

#ifdef EBI_DEBUG
#define eDBG( debugModeOnly_code ) if(DBG_CHK){DBG_MSG; debugModeOnly_code; }
#else
#define eDBG( debugModeOnly_code )
#endif

/* ********************************************************************** */
#define ebiMalloc(newA,A,n,obj) 0;{                                            \
                                                                              \
  if( (n)!=0){                                                                \
    ebi_globals.ierr = PetscMalloc( sizeof(A)*(n), (void**)&(newA));            \
    CHKERRQ(ebi_globals.ierr);                                                  \
                                                                              \
    newA  = (A*) newA;                                                        \
    ebi_globals.ierr = PetscMemzero( newA, sizeof(A)*(n));                      \
    CHKERRQ(ebi_globals.ierr);                                                  \
  }                                                                           \
  else { newA = PETSC_NULL; }                                                 \
  if(obj) PetscLogObjectMemory( obj, sizeof(A));                              \
}

/* ********************************************************************** */
#define ebiFree(x) 0; {                                            \
  if( (x)!=PETSC_NULL ){                                          \
    ebi_globals.ierr = PetscFree(x); CHKERRQ(ebi_globals.ierr)        \
  }                                                               \
  x = PETSC_NULL;                                                 \
}

/* ********************************************************************** */

#define ebiFreeNoReturn(x)  {                                      \
  if( (x)!=PETSC_NULL ){                                          \
    ebi_globals.ierr = PetscFree(x);                                \
  }                                                               \
  x = PETSC_NULL;                                                 \
}


/* ********************************************************************** */
/* ********************** ERROR ROUTINE ******************************* */
#include "ebi_petsc_error.hpp"


/// Users have  to call this routine
#undef __FUNCT__
#define __FUNCT__ "ebiInit"
inline ebiEC ebiInit(void)
{
  ebiFunctionBegin;
  ebiAssert(1);
  iC( PetscOptionsHasName   (NULL, PETSC_NULL,(char *) "-debug",        &ebi_globals.debug));  
  iC( PetscOptionsGetString (NULL, PETSC_NULL,(char *) "-dbg_function", ebi_globals.dbg_functionname, STRING_LENGTH, PETSC_NULL));
  iC( PetscOptionsGetString (NULL, PETSC_NULL,(char *) "-dbg_file",     ebi_globals.dbg_filename,     STRING_LENGTH, PETSC_NULL));
  iC( PetscOptionsGetString (NULL, PETSC_NULL,(char *) "-dbg_modulename",ebi_globals.dbg_modulename,  STRING_LENGTH, PETSC_NULL));
  ebiFunctionReturn(0);
}

/* ********************************************************************** */
#define ebiLogInfo(msg) 0;iC(PetscLogInfo(0, "\tEbi::%s() %d:" msg "\n",__FUNCT__,__LINE__))
#define ebiLogInfo1(msg,a) 0;iC(PetscLogInfo(0,"\tEbi::%s() %d:" msg "\n",__FUNCT__,__LINE__,a))
#define ebiLogInfo2(msg,a,b) 0;iC(PetscLogInfo(0,"\tEbi::%s() %d:" msg "\n",__FUNCT__,__LINE__,a,b))
#define ebiLogInfo3(msg,a,b,c) 0;iC(PetscLogInfo(0,"\tEbi::%s() %d:" msg "\n",__FUNCT__,__LINE__,a,b,c))
#define ebiCoutRank(rank) cout << "[" << rank << "] "

END_EBI_NAMESPACE
//#include "common/petsc/schur.h"

#endif 


