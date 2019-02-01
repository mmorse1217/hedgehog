/*! \file */
#ifndef _EBI_PETSC_ERROR_HPP_EBI_
#define _EBI_PETSC_ERROR_HPP_EBI_

BEGIN_EBI_NAMESPACE

#if defined(EBI_HAVE_ERROR_TRACING)

/* ********************************************************************** */
/* *******************   ERROR CHECKING ********************************* */
/* ********************************************************************** */

/* ********************************************************************** */
/// Similar to the ANCI C assert but this one traces the call stack. call ebiAssertAbrt if you want no trace
#define ebiAssert(expr) 																		\
  if(!(expr)){																						\
	 PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,1,1,"Assertion: "#expr" failed!");	\
	 ebi_globals.ierr = 1;																		\
	 goto _local_function_ebi_return_error_just_occured;								\
  }																									

#define iA(expr) { if((expr)==0) { std::cerr<<"wrong"<<std::endl; assert(expr); exit(-1); } } //assert

/// Prints a warning message in the standard PETSc error format  msg: ANSI C string
#define ebiWarning(msg) PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,1,1,(msg))

// could check if ierr is on and abort but it's ok (like a forced error checking)
/// Should be called in the begining of a function with ebiFunctionReturn(), ebiFunctionReturnVoid  at the end.
///
/// If a function uses ebiFunctionReturnValue(val), then val MUST be declared
/// \e before ebiFunctionBegin.
#define ebiFunctionBegin  ebi_globals.ierr = 0; PetscFunctionBegin;{

/*! \def  ebiAssertAbrt(expr,comm,msg)
   \brief Aborts execution to program (with MPI_Abort())

   \param expr if int(expr) !=0 the program will abort
	\param comm the MPI communicator. All processors in this communicator will abort
	\param msg  an ANSI C string to be messaged before termination
 */
#define ebiAssertAbrt(expr,comm,msg) if(!(expr)) SETERRQABORT(comm,1,msg)


/// \def ebiAssertMsg(expr,msg) 
/// \brief Like ebiAssert() but it sends msg to stderr
/// \sa eC
#define ebiAssertMsg(expr, msg)  if(!(expr)){	 ebiWarning(msg);	 ebiAssert(expr);  }


/*! \brief eC : Used for error checking

   Use after every call to an ebi funtion.
	To use it with your functions.

	

	Every function should:
	- begin with   ebiFunctionBegin;
	- if return with   ebiFunctionReturn(return_value);
	- or return with  ebiFunctionReturnVoid;
	- every call to this function should be followed by eC;
	- As with PETSc use __FUNCT__ to trace function names

	Example:
	\verbatim
	 #undef __FUNCT__
	 #define __FUNCT__ "function1"
	 double function1(args){
	   ebiFunctionBegin;
 	   .
	   . //some error
		ebiAssert(some_error==true)
		.
	   ebiFunctionReturn(some_double);
	 }

	 #undef __FUNCT__
	 #define __FUNCT__ "function1"	 
	 double function2(args){
	 	ebiFunctionBegin;
 	   .
	   . //some code
		.
		function1(some_args); eC;
		ebifunction(); eC;
      .
		.
		.
	   ebiFunctionReturn(some_double);
	 }
	 \endverbatim

	 \sa ebiAssertAbrt() ebiAssertMsg() ebiAssert() ebiWarning() ebiFunctionReturn() ebiFunctionReturnVoid
 */
#define eC {												\
  if( ebi_globals.ierr == 1) {						\
	 ebi_globals.line = __LINE__;						\
	 goto _local_function_ebi_return_error;		\
  }															\
}

/* ********************************************************************** */
/*! Needed fore error checking

    All the functions are required to use ebiFunctionReturn to terminate
	 You should start with ebiFunctionBegin. Example:
	 \verbatim
	 double function(args){
	   ebiFunctionBegin;
 	   .
	   . //some code
		.
	   ebiFunctionReturnValue(some_double);
	 }
	 \endverbatim

	 
	 \param val is the return value
	 \warning val will NOT be invalid in case an error has occured!

	 \remark Use ebiFunctionReturnVoid if you have a void function
	 \remark All calls to ebi functions should be followed by eC;
	 \sa ebiAssert() ebiAssertAbrt() ebiAssertMsg() eC 
	 
*/	 
#define ebiFunctionReturnValue(val) 												   \
  }PetscStackPop;																				\
  return val;																					\
 _local_function_ebi_return_error_just_occured:										\
  PetscStackPop;	 																			\
  return val;																					\
 _local_function_ebi_return_error:	 								   				\
  PetscError( ebi_globals.line ,__FUNCT__,__FILE__,__SDIR__,1,0," ");		\
  PetscStackPop;																				\
  return val;																					


/*! Needed for error checking */
#define ebiFunctionReturn(val) 												         \
  }PetscStackPop;																				\
  return 0;																					   \
 _local_function_ebi_return_error_just_occured:										\
  PetscStackPop;	 																			\
  return 1;																					   \
 _local_function_ebi_return_error:	 								   				\
  PetscError( ebi_globals.line ,__FUNCT__,__FILE__,__SDIR__,1,0," ");		\
  PetscStackPop;																				\
  return 1;																					   



/// This should be used for premature returns. Used with if statements.
///
/// Can be used with both error returning routines or other-type-retuning routines
#define ebiFunctionReturnEarly(val) { PetscStackPop; return val; }


/// This should be used for premature returns.
#define ebiFunctionReturnVoidEarly   {PetscStackPop; return;}


/// Must be used with void-returning functions.
/// \sa ebiReturn() eC
#define ebiFunctionReturnVoid 															\
  }PetscStackPop;                                                          \
  return;																						\
 _local_function_ebi_return_error_just_occured:										\
  PetscStackPop;	 																			\
  return;																						\
 _local_function_ebi_return_error:	 											      \
  PetscError( ebi_globals.line ,__FUNCT__,__FILE__,__SDIR__,1,0," ");		\
  PetscStackPop;																				\
  return;																						

/// iC(fun)  must be used with all PETSc objects only! ebi objects should use ebi routines
#define iC(fun) {ebi_globals.ierr=fun; eC;}

// ebiReturn* has to maintain consistency with PetscFunctionReturn() (which supports only error)

#else // ERROR TRACINGI IS NOT DEFINED

#define iA(expr) { if((expr)==0) { std::cerr<<"wrong"<<std::endl; assert(expr); exit(-1); } } //assert
#define eC
#define ebiFunctionBegin
#define ebiFunctionReturnVoid return
#define ebiFunctionReturn(val) return val
#define ebiFunctionReturnValue(val) return val
#define ebiAssert(expr) assert(expr)
#define ebiFunctionReturnEarly(val) return val
#define ebiAssertMsg(expr,msg) assert(expr)
#define ebiAssertAbrt(expr)  
#define ebiWarning(msg) cerr<<msg<<endl
#define iC(fun) fun

#endif // EBI_HAVE_ERROR_TRACING

END_EBI_NAMESPACE

#endif 


