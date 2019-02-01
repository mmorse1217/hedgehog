/*! \file */
#ifndef _EBIOBJECT_HPP_EBI_
#define _EBIOBJECT_HPP_EBI_

#include "ebi.hpp"

BEGIN_EBI_NAMESPACE

using std::cout;
using std::endl;

/// Implements distributed object
class MpiObject{
public:
  MpiObject(MPI_Comm incomm) : comm(incomm) {}
  virtual ~MpiObject() {;}
  
  const MPI_Comm& mpiComm( ) const { return  this->comm; }   ///get the communicator
  int mpiRank() const { int rank; MPI_Comm_rank( this->comm, &rank); return rank; }   ///Returns the rank (id)
  int mpiSize() const { int size; MPI_Comm_size( this->comm, &size); return size; }   ///Returns the size
protected:
  MPI_Comm comm; ///< The MPI comunicator
};

/// Implementes basic info for an object
class InfoObject{
public:
  InfoObject(const string& n, const string& p) : name_(n), prefix_(p), view_(false), user_context(NULL) {}
  virtual ~InfoObject() {}
  virtual int setFromOptions(); ///< Set object parameters from options
  ebiEC setOptionsPrefix(string &prefix_val);      ///< Set prefix for command line options
  ebiEC setOptionsPrefix(const char* prefix_val);  ///< Set prefix for command line options  
  ebiEC appendOptionsPrefix(string &prefix_val);   ///< Append a string to the prefix values
  const string & prefix()   const   { return this->prefix_;}
  bool hasSetFromOptions()  const   { return set_from_options_;}
  bool hasViewFromOptions() const   { return view_;}
  const string& name() const  { return this->name_;} ///< Get class name -- default "EbiObject"
  ebiEC setName(string input) { this->name_ = input; return 0;};///< Set the name of the class; used in profiling
  ebiEC setName(char *input)  { this->name_ = input; return 0;};///< Set the name of the class; used in profiling  
  ebiEC setUserContext(void*c){ this->user_context=c;return 0;} ///< allows users to get objects with their data 
  void *userContext() const  { return this->user_context;}  ///< get user object
protected:
  string name_;
  string prefix_; ///< prefix used to set the objects properties using line command options
  bool set_from_options_; ///< set to true if setFromOptions is called;
  bool view_;
  void *user_context; 
};
//LEXING: InfoObject(char* n, char* p) : name_(n), prefix_(p), view_(false), user_context(NULL) {}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
/// heavy  class for all objects in the EmbeddedBoundaryIntegral namespace
class EbiObject: public MpiObject, public InfoObject {
public:
  map<string, double> _function_timing;
  EbiObject(const string& name, const string& prefix): MpiObject(PETSC_COMM_WORLD), InfoObject(name, prefix) { }
  EbiObject(): MpiObject(PETSC_COMM_WORLD), InfoObject("", "") { }
  virtual ~EbiObject() {}
  virtual int view() { if(mpiRank() == 0) cout<<name_<<endl; return 0; }/// Viewing information about the object
};
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------


//EbiObject(MPI_Comm comm):MpiObject(comm), InfoObject("","") {}
//EbiObject(MPI_Comm comm, const string& name, const string& prefix): MpiObject(comm), InfoObject(name, prefix) { }
//protected:
//EbiObject& operator=( const EbiObject & rhs); 
//LEXING: EbiObject(MPI_Comm incomm, char *in_name=(char *)"EO", char *in_prefix=(char*)"eo_" ):  MpiObject(incomm), InfoObject(in_name, in_prefix) { }
//EbiObject(MPI_Comm comm, const string& name="EO_", const string& prefix="eo_"):
//  EbiObject(MPI_Comm incomm) : MpiObject(incomm), InfoObject("eo", "eo") { }
//  EbiObject(MPI_Comm incomm, char* n, char* p) : MpiObject(incomm), InfoObject(n, p) { }
/*
/// light parrent  (cannot be EbiObject and LightObject at the same time.
class LightEbiObject{
public:
  LightEbiObject () {}
  virtual ~LightEbiObject() {}
protected:
  LightEbiObject& operator=( const LightEbiObject & rhs); 
};
*/

#include "ebiobject_inline.hpp"

END_EBI_NAMESPACE

#endif /* _MOBOOBJECT_H */












