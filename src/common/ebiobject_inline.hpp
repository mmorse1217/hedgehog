

/* ********************************************************************** */
/* ********************************************************************** */
/* Implementations */


/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "InfoObject::setFromOptions"
///For derived classes: Should always include calls to the parent setFromOptions();  
inline ebiEC InfoObject::setFromOptions()
{
  ebiFunctionBegin;
  this->set_from_options_ = true;
  PetscBool flg=PETSC_FALSE;
  iC( PetscOptionsHasName(NULL, this->prefix().c_str(),"-view",&flg));
  if(flg==PETSC_TRUE) this->view_ = true;
  ebiFunctionReturn(0);
}


  
/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "InfoObject::setOptionsPrefix"
/// /param prefix_val used for the command line. Follow PETSc conventions
inline ebiEC InfoObject::setOptionsPrefix(string &prefix_val)
{
  this->prefix_ = prefix_val;
  return 0;
}

inline ebiEC InfoObject::setOptionsPrefix(const char *prefix_val)
{
  this->prefix_ = prefix_val;
  return 0;
}



/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "InfoObject::appendOptionsPrefix"
/// /param prefix_val used for the command line. Follow PETSc conventions
inline ebiEC InfoObject::appendOptionsPrefix(string &prefix_val)
{
  this->prefix_ += prefix_val;
  return 0;
}





