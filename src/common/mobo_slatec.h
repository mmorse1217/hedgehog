#ifndef _MOBO_SLATEC_H_
#define _MOBO_SLATEC_H_

#include "ebi_namespace.hpp"

BEGIN_EBI_NAMESPACE

#define DBESK0 dbesk0_
#define DBESK1 dbesk1_

extern "C"
{
  double DBESK0(double*);
  double DBESK1(double*);
}

END_EBI_NAMESPACE

#endif
