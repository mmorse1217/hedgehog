/* Various Memory Allocation routines */
#ifndef __MEMALLOC_HPP__
#define __MEMALLOC_HPP__

#include "ebi.hpp"

BEGIN_EBI_NAMESPACE

double *vecalloc(int nr);
long *lvecalloc(int nr);
char **cvecalloc(int nr);
char *charalloc(int nr);
float *fvecalloc(int nr);
int *ivecalloc(int nr);
double **matalloc(int nr, int nc);
double **matalloc(int nr, int nc);
float **fmatalloc(int nr, int nc);
int  **imatalloc(int nr, int nc);

END_EBI_NAMESPACE

#endif
