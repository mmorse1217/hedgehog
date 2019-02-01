/* mem.c   LAST MODIFIED:  10/07/98       */
/* Memory allocation subroutines from the */
/* test Numerical Computation Using C by  */
/* R. Glassey and elsewhere...            */
#include <stdio.h>
#include <stdlib.h>
#include "memAlloc.hpp"

BEGIN_EBI_NAMESPACE

double *vecalloc(int nr)
{
  double *x;
  x=(double *) calloc((unsigned) nr,sizeof(double));
  if (x==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
  return x;
}

long *lvecalloc(int nr)
{
  long *x;
  x=(long *) calloc((unsigned) nr,sizeof(long));
  if (x==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
  return x;
}

/* See page 50 of Programs and Data Structures in C by Ammeraal */
/* for example from which the two following routines where drawn. */
char **cvecalloc(int nr)
{
  char **x;
  x=(char **) malloc(nr * sizeof(char *));
  if (x==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
  return x;
}

char *charalloc(int nr)
{
  char *x;
  x=(char *) malloc(nr);
  if (x==NULL){
    fprintf(stderr, "Unable to allocate dynamic memory");
    exit(-1);
  }
  return x;
}


float *fvecalloc(int nr)
{
  float *x;
  x=(float *) calloc((unsigned) nr,sizeof(float));
  if (x==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
  return x;
}

int *ivecalloc(int nr)
{
  int *x;
  x=(int *) calloc((unsigned) nr,sizeof(int));
  if (x==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
  return x;
}

double **matalloc(int nr, int nc)
{
int k;
double **x;
x=(double **) calloc((unsigned) nr, sizeof(double *));
if (x==NULL){
  fprintf(stderr,"Unable to allocate dynamic memory");
  exit(-1);
}
for (k=0;k<nr;k++){
  x[k]=(double *) calloc((unsigned) nc, sizeof(double));
  if (x[k]==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
}
return x;
}

float **fmatalloc(int nr, int nc)
{
int k;
float **x;
x=(float **) calloc((unsigned) nr, sizeof(float *));
if (x==NULL){
  fprintf(stderr,"Unable to allocate dynamic memory");
  exit(-1);
}
for (k=0;k<nr;k++){
  x[k]=(float *) calloc((unsigned) nc, sizeof(float));
  if (x[k]==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
}
return x;
}


int  **imatalloc(int nr, int nc)
{
int k;
int **x;
x=(int **) calloc((unsigned) nr, sizeof(int *));
if (x==NULL){
  fprintf(stderr,"Unable to allocate dynamic memory");
  exit(-1);
}
for (k=0;k<nr;k++){
  x[k]=(int *) calloc((unsigned) nc, sizeof(int));
  if (x[k]==NULL){
    fprintf(stderr,"Unable to allocate dynamic memory");
    exit(-1);
  }
}
return x;
}

END_EBI_NAMESPACE
