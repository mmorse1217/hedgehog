#ifndef _EBI_INCLUDE_HPP_
#define _EBI_INCLUDE_HPP_

#include "ebi_namespace.hpp"

#ifdef __DECCXX
#define __USE_STD_IOSTREAM
#endif

#if defined(__GNUG__) || defined(_STANDARD_C_PLUS_PLUS) || defined(__DECCXX) || defined(sun) || defined(__linux)
#include <iostream>
#include <fstream>
#include <sstream>
#else
//#include <iostream.h>
//#include <fstream.h>
//#include <sstream.h>
#endif


#include <cfloat>
#include <cassert>
#include <cmath>
#include <string>
#include <complex>

#include <vector>
#include <set>
#include <map>
#include <deque>
#include <queue>
#include <utility>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

//FFTW stuff
#ifdef FFTW3
#include "fftw3.h"
#else
//#include "rfftw.h"
#endif

#endif
