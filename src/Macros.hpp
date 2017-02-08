#ifndef _MACROS_HPP_
#define _MACROS_HPP_

#include <iostream>
#include <cfloat>
#include <climits>
#include <omp.h>

#ifndef NDEBUG
#define _assert(ex) {if (!(ex)){fflush(stdout);fprintf(stderr,"Assertion failed: file %s, line %d\n", __FILE__, __LINE__);abort();}}
#define assert(ex) {if (!(ex)){fflush(stdout);fprintf(stderr,"Assertion failed: file %s, line %d\n", __FILE__, __LINE__);abort();}}
#else
#define _assert(ex);
#define assert(ex);
#endif

#define ALIGN_WIDTH	64
#define	uint		unsigned

#define	PRECISION	float

#endif // #ifdef _MACROS_HPP_ 
