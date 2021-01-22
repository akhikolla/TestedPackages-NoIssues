// All source files should include this single header file. The idea here
// is to make it easy to switch to using precompiled headers. Currently,
// the full compile is fast enough that I haven't bothered.

#ifdef EIGEN_WORLD_VERSION
#error "rpf.h must be included before Eigen"
#endif

#define EIGEN_NO_DEBUG 1
#define EIGEN_DONT_PARALLELIZE

#ifdef DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_NAN
#undef EIGEN_NO_DEBUG
#endif

#include <Rcpp.h>
using namespace Rcpp;

#include "Eigen/Core"

static inline int triangleLoc1(int diag)
{
	return (diag) * (diag+1) / 2;   // 0 1 3 6 10 15 ..
}

static inline int triangleLoc0(int diag)
{
	return triangleLoc1(diag+1) - 1;  // 0 2 5 9 14 ..
}

#include "../inst/include/libifa-rpf.h"
#include "dmvnorm.h"
#include "ba81quad.h"

extern int GlobalNumberOfCores;

static inline bool strEQ(const char *s1, const char *s2) { return strcmp(s1,s2)==0; }

static inline void
pda(const double *ar, int rows, int cols) {   // column major order
	for (int rx=0; rx < rows; rx++) {
		for (int cx=0; cx < cols; cx++) {
			Rprintf("%.6g, ", ar[cx * rows + rx]);
		}
		Rprintf("\n");
	}
}
