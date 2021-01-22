#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ICcalib_Calcb(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CalcbZ(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CalcNablabeetaUbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CalcNablabeetaUgamma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CalcSurvFromNPMLE(SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CalcUbetabeeta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CalcUbetabeetaRS(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CoxLogLik(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CoxLogLikGrad(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CoxLogLikHess(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CoxLogLikNoBeta(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_CoxLogLikX(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICcalib_FindIntervalCalibCPP(SEXP, SEXP);
extern SEXP _ICcalib_FindIntervalCalibCPPvec(SEXP, SEXP);
extern SEXP _ICcalib_FindIntervalCPP(SEXP, SEXP);
extern SEXP _ICcalib_myFmyHess(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_ICcalib_Calcb",                   (DL_FUNC) &_ICcalib_Calcb,                   4},
  {"_ICcalib_CalcbZ",                  (DL_FUNC) &_ICcalib_CalcbZ,                  5},
  {"_ICcalib_CalcNablabeetaUbeta",     (DL_FUNC) &_ICcalib_CalcNablabeetaUbeta,     6},
  {"_ICcalib_CalcNablabeetaUgamma",    (DL_FUNC) &_ICcalib_CalcNablabeetaUgamma,    6},
  {"_ICcalib_CalcSurvFromNPMLE",       (DL_FUNC) &_ICcalib_CalcSurvFromNPMLE,       3},
  {"_ICcalib_CalcUbetabeeta",          (DL_FUNC) &_ICcalib_CalcUbetabeeta,          5},
  {"_ICcalib_CalcUbetabeetaRS",        (DL_FUNC) &_ICcalib_CalcUbetabeetaRS,        5},
  {"_ICcalib_CoxLogLik",               (DL_FUNC) &_ICcalib_CoxLogLik,               5},
  {"_ICcalib_CoxLogLikGrad",           (DL_FUNC) &_ICcalib_CoxLogLikGrad,           5},
  {"_ICcalib_CoxLogLikHess",           (DL_FUNC) &_ICcalib_CoxLogLikHess,           5},
  {"_ICcalib_CoxLogLikNoBeta",         (DL_FUNC) &_ICcalib_CoxLogLikNoBeta,         4},
  {"_ICcalib_CoxLogLikX",              (DL_FUNC) &_ICcalib_CoxLogLikX,              4},
  {"_ICcalib_FindIntervalCalibCPP",    (DL_FUNC) &_ICcalib_FindIntervalCalibCPP,    2},
  {"_ICcalib_FindIntervalCalibCPPvec", (DL_FUNC) &_ICcalib_FindIntervalCalibCPPvec, 2},
  {"_ICcalib_FindIntervalCPP",         (DL_FUNC) &_ICcalib_FindIntervalCPP,         2},
  {"_ICcalib_myFmyHess",               (DL_FUNC) &_ICcalib_myFmyHess,               4},
  {NULL, NULL, 0}
};

void R_init_ICcalib(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
