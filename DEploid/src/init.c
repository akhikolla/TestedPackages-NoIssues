#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _DEploid_dEploid(SEXP);
extern SEXP _DEploid_extractVcf(SEXP);
extern SEXP _DEploid_test_RRG_sample();
extern SEXP _DEploid_test_RRG_sampleExpoExpoLimit(SEXP, SEXP, SEXP);
extern SEXP _DEploid_test_RRG_sampleUnitExpo();

static const R_CallMethodDef CallEntries[] = {
    {"_DEploid_dEploid",                      (DL_FUNC) &_DEploid_dEploid,                      1},
    {"_DEploid_extractVcf",                   (DL_FUNC) &_DEploid_extractVcf,                   1},
    {"_DEploid_test_RRG_sample",              (DL_FUNC) &_DEploid_test_RRG_sample,              0},
    {"_DEploid_test_RRG_sampleExpoExpoLimit", (DL_FUNC) &_DEploid_test_RRG_sampleExpoExpoLimit, 3},
    {"_DEploid_test_RRG_sampleUnitExpo",      (DL_FUNC) &_DEploid_test_RRG_sampleUnitExpo,      0},
    {NULL, NULL, 0}
};

void R_init_DEploid(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
