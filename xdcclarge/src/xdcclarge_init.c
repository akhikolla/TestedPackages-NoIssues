#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _xdcclarge_cdcc_compositelik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xdcclarge_cdcc_construct(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xdcclarge_dcc_compositelik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _xdcclarge_dcc_construct(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_xdcclarge_cdcc_compositelik", (DL_FUNC) &_xdcclarge_cdcc_compositelik, 8},
  {"_xdcclarge_cdcc_construct",    (DL_FUNC) &_xdcclarge_cdcc_construct,    7},
  {"_xdcclarge_dcc_compositelik",  (DL_FUNC) &_xdcclarge_dcc_compositelik,  8},
  {"_xdcclarge_dcc_construct",     (DL_FUNC) &_xdcclarge_dcc_construct,     7},
  {NULL, NULL, 0}
};

void R_init_xdcclarge(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
