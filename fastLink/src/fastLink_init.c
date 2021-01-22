#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _fastLink_calcPWDcpp(SEXP, SEXP);
extern SEXP _fastLink_m_func_par(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_fastLink_calcPWDcpp", (DL_FUNC) &_fastLink_calcPWDcpp,  2},
    {"_fastLink_m_func_par", (DL_FUNC) &_fastLink_m_func_par, 11},
    {NULL, NULL, 0}
};

void R_init_fastLink(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
