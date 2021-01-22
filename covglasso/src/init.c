#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.

 Routines registration obtained with
 tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
 */

/* .Call calls */
extern SEXP _covglasso_covglasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _covglasso_covglassopath_bic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _covglasso_profileloglik(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_covglasso_covglasso",         (DL_FUNC) &_covglasso_covglasso,         8},
    {"_covglasso_covglassopath_bic", (DL_FUNC) &_covglasso_covglassopath_bic, 9},
    {"_covglasso_profileloglik",     (DL_FUNC) &_covglasso_profileloglik,     3},
    {NULL, NULL, 0}
};

void R_init_covglasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
