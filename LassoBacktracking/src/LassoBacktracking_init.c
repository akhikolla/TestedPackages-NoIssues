#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
/* .Call calls */
extern SEXP LassoBacktracking_add_inter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LassoBacktracking_add_inter_orig(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LassoBacktracking_any_ind(SEXP, SEXP);
extern SEXP LassoBacktracking_any_indmax(SEXP, SEXP);
extern SEXP LassoBacktracking_beta_active(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LassoBacktracking_change_dim(SEXP, SEXP);
extern SEXP LassoBacktracking_find_l0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LassoBacktracking_in_log(SEXP, SEXP);
extern SEXP LassoBacktracking_inner_prod_abs_comp(SEXP, SEXP, SEXP, SEXP);
extern SEXP LassoBacktracking_inner_prod_abs_comp2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LassoBacktracking_scale_cen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LassoBacktracking_which_indmax(SEXP, SEXP);
extern SEXP LassoBacktracking_zero(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"LassoBacktracking_add_inter",            (DL_FUNC) &LassoBacktracking_add_inter,             9},
  {"LassoBacktracking_add_inter_orig",       (DL_FUNC) &LassoBacktracking_add_inter_orig,        5},
  {"LassoBacktracking_any_ind",              (DL_FUNC) &LassoBacktracking_any_ind,               2},
  {"LassoBacktracking_any_indmax",           (DL_FUNC) &LassoBacktracking_any_indmax,            2},
  {"LassoBacktracking_beta_active",          (DL_FUNC) &LassoBacktracking_beta_active,          12},
  {"LassoBacktracking_change_dim",           (DL_FUNC) &LassoBacktracking_change_dim,            2},
  {"LassoBacktracking_find_l0",              (DL_FUNC) &LassoBacktracking_find_l0,               6},
  {"LassoBacktracking_in_log",               (DL_FUNC) &LassoBacktracking_in_log,                2},
  {"LassoBacktracking_inner_prod_abs_comp",  (DL_FUNC) &LassoBacktracking_inner_prod_abs_comp,   4},
  {"LassoBacktracking_inner_prod_abs_comp2", (DL_FUNC) &LassoBacktracking_inner_prod_abs_comp2,  5},
  {"LassoBacktracking_scale_cen",            (DL_FUNC) &LassoBacktracking_scale_cen,             6},
  {"LassoBacktracking_which_indmax",         (DL_FUNC) &LassoBacktracking_which_indmax,          2},
  {"LassoBacktracking_zero",                 (DL_FUNC) &LassoBacktracking_zero,                  2},
  {NULL, NULL, 0}
};

void R_init_LassoBacktracking(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
