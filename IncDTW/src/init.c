#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _IncDTW_BACKTRACK_cpp(SEXP);
extern SEXP _IncDTW_BACKTRACK2II_cpp(SEXP, SEXP);
extern SEXP _IncDTW_BACKTRACK2IN_cpp(SEXP, SEXP);
extern SEXP _IncDTW_cpp_cm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_diffm(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_cm(SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_cm_inc(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_cm_ws_ea(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_cm_ws_inc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_ea(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_inc(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_inc_mv(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_inc_mv_ws(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_inc_ws(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_mv(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_mv_ws_ea(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_v32(SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_ws(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_ws_ea(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_get_tube(SEXP, SEXP);
extern SEXP _IncDTW_cpp_get_tube_mv(SEXP, SEXP);
extern SEXP _IncDTW_cpp_kNN_rev(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_local_min(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_norm01(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw_lot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw_mv_lot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw_znorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw_znorm_lot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw_znorm_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_rundtw_znorm_mv_lot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_set_tube(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_znorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_GCM_cpp(SEXP, SEXP);
extern SEXP _IncDTW_GCM_Sakoe_cpp(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_get_lb(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_get_lb_mv1(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_get_lb_mv2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_get_lb_mv22(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_get_lb_znorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_IGCM_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_IGCM_Sakoe_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_normmat(SEXP);
extern SEXP _IncDTW_parallel_dm_dtw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_parallel_dm_dtw_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_parallel_dv_dtw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_parallel_dv_dtw_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_IncDTW_BACKTRACK_cpp",           (DL_FUNC) &_IncDTW_BACKTRACK_cpp,            1},
    {"_IncDTW_BACKTRACK2II_cpp",        (DL_FUNC) &_IncDTW_BACKTRACK2II_cpp,         2},
    {"_IncDTW_BACKTRACK2IN_cpp",        (DL_FUNC) &_IncDTW_BACKTRACK2IN_cpp,         2},
    {"_IncDTW_cpp_cm",                  (DL_FUNC) &_IncDTW_cpp_cm,                   5},
    {"_IncDTW_cpp_diffm",               (DL_FUNC) &_IncDTW_cpp_diffm,                4},
    {"_IncDTW_cpp_dtw2vec",             (DL_FUNC) &_IncDTW_cpp_dtw2vec,              3},
    {"_IncDTW_cpp_dtw2vec_cm",          (DL_FUNC) &_IncDTW_cpp_dtw2vec_cm,           2},
    {"_IncDTW_cpp_dtw2vec_cm_inc",      (DL_FUNC) &_IncDTW_cpp_dtw2vec_cm_inc,       3},
    {"_IncDTW_cpp_dtw2vec_cm_ws_ea",    (DL_FUNC) &_IncDTW_cpp_dtw2vec_cm_ws_ea,     4},
    {"_IncDTW_cpp_dtw2vec_cm_ws_inc",   (DL_FUNC) &_IncDTW_cpp_dtw2vec_cm_ws_inc,    5},
    {"_IncDTW_cpp_dtw2vec_ea",          (DL_FUNC) &_IncDTW_cpp_dtw2vec_ea,           4},
    {"_IncDTW_cpp_dtw2vec_inc",         (DL_FUNC) &_IncDTW_cpp_dtw2vec_inc,          4},
    {"_IncDTW_cpp_dtw2vec_inc_mv",      (DL_FUNC) &_IncDTW_cpp_dtw2vec_inc_mv,       5},
    {"_IncDTW_cpp_dtw2vec_inc_mv_ws",   (DL_FUNC) &_IncDTW_cpp_dtw2vec_inc_mv_ws,    7},
    {"_IncDTW_cpp_dtw2vec_inc_ws",      (DL_FUNC) &_IncDTW_cpp_dtw2vec_inc_ws,       6},
    {"_IncDTW_cpp_dtw2vec_mv",          (DL_FUNC) &_IncDTW_cpp_dtw2vec_mv,           4},
    {"_IncDTW_cpp_dtw2vec_mv_ws_ea",    (DL_FUNC) &_IncDTW_cpp_dtw2vec_mv_ws_ea,     6},
    {"_IncDTW_cpp_dtw2vec_v32",         (DL_FUNC) &_IncDTW_cpp_dtw2vec_v32,          2},
    {"_IncDTW_cpp_dtw2vec_ws",          (DL_FUNC) &_IncDTW_cpp_dtw2vec_ws,           4},
    {"_IncDTW_cpp_dtw2vec_ws_ea",       (DL_FUNC) &_IncDTW_cpp_dtw2vec_ws_ea,        5},
    {"_IncDTW_cpp_get_tube",            (DL_FUNC) &_IncDTW_cpp_get_tube,             2},
    {"_IncDTW_cpp_get_tube_mv",         (DL_FUNC) &_IncDTW_cpp_get_tube_mv,          2},
    {"_IncDTW_cpp_kNN_rev",             (DL_FUNC) &_IncDTW_cpp_kNN_rev,              3},
    {"_IncDTW_cpp_local_min",           (DL_FUNC) &_IncDTW_cpp_local_min,            3},
    {"_IncDTW_cpp_norm01",              (DL_FUNC) &_IncDTW_cpp_norm01,               4},
    {"_IncDTW_cpp_rundtw",              (DL_FUNC) &_IncDTW_cpp_rundtw,              11},
    {"_IncDTW_cpp_rundtw_lot",          (DL_FUNC) &_IncDTW_cpp_rundtw_lot,          15},
    {"_IncDTW_cpp_rundtw_mv",           (DL_FUNC) &_IncDTW_cpp_rundtw_mv,           12},
    {"_IncDTW_cpp_rundtw_mv_lot",       (DL_FUNC) &_IncDTW_cpp_rundtw_mv_lot,       16},
    {"_IncDTW_cpp_rundtw_znorm",        (DL_FUNC) &_IncDTW_cpp_rundtw_znorm,        10},
    {"_IncDTW_cpp_rundtw_znorm_lot",    (DL_FUNC) &_IncDTW_cpp_rundtw_znorm_lot,    14},
    {"_IncDTW_cpp_rundtw_znorm_mv",     (DL_FUNC) &_IncDTW_cpp_rundtw_znorm_mv,     11},
    {"_IncDTW_cpp_rundtw_znorm_mv_lot", (DL_FUNC) &_IncDTW_cpp_rundtw_znorm_mv_lot, 15},
    {"_IncDTW_cpp_set_tube",            (DL_FUNC) &_IncDTW_cpp_set_tube,             3},
    {"_IncDTW_cpp_znorm",               (DL_FUNC) &_IncDTW_cpp_znorm,                4},
    {"_IncDTW_GCM_cpp",                 (DL_FUNC) &_IncDTW_GCM_cpp,                  2},
    {"_IncDTW_GCM_Sakoe_cpp",           (DL_FUNC) &_IncDTW_GCM_Sakoe_cpp,            3},
    {"_IncDTW_get_lb",                  (DL_FUNC) &_IncDTW_get_lb,                   4},
    {"_IncDTW_get_lb_mv1",              (DL_FUNC) &_IncDTW_get_lb_mv1,               5},
    {"_IncDTW_get_lb_mv2",              (DL_FUNC) &_IncDTW_get_lb_mv2,               5},
    {"_IncDTW_get_lb_mv22",             (DL_FUNC) &_IncDTW_get_lb_mv22,              5},
    {"_IncDTW_get_lb_znorm",            (DL_FUNC) &_IncDTW_get_lb_znorm,             7},
    {"_IncDTW_IGCM_cpp",                (DL_FUNC) &_IncDTW_IGCM_cpp,                 4},
    {"_IncDTW_IGCM_Sakoe_cpp",          (DL_FUNC) &_IncDTW_IGCM_Sakoe_cpp,           5},
    {"_IncDTW_normmat",                 (DL_FUNC) &_IncDTW_normmat,                  1},
    {"_IncDTW_parallel_dm_dtw",         (DL_FUNC) &_IncDTW_parallel_dm_dtw,          7},
    {"_IncDTW_parallel_dm_dtw_mv",      (DL_FUNC) &_IncDTW_parallel_dm_dtw_mv,       8},
    {"_IncDTW_parallel_dv_dtw",         (DL_FUNC) &_IncDTW_parallel_dv_dtw,          6},
    {"_IncDTW_parallel_dv_dtw_mv",      (DL_FUNC) &_IncDTW_parallel_dv_dtw_mv,       7},
    {NULL, NULL, 0}
};

void R_init_IncDTW(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
