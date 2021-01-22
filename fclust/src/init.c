#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _fclust_centroids_FKM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_centroids_FKM_ent(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_centroids_FKM_pf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_distCheck(SEXP, SEXP, SEXP);
extern SEXP _fclust_euclidean_distance(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_euclidean_distance_gk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_euclidean_distance_gkb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_euclidean_distance_medoid(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_F_gk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_F_gk_ent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_F_gkb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_F_gkb_ent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_indices(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_InvCheck(SEXP);
extern SEXP _fclust_mainFKM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_ent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_ent_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_ent_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_ent_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk_ent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk_ent_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk_ent_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk_ent_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gk_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb_ent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb_ent_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb_ent_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb_ent_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_gkb_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_med(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_med_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_med_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_med_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_pf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_pf_noise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_pf_noise_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_pf_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainFKM_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainnefrc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainnefrc_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainrnefrc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_mainrnefrc_U(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_Match(SEXP, SEXP);
extern SEXP _fclust_medoids_FKM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_memb_degree(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_memb_degree_ent(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_memb_degree_medoid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_memb_degree_pf(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_partCoef(SEXP, SEXP);
extern SEXP _fclust_partCoef_mod(SEXP, SEXP, SEXP);
extern SEXP _fclust_partEntropy(SEXP, SEXP, SEXP);
extern SEXP _fclust_partition_comp(SEXP, SEXP, SEXP);
extern SEXP _fclust_replace(SEXP, SEXP, SEXP);
extern SEXP _fclust_silhouette(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_silhouette_internal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_silhouetteFuzzy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fclust_unifInit(SEXP, SEXP);
extern SEXP _fclust_xie_beni(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_fclust_centroids_FKM",             (DL_FUNC) &_fclust_centroids_FKM,              6},
    {"_fclust_centroids_FKM_ent",         (DL_FUNC) &_fclust_centroids_FKM_ent,          5},
    {"_fclust_centroids_FKM_pf",          (DL_FUNC) &_fclust_centroids_FKM_pf,           6},
    {"_fclust_distCheck",                 (DL_FUNC) &_fclust_distCheck,                  3},
    {"_fclust_euclidean_distance",        (DL_FUNC) &_fclust_euclidean_distance,         5},
    {"_fclust_euclidean_distance_gk",     (DL_FUNC) &_fclust_euclidean_distance_gk,      8},
    {"_fclust_euclidean_distance_gkb",    (DL_FUNC) &_fclust_euclidean_distance_gkb,     6},
    {"_fclust_euclidean_distance_medoid", (DL_FUNC) &_fclust_euclidean_distance_medoid,  5},
    {"_fclust_F_gk",                      (DL_FUNC) &_fclust_F_gk,                       8},
    {"_fclust_F_gk_ent",                  (DL_FUNC) &_fclust_F_gk_ent,                   7},
    {"_fclust_F_gkb",                     (DL_FUNC) &_fclust_F_gkb,                     11},
    {"_fclust_F_gkb_ent",                 (DL_FUNC) &_fclust_F_gkb_ent,                 10},
    {"_fclust_indices",                   (DL_FUNC) &_fclust_indices,                   11},
    {"_fclust_InvCheck",                  (DL_FUNC) &_fclust_InvCheck,                   1},
    {"_fclust_mainFKM",                   (DL_FUNC) &_fclust_mainFKM,                   10},
    {"_fclust_mainFKM_ent",               (DL_FUNC) &_fclust_mainFKM_ent,               10},
    {"_fclust_mainFKM_ent_noise",         (DL_FUNC) &_fclust_mainFKM_ent_noise,         11},
    {"_fclust_mainFKM_ent_noise_U",       (DL_FUNC) &_fclust_mainFKM_ent_noise_U,       11},
    {"_fclust_mainFKM_ent_U",             (DL_FUNC) &_fclust_mainFKM_ent_U,             10},
    {"_fclust_mainFKM_gk",                (DL_FUNC) &_fclust_mainFKM_gk,                11},
    {"_fclust_mainFKM_gk_ent",            (DL_FUNC) &_fclust_mainFKM_gk_ent,            11},
    {"_fclust_mainFKM_gk_ent_noise",      (DL_FUNC) &_fclust_mainFKM_gk_ent_noise,      12},
    {"_fclust_mainFKM_gk_ent_noise_U",    (DL_FUNC) &_fclust_mainFKM_gk_ent_noise_U,    12},
    {"_fclust_mainFKM_gk_ent_U",          (DL_FUNC) &_fclust_mainFKM_gk_ent_U,          11},
    {"_fclust_mainFKM_gk_noise",          (DL_FUNC) &_fclust_mainFKM_gk_noise,          12},
    {"_fclust_mainFKM_gk_noise_U",        (DL_FUNC) &_fclust_mainFKM_gk_noise_U,        12},
    {"_fclust_mainFKM_gk_U",              (DL_FUNC) &_fclust_mainFKM_gk_U,              11},
    {"_fclust_mainFKM_gkb",               (DL_FUNC) &_fclust_mainFKM_gkb,               13},
    {"_fclust_mainFKM_gkb_ent",           (DL_FUNC) &_fclust_mainFKM_gkb_ent,           13},
    {"_fclust_mainFKM_gkb_ent_noise",     (DL_FUNC) &_fclust_mainFKM_gkb_ent_noise,     14},
    {"_fclust_mainFKM_gkb_ent_noise_U",   (DL_FUNC) &_fclust_mainFKM_gkb_ent_noise_U,   14},
    {"_fclust_mainFKM_gkb_ent_U",         (DL_FUNC) &_fclust_mainFKM_gkb_ent_U,         13},
    {"_fclust_mainFKM_gkb_noise",         (DL_FUNC) &_fclust_mainFKM_gkb_noise,         14},
    {"_fclust_mainFKM_gkb_noise_U",       (DL_FUNC) &_fclust_mainFKM_gkb_noise_U,       14},
    {"_fclust_mainFKM_gkb_U",             (DL_FUNC) &_fclust_mainFKM_gkb_U,             13},
    {"_fclust_mainFKM_med",               (DL_FUNC) &_fclust_mainFKM_med,               10},
    {"_fclust_mainFKM_med_noise",         (DL_FUNC) &_fclust_mainFKM_med_noise,         11},
    {"_fclust_mainFKM_med_noise_U",       (DL_FUNC) &_fclust_mainFKM_med_noise_U,       11},
    {"_fclust_mainFKM_med_U",             (DL_FUNC) &_fclust_mainFKM_med_U,             10},
    {"_fclust_mainFKM_noise",             (DL_FUNC) &_fclust_mainFKM_noise,             11},
    {"_fclust_mainFKM_noise_U",           (DL_FUNC) &_fclust_mainFKM_noise_U,           11},
    {"_fclust_mainFKM_pf",                (DL_FUNC) &_fclust_mainFKM_pf,                10},
    {"_fclust_mainFKM_pf_noise",          (DL_FUNC) &_fclust_mainFKM_pf_noise,          11},
    {"_fclust_mainFKM_pf_noise_U",        (DL_FUNC) &_fclust_mainFKM_pf_noise_U,        11},
    {"_fclust_mainFKM_pf_U",              (DL_FUNC) &_fclust_mainFKM_pf_U,              10},
    {"_fclust_mainFKM_U",                 (DL_FUNC) &_fclust_mainFKM_U,                 10},
    {"_fclust_mainnefrc",                 (DL_FUNC) &_fclust_mainnefrc,                  9},
    {"_fclust_mainnefrc_U",               (DL_FUNC) &_fclust_mainnefrc_U,                9},
    {"_fclust_mainrnefrc",                (DL_FUNC) &_fclust_mainrnefrc,                10},
    {"_fclust_mainrnefrc_U",              (DL_FUNC) &_fclust_mainrnefrc_U,              10},
    {"_fclust_Match",                     (DL_FUNC) &_fclust_Match,                      2},
    {"_fclust_medoids_FKM",               (DL_FUNC) &_fclust_medoids_FKM,                6},
    {"_fclust_memb_degree",               (DL_FUNC) &_fclust_memb_degree,                5},
    {"_fclust_memb_degree_ent",           (DL_FUNC) &_fclust_memb_degree_ent,            5},
    {"_fclust_memb_degree_medoid",        (DL_FUNC) &_fclust_memb_degree_medoid,         6},
    {"_fclust_memb_degree_pf",            (DL_FUNC) &_fclust_memb_degree_pf,             5},
    {"_fclust_partCoef",                  (DL_FUNC) &_fclust_partCoef,                   2},
    {"_fclust_partCoef_mod",              (DL_FUNC) &_fclust_partCoef_mod,               3},
    {"_fclust_partEntropy",               (DL_FUNC) &_fclust_partEntropy,                3},
    {"_fclust_partition_comp",            (DL_FUNC) &_fclust_partition_comp,             3},
    {"_fclust_replace",                   (DL_FUNC) &_fclust_replace,                    3},
    {"_fclust_silhouette",                (DL_FUNC) &_fclust_silhouette,                 6},
    {"_fclust_silhouette_internal",       (DL_FUNC) &_fclust_silhouette_internal,        6},
    {"_fclust_silhouetteFuzzy",           (DL_FUNC) &_fclust_silhouetteFuzzy,            7},
    {"_fclust_unifInit",                  (DL_FUNC) &_fclust_unifInit,                   2},
    {"_fclust_xie_beni",                  (DL_FUNC) &_fclust_xie_beni,                   6},
    {NULL, NULL, 0}
};

void R_init_fclust(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
