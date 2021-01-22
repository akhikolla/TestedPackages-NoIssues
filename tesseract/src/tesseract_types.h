#include <tesseract/baseapi.h>
#include <leptonica/allheaders.h>

#define R_NO_REMAP
#define STRICT_R_HEADERS

#include <Rcpp.h>

inline void tess_finalizer(tesseract::TessBaseAPI *engine) {
  engine->End();
  delete engine;
}

//double check C++ API version for Xptr finalizer
#if RCPP_VERSION < Rcpp_Version(0,12,10)
#error RCPP too old. Need at least 0.12.10.
#endif

typedef Rcpp::XPtr<tesseract::TessBaseAPI, Rcpp::PreserveStorage, tess_finalizer, true> TessPtr;
