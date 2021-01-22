#include "rpf.h"

// [[Rcpp::export]]
int get_model_id(const StringVector &str)
{
	const char *target = str[0];
  for (int sx=0; sx < Glibrpf_numModels; sx++) {
    if (strcmp(Glibrpf_model[sx].name, target) == 0) {
			return sx;
    }
  }
  return NA_REAL;
}

static int getSpecID(const NumericVector &spec)
{
  if (spec.size() < RPF_ISpecCount)
    stop("Item spec must be of length %d, not %d",
				 RPF_ISpecCount, int(spec.size()));

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    stop("Item model %d out of range", id);

	return id;
}

// [[Rcpp::export]]
int numSpec(const NumericVector &spec)
{
	int id = getSpecID(spec);
  int numSpec = (*Glibrpf_model[id].numSpec)(spec.begin());
	return numSpec;
}

// [[Rcpp::export]]
int numParam(const NumericVector &spec)
{
	int id = getSpecID(spec);
  int numParam = (*Glibrpf_model[id].numParam)(spec.begin());
	return numParam;
}

// [[Rcpp::export]]
SEXP paramInfo(const NumericVector &spec, int pnum)
{
	int id = getSpecID(spec);

  int numParam = (*Glibrpf_model[id].numParam)(spec.begin());
  if (pnum < 0 || pnum >= numParam) stop("Item model %d has %d parameters", id, numParam);

  const char *type;
  double upper, lower;
  (*Glibrpf_model[id].paramInfo)(spec.begin(), pnum, &type, &upper, &lower);

  int len = 3;
  SEXP names, ans;
  Rf_protect(names = Rf_allocVector(STRSXP, len));
  Rf_protect(ans = Rf_allocVector(VECSXP, len));
  int lx = 0;
  SET_STRING_ELT(names, lx, Rf_mkChar("type"));
  SET_VECTOR_ELT(ans,   lx, Rf_ScalarString(Rf_mkChar(type)));
  SET_STRING_ELT(names, ++lx, Rf_mkChar("upper"));
  SET_VECTOR_ELT(ans,   lx, Rf_ScalarReal(std::isfinite(upper)? upper : NA_REAL));
  SET_STRING_ELT(names, ++lx, Rf_mkChar("lower"));
  SET_VECTOR_ELT(ans,   lx, Rf_ScalarReal(std::isfinite(lower)? lower : NA_REAL));
  Rf_namesgets(ans, names);
  UNPROTECT(2);

  return ans;
}

int unpack_theta(int dims, double *param, int numAbilities, double *theta, double *out)
{
  if (numAbilities == dims) {
    for (int dx=0; dx < dims; ++dx) {
      double th = theta[dx];
      if (!std::isfinite(th)) return 0;
      out[dx] = th;
    }
  } else {
    int ax = 0;
    for (int dx=0; dx < dims; ++dx) {
      if (param[dx] == 0) continue;
      double th = theta[ax]; // could read uninitialized memory, but we detect below
      if (!std::isfinite(th)) return 0;
      out[dx] = th;
      ++ax;
    }
    if (ax != numAbilities) {
      stop("Item has %d nonzero dims but given %d abilities", ax, numAbilities);
    }
  }
  return 1;
}

// [[Rcpp::export]]
NumericMatrix prob(const NumericVector &spec, SEXP r_param, RObject r_theta)
{
	int id = getSpecID(spec);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec.begin());
  if (spec.size() < numSpec)
    stop("Item spec must be of length %d, not %d", numSpec, spec.size());
    
  int numParam = (*Glibrpf_model[id].numParam)(spec.begin());
  if (Rf_length(r_param) < numParam)
    stop("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const int dims = spec[RPF_ISpecDims];
  double *param = REAL(r_param);
  int numPeople = 1;
  int numAbilities = 1;
  double *theta = 0;
  if (dims == 0) {
	  if (r_theta != R_NilValue) {
			NumericVector tvec = as<NumericVector>(r_theta);
			numPeople = tvec.size();
		}
  } else if (dims == 1) {
		NumericVector tvec = as<NumericVector>(r_theta);
    numPeople = Rf_length(tvec);
		theta = tvec.begin();
  } else {
		NumericMatrix tmat = as<NumericMatrix>(r_theta);
		numAbilities = tmat.nrow();
		numPeople = tmat.ncol();
		theta = tmat.begin();
  }

	NumericMatrix out(numOutcomes, numPeople);
  Eigen::VectorXd thBuf(dims);
  for (int px=0; px < numPeople; px++) {
	  if (dims && !unpack_theta(dims, param, numAbilities, theta + px*numAbilities, thBuf.data())) {
		  for (int ox=0; ox < numOutcomes; ox++) {
			  out(ox, px) = NA_REAL;
		  }
		  continue;
	  }
	  (*Glibrpf_model[id].prob)(spec.begin(), param, thBuf.data(), out.begin()+px*numOutcomes);
    for (int ox=0; ox < numOutcomes; ox++) {
      double prob = out(ox, px);
      if (!std::isfinite(prob)) {
				out(ox, px) = NA_REAL;  // legitimate (e.g., grm thresholds misordered)
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
SEXP logprob(const NumericVector &spec, SEXP r_param, RObject r_theta)
{
	int id = getSpecID(spec);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec.begin());
  if (spec.size() < numSpec)
    stop("Item spec must be of length %d, not %d", numSpec, spec.size());
    
  int numParam = (*Glibrpf_model[id].numParam)(spec.begin());
  if (Rf_length(r_param) < numParam)
    stop("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const int dims = spec[RPF_ISpecDims];
  int numPeople = 1;
  double *param = REAL(r_param);
  int numAbilities = 1;
  double *theta = 0;
  if (dims == 0) {
	  if (r_theta != R_NilValue) {
			NumericVector tvec = as<NumericVector>(r_theta);
			numPeople = tvec.size();
		}
  } else if (dims == 1) {
		NumericVector tvec = as<NumericVector>(r_theta);
    numPeople = Rf_length(tvec);
		theta = tvec.begin();
  } else {
		NumericMatrix tmat = as<NumericMatrix>(r_theta);
		numAbilities = tmat.nrow();
		numPeople = tmat.ncol();
		theta = tmat.begin();
  }

	NumericMatrix out(numOutcomes, numPeople);
  Eigen::VectorXd thBuf(dims);
  for (int px=0; px < numPeople; px++) {
	  if (dims && !unpack_theta(dims, param, numAbilities, theta + px*numAbilities, thBuf.data())) {
		  for (int ox=0; ox < numOutcomes; ox++) {
			  out(ox, px) = NA_REAL;
		  }
		  continue;
	  }
	  (*Glibrpf_model[id].logprob)(spec.begin(), param, thBuf.data(), out.begin()+px*numOutcomes);
    for (int ox=0; ox < numOutcomes; ox++) {
      double prob = out(ox, px);
      if (!std::isfinite(prob)) {
				out(ox, px) = NA_REAL;  // legitimate (e.g., grm thresholds misordered)
      }
    }
  }
	return out;
}

// [[Rcpp::export]]
SEXP dLL(const NumericVector &spec, SEXP r_param, SEXP r_where, SEXP r_weight)
{
	int id = getSpecID(spec);
  int numSpec = (*Glibrpf_model[id].numSpec)(spec.begin());
  if (spec.size() < numSpec)
    stop("Item spec must be of length %d, not %d", numSpec, spec.size());
    
  int numParam = (*Glibrpf_model[id].numParam)(spec.begin());
  if (Rf_length(r_param) < numParam)
    stop("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (Rf_length(r_where) != dims)
    stop("Item has %d dimensions, but where is of length %d",
	  dims, Rf_length(r_where));

  int outcomes = spec[RPF_ISpecOutcomes];
  if (Rf_length(r_weight) != outcomes)
    stop("Item has %d outcomes, but weight is of length %d",
	  outcomes, Rf_length(r_weight));

  double *where = NULL;
  if (dims) where = REAL(r_where);

  const int numDeriv = numParam + numParam*(numParam+1)/2;
  SEXP ret;
  Rf_protect(ret = Rf_allocVector(REALSXP, numDeriv));
  memset(REAL(ret), 0, sizeof(double) * numDeriv);
  (*Glibrpf_model[id].dLL1)(spec.begin(), REAL(r_param),
			    where, REAL(r_weight), REAL(ret));
  for (int px=0; px < numDeriv; px++) {
    if (!std::isfinite(REAL(ret)[px])) stop("Deriv %d not finite at step 1", px);
  }
  (*Glibrpf_model[id].dLL2)(spec.begin(), REAL(r_param), REAL(ret));
  for (int px=0; px < numDeriv; px++) {
	  //if (!std::isfinite(REAL(ret)[px])) stop("Deriv %d not finite at step 2", px);
  }
  UNPROTECT(1);
  return ret;
}

// [[Rcpp::export]]
SEXP dTheta(const NumericVector &spec, SEXP r_param, SEXP r_where, SEXP r_dir)
{
	int id = getSpecID(spec);
  int numSpec = (*Glibrpf_model[id].numSpec)(spec.begin());
  if (spec.size() < numSpec)
    stop("Item spec must be of length %d, not %d", numSpec, spec.size());
    
  int numParam = (*Glibrpf_model[id].numParam)(spec.begin());
  if (Rf_length(r_param) < numParam)
    stop("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (dims == 0) stop("Item has no factors");
  if (Rf_length(r_dir) != dims)
    stop("Item has %d dimensions, but dir is of length %d",
	  dims, Rf_length(r_dir));
  if (Rf_length(r_where) != dims)
    stop("Item has %d dimensions, but where is of length %d",
	  dims, Rf_length(r_where));

  SEXP ret, names;
  Rf_protect(ret = Rf_allocVector(VECSXP, 2));
  Rf_protect(names = Rf_allocVector(STRSXP, 2));

  int outcomes = spec[RPF_ISpecOutcomes];
  SEXP grad, hess;
  Rf_protect(grad = Rf_allocVector(REALSXP, outcomes));
  Rf_protect(hess = Rf_allocVector(REALSXP, outcomes));
  memset(REAL(grad), 0, sizeof(double) * outcomes);
  memset(REAL(hess), 0, sizeof(double) * outcomes);
  (*Glibrpf_model[id].dTheta)(spec.begin(), REAL(r_param), REAL(r_where), REAL(r_dir),
			     REAL(grad), REAL(hess));
  SET_VECTOR_ELT(ret, 0, grad);
  SET_VECTOR_ELT(ret, 1, hess);
  SET_STRING_ELT(names, 0, Rf_mkChar("gradient"));
  SET_STRING_ELT(names, 1, Rf_mkChar("hessian"));
  Rf_namesgets(ret, names);
  UNPROTECT(4);
  return ret;
}

// [[Rcpp::export]]
NumericVector rescale(const NumericVector &spec, SEXP r_param, SEXP r_mean, const NumericMatrix &r_cov)
{
	int id = getSpecID(spec);
  int numSpec = (*Glibrpf_model[id].numSpec)(spec.begin());
  if (spec.size() < numSpec)
    stop("Item spec must be of length %d, not %d", numSpec, spec.size());
    
  int numParam = (*Glibrpf_model[id].numParam)(spec.begin());
  if (Rf_length(r_param) < numParam)
    stop("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (dims == 0) stop("Item has no factors");
  if (Rf_length(r_mean) != dims)
    stop("Item has %d dimensions, but mean is of length %d",
	  dims, Rf_length(r_mean));

  int cov_rows = r_cov.nrow();
	int cov_cols = r_cov.ncol();
  if (cov_rows != dims || cov_rows != dims)
    stop("Item has %d dimensions, but cov is %dx%d",
	  dims, cov_rows, cov_cols);

  Eigen::VectorXi mask(numParam);
  mask.setZero();

  NumericVector ret = clone(r_param);
  (*Glibrpf_model[id].rescale)(spec.begin(), ret.begin(), mask.data(),
															 REAL(r_mean), r_cov.begin());
  return ret;
}

// [[Rcpp::export]]
bool has_openmp()
{
#if defined(_OPENMP)
	return true;
#else
	return false;
#endif
}

int GlobalNumberOfCores = 1;

// [[Rcpp::export]]
int setNumberOfCores(int num)
{
#if defined(_OPENMP)
	GlobalNumberOfCores = num;
#endif
	return num;
}

extern const struct rpf librpf_model[];
extern const int librpf_numModels;

const struct rpf *Glibrpf_model;
int Glibrpf_numModels;

// Called by OpenMx
void get_librpf_models(int version, int *numModels, const struct rpf **model) // nocov start
{
  if (version != LIBIFA_RPF_API_VERSION) stop("LIBIFA_RPF binary API version mismatch");
  *numModels = librpf_numModels;
  *model = librpf_model;
} // nocov end

// [[Rcpp::export]]
void registerCCallable()
{
    Glibrpf_numModels = librpf_numModels;
    Glibrpf_model = librpf_model;
    R_RegisterCCallable("rpf", "get_librpf_model_GPL", (DL_FUNC) get_librpf_models);
}
