#include <cstdlib>
#include <utility>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rcpp.h>
#include <stdlib.h>
#include <string.h>
// [[Rcpp::plugins(cpp11)]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

enum Timing {exact, interval, none};

typedef struct {
  int *val;
  int nlevels;
  double *x;
} FACTOR;

typedef struct {
  double *par;
  int nlevels;
} FACTORPAR;


static bool debug = false;
// [[Rcpp::export]]
void cdebug(bool dodebug) {
  debug=dodebug;
}

// return log(sum(exp(x+y)))
inline double logsumofexp(const int n,const double *x, const double *y) {
  if(n < 1) stop("logsumofexp called with n < 1 (%d)",n);
  double cur = x[0]+y[0];
  for(int i = 1; i < n; i++) cur = R::logspace_add(cur,x[i]+y[i]);
  return cur;
}

// log(1+sum(x))
inline double log1sumx(int n, const double *x) {
  double lsum = 0;
  for(int i = 0; i < n; i++) lsum = R::logspace_add(lsum,x[i]);
  return lsum;
}



// convert parametrized probabilities Pi = exp(ai)/(1+sum(exp(ai))) to log probabilities
// a2logp <- function(a) {b <- c(0,a); logp <- b - logsumofexp(b,0); ifelse(is.na(logp),0,logp)}
void a2logp(const int n, const double *a, double *logp) {
  if(n < 0) stop("a2logp called with negative n (%d)",n);
  
  double lsum = log1sumx(n, a);
  logp[0] = -lsum;
  for(int i = 0; i < n; i++) logp[i+1] = a[i]-lsum;
  return;
}


// Taken from Martin Maechler https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
// "Accurately Computing log(1-exp(-a)), Assessed by the Rmpfr package"
// compute log(1-exp(-x))
inline double log1mexp(double x) {
  return (x <= M_LN2) ? log(-expm1(-x)) : log1p(-exp(-x));
}

/*
// log(1+exp(x))
inline double log1pexp(double x) {
  if(x <= -37) return exp(x);
  if(x <= 18) return log1p(exp(x));
  if(x <= 33.3) return x+exp(-x);
  return x;
}
*/

inline double log0(double x) { return (x < 0) ? -DBL_MAX : log(x); }

inline void obsloglik(const int tr, const Timing timing, const double *lh, double dur,
		      double **mup, const int npoints, int transitions,
		      const bool *riskmask,
		      double *llspell) {

  // If there is a transition, add the loghazard
  // divide by the survival prob up until now, i.e. subtract its log
  for(int j = 0; j < npoints; j++) {
    double logsumhaz = -DBL_MAX;
    for(int t = 0; t < transitions; t++) {
      if(riskmask && !riskmask[t]) continue;
      logsumhaz = R::logspace_add(logsumhaz, lh[t]+mup[t][j]);
    }
    switch(timing) {
    case exact:
      llspell[j] -= dur*exp(logsumhaz);
      break;
    case interval:
      if(tr <= 0) llspell[j] -= dur*exp(logsumhaz);
      break;
    case none:
      llspell[j] -= logsumhaz;
      break;
    }
    if(tr > 0) {
      switch(timing) {
      case interval:
	llspell[j] += log1mexp(dur*exp(logsumhaz)) - logsumhaz;
	// NOTE: no break here!
      case exact:
      case none:
	llspell[j] += lh[tr-1] + mup[tr-1][j];
	break;
      }
    } 
  }
}

// compute gradient
inline void gobsloglik(const int tr, const Timing timing, const double *lh, double dur, int obs,
		       double **mup, const int npoints, int transitions,
		       int npars, int *nfacs, int *faclevels,
		       const bool *riskmask,
		       int *nvars,
		       double **matp,
		       FACTOR **factors,
		       double *dllspell, const bool onlydist) {
  const int t = tr-1;
  for(int j = 0; j < npoints; j++) {
    double sumhaz = 0;
    double logsumhaz = -DBL_MAX;
    if(timing == none || (timing == interval && t >= 0)) {
      for(int tt = 0; tt < transitions; tt++) {
	if(riskmask && !riskmask[tt]) continue;
	logsumhaz = R::logspace_add(logsumhaz, lh[tt] + mup[tt][j]);
      }
      sumhaz = exp(logsumhaz);
    }

    double *dll = &dllspell[j*npars];
    int pos = 0;
    for(int tt = 0; tt < transitions; tt++) {
      if(riskmask && !riskmask[tt]) {pos += nvars[tt]+faclevels[tt]+npoints;continue;}
      const double haz = exp(lh[tt] + mup[tt][j]);
      const double funnyexpr = (timing==interval) ? 
	haz*(-dur*exp(-dur*sumhaz)/expm1(-dur*sumhaz) - 
	     1.0/sumhaz) + (tt==t) : 0.0;
      if(!onlydist) {
	const double *mat = &matp[tt][obs*nvars[tt]];

	for(int k = 0; k < nvars[tt]; k++) {
	  switch(timing) {
	  case exact:
	    dll[pos++] -= dur*haz*mat[k] - ((tt==t) ? mat[k] : 0.0);
	    break;
	  case interval:
	    dll[pos++] -= (t < 0) ? dur*haz*mat[k] : -mat[k]*funnyexpr;
	    break;
	  case none:
	    dll[pos++] -= haz*mat[k]/sumhaz - ((tt==t) ? mat[k] : 0.0);
	    break;
	  }
	}
	const FACTOR *fac = factors[tt];
	// loop through the factors for this transition
	for(int k = 0; k < nfacs[tt]; k++) {
	  const int fval = fac[k].val[obs];
	  if(fval <= 0 || ISNAN(fval)) {pos += fac[k].nlevels; continue;}
	  const double f = (fac[k].x != 0) ? fac[k].x[obs] : 1.0;
	  switch(timing) {
	  case exact:
	    dll[pos + fval-1] -= dur*haz*f - ((tt==t) ? f : 0);
	    break;
	  case interval:
	    dll[pos + fval-1] -= (t < 0) ? dur*haz*f : -f*funnyexpr;
	    break;
	  case none:
	    dll[pos + fval-1] -= haz*f/sumhaz - ((tt==t) ? f : 0);
	    break;
	  }
	  pos += fac[k].nlevels;
	}
      } else {
	pos += nvars[tt];
	for(int k = 0; k < nfacs[tt]; k++) pos += factors[tt][k].nlevels;
      }
      // The mu, the covariate is 1
      
      switch(timing) {
      case exact:
	dll[pos+j] -= dur*haz - ((tt==t) ? 1 : 0);
	break;
      case interval:
	dll[pos+j] -= (t < 0) ? dur*haz : -funnyexpr;
	break;
      case none:
	dll[pos+j] -= haz/sumhaz - ((tt==t) ? 1 : 0);
	break;
      }
      pos += npoints;
    }
  }
}

inline void updategradient(const int npoints, const double *dllspell, const double *llspell, 
			   const double *logprobs, const double ll,
			   const int transitions, const int npars, const int *nvars, const int *faclevels, 
			   const double *pargs, const int totalpars, double *spellgrad, const bool onlyprobs) {

  (void) memset(spellgrad, 0, totalpars*sizeof(*spellgrad));
  //  (void) memset(gkahanc, 0, totalpars*sizeof(*gkahanc));
  // compute gradient
  // we have the gradient of each llspell component in dllspell
  // Now, ll = log(sum(p_j exp(lh_j)))
  // where j ranges over 1..masspoints
  // We have lh_j in llspell
  // we have d lh_j / dx in dllspell
  // we have d ll / dx = 1/sum(p_j exp(lh_j)) * d sum(p_j exp(lh_j)) / dx.
  //
  // For x not a probability parameter, i.e. only occuring in lh_j
  // we have d sum(p_j exp(lh_j)) / dx = sum(p_j exp(lh_j) dlh_j/dx).
  // 
  // For the probability parameters y (which do not occur in each of the lh_j) we have
  // d sum(p_j exp(lh_j)) / dy = sum(dp_j/dy exp(lh_j))
  
  // then compute d sum(p_j exp(lh_j)) / dx.
  // First for the ordinary covariates where it equals sum(p_j exp(lh_j) dlh_j/dx)

  for(int j = 0; j < npoints; j++) {
    if(!onlyprobs) {
      const double *dll = &dllspell[j*npars];
      const double scale = exp(logprobs[j] + llspell[j] - ll);
    
      int pos = 0;
      for(int t = 0; t < transitions; t++) {
	for(int k = 0; k < nvars[t]+faclevels[t]; k++) {
	  spellgrad[pos] += scale * dll[pos];
	  pos++;
	}
	// the mu for this masspoint in this transition
	spellgrad[pos+j] += scale*dll[pos+j];
	pos += npoints; // next transition
      }
    }
    // The probability-parameters ak occur in all the probabilities
    // compute dPj / dak
    // The a's are in pargs, assume a0==0 so that Pk=exp(ak)/sum(exp(aj))
    // let ld = log(1+sum(exp(pargs))), the log denominator.
    // Can be moved out of j-loop. We trust the compiler to do that.
    // Stick to logs, 

    const double ld = log1sumx(npoints-1, pargs);
    const double lscale = llspell[j] - ll;
    for(int k = 0; k < npoints-1; k++) {
      const double ak = pargs[k];
      double dPdak;

      if(j == 0) {
	// for j=0, special case, it's dP0/dak
	dPdak = -exp(ak-2*ld + lscale);
      } else if(j == k+1) {
	// dPk/dak
	dPdak = exp(ak-ld + lscale) - exp(2*(ak-ld) + lscale);
      } else {
	// dPj/dak
	dPdak = -exp(logprobs[k+1] + logprobs[j] + lscale);
      }

      spellgrad[npars+k] += dPdak; 
    }
  }
}

inline void  updatefisher(int *gradfill, int fishblock, int totalpars, double *gradblock, 
			  int *nonzero, double *spellgrad, double *fisher, int *memfail) {
  // When computing the fisher matrix, we have individual spell
  // gradients in the gradblock matrix, we should add this
  // spell's gradient as a column.  When the gradblock matrix is
  // full, we dsyrk it into the fisher matrix. We do this in a
  // critical region, since all the threads use the same
  // gradblock and fisher matrix. This has the additional
  // bonus that only one thread runs the dsyrk at any time, so
  // we can use a parallel blas. The reason we collect gradients
  // in gradblock is that dsyrk is a lot faster than dsyr (rank
  // 1 update).

  if(*gradfill == fishblock) {
    // dsyrk the gradblock into fisher
    //       subroutine dsyrk (UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)
    // C = alpha A * A' + beta C
    
    // how many rows are nonzero?
    int nnrank = 0;
    for(int j = 0; j < totalpars; j++) {
      bool nonz = false;
      for(int k = 0; k < fishblock && !nonz; k++) {
	nonz |= (gradblock[k*totalpars + j] != 0.0);
      }
      if(nonz) nonzero[nnrank++] = j;
    }
    if(4*nnrank < 3*totalpars) {
      // Less than some limit, we dsyrk a smaller one
      
      double *smallblock = new (std::nothrow) double[fishblock*nnrank];
      if(smallblock == 0) {(*memfail)++; return;}

      for(int b = 0; b < fishblock; b++) {
	for(int k = 0; k < nnrank; k++) {
	  smallblock[b*nnrank + k] = gradblock[b*totalpars + nonzero[k]];
	}
      }
      // dsyrk it into smallfish
      const double alpha=-1, beta=0;
      
      double *smallfish = new (std::nothrow) double[nnrank*nnrank];
      if(smallfish == 0) {(*memfail)++; delete [] smallblock; return;}
      F77_CALL(dsyrk)("U","N",&nnrank,&fishblock,&alpha,smallblock,&nnrank,&beta,
		      smallfish, &nnrank);
      delete [] smallblock;
      // update fisher from smallfish
      for(int j = 0; j < nnrank; j++) {
	for(int k = 0; k <= j; k++) {
	  fisher[nonzero[j]*totalpars + nonzero[k]] += smallfish[j*nnrank + k];
	}
      }
      delete [] smallfish;
    } else {
      const double alpha = -1, beta = 1;
      F77_CALL(dsyrk)("U","N", &totalpars, &fishblock, &alpha, gradblock, 
		      &totalpars, &beta, fisher, &totalpars);
    }
    *gradfill = 0;
  }
  // append the spell gradient to the block
  for(int k = 0; k < totalpars; k++) gradblock[*gradfill * totalpars + k] = spellgrad[k];
  (*gradfill)++;
}

// [[Rcpp::export]]
NumericVector cloglik(List dataset, List pset, List control,
		      const bool gdiff=false, const bool dogradient=false, 
		      const bool dofisher=false, const bool onlyprobs=false, const bool onlydist=false) {
  const IntegerVector d = as<IntegerVector>(dataset["d"]);
  const IntegerVector id = as<IntegerVector>(dataset["id"]);
  const NumericVector duration = as<NumericVector>(dataset["duration"]);
  const NumericVector spellidx = as<NumericVector>(dataset["spellidx"]);
  const CharacterVector ctiming = as<CharacterVector>(dataset["timing"]);
  const Timing timing = (ctiming[0] == "exact") ? exact : (ctiming[0] == "interval" ? interval : none);

  const int nthreads=as<IntegerVector>(control["threads"])[0];
  const int nspells = spellidx.size() - 1;
  List parset = as<List>(pset["parset"]);
  const double *pargs = REAL(as<NumericVector>(pset["pargs"]));
  const int npoints = 1 + as<NumericVector>(pset["pargs"]).size();
  double *logprobs = (double*) R_alloc(npoints, sizeof(double)); 
  
  a2logp(npoints-1,pargs,logprobs);
  const int N = d.size();
  const int transitions = parset.size();

  List risklist = as<List>(dataset["risksets"]);
  const int nrisks = risklist.size();
  IntegerVector state = as<IntegerVector>(dataset["state"]);
  
  const int fishblock = as<IntegerVector>(control["fishblock"])[0];

  bool *riskmasks = 0;
  if(nrisks > 0) {
    riskmasks = (bool*) R_alloc(nrisks*transitions, sizeof(bool)); 
    memset(riskmasks,0,nrisks*transitions*sizeof(bool));
  }
  
  for(int i = 0; i < nrisks; i++) {
    IntegerVector v = as<IntegerVector>(risklist[i]);
    for(int t = 0; t < v.size(); t++) {
      riskmasks[i*transitions + v[t]-1] = true;
    }
  }
  
  // pointers into the data

  double **matp = (double **) R_alloc(transitions, sizeof(double*));
  
  int *nvars = (int*) R_alloc(transitions, sizeof(int)); 
  
  // number of factors in each transition

  int *nfacs = (int*) R_alloc(transitions, sizeof(int)); 

  // pointers to factor lists
  FACTOR **factors = (FACTOR**) R_alloc(transitions, sizeof(FACTOR*)); 
  List data = as<List>(dataset["data"]);
  for(int i = 0; i < transitions; i++) {
    NumericMatrix smat = as<NumericMatrix>(as<List>(data[i])["mat"]);
    matp[i] = REAL(smat);
    nvars[i] = smat.nrow();

    List facs = as<List>(data[i])["faclist"];
    nfacs[i] = facs.size();
    //    printf("nfacs[%d] = %ld\n",i,nfacs[i]);
    if(nfacs[i] == 0) continue;
    factors[i] = (FACTOR*) R_alloc(nfacs[i], sizeof(FACTOR)); 

    for(int j = 0; j < nfacs[i]; j++) {
      IntegerVector fac = as<IntegerVector>(facs[j]);
      factors[i][j].val = INTEGER(fac);
      factors[i][j].nlevels = as<CharacterVector>(fac.attr("levels")).size();
      NumericVector xv = as<NumericVector>(fac.attr("x"));
      if(xv.size() > 0) {
	if(xv.size() != N) stop("Factor interaction term length(%d) does not match dataset(%d)",
				xv.size(), N);
	factors[i][j].x = REAL(xv);
      } else {
	factors[i][j].x = 0;
      }
    }
  }

  // pointers into the parameters

  double **betap = (double **) R_alloc(transitions, sizeof(double *)); 
  double **mup = (double **) R_alloc(transitions, sizeof(double *)); 

  FACTORPAR **facpars = (FACTORPAR **) R_alloc(transitions, sizeof(FACTORPAR*)); 
  int *faclevels = (int*) R_alloc(transitions, sizeof(int)); 
  for(int i = 0; i < transitions; i++) {
    betap[i] = REAL(as<List>(parset[i])["pars"]);
    mup[i] = REAL(as<List>(parset[i])["mu"]);
    List Rfacs = as<List>(parset[i])["facs"];
    if(nfacs[i] != Rfacs.size()) stop("number of factor parameters(%d) does no match data(%d)",
				      Rfacs.size(),nfacs[i]);
    faclevels[i] = 0;
    if(nfacs[i] > 0) {
      facpars[i] = (FACTORPAR *) R_alloc(nfacs[i], sizeof(FACTORPAR)); 
      for(int j = 0; j < nfacs[i]; j++) {
	facpars[i][j].par = REAL(as<NumericVector>(Rfacs[j]));
	facpars[i][j].nlevels = as<NumericVector>(Rfacs[j]).size();
	faclevels[i] += facpars[i][j].nlevels;
      }
    }
  }

  const bool dograd = dogradient ? true : dofisher;

  // find the number of parameters
  int npars = 0;
  for(int t = 0; t < transitions; t++) npars += nvars[t] + faclevels[t];
// and mu's
  npars += npoints*transitions;

  const int totalpars = npars+npoints-1;

  double LL = 0.0;
  const int gradsize = dograd ? totalpars : 1;

  double *gradblock = 0;
  if(dofisher) gradblock = (double *) R_alloc(fishblock*gradsize, sizeof(double)); 
  int gradfill = 0;

  // Our fisher matrix is global, and updated inside a critial region
  // We could do a reduction(+) on it instead, but that may consume too much memory
  NumericMatrix *retfisher = 0;
  double *fisher;
  if(dofisher) {
    retfisher = new NumericMatrix(gradsize,gradsize);
    fisher = REAL(*retfisher);
  }

  // Then some thread private storage which are static in the parallel for below.
  // We don't want to allocate in the loop.
  // must be static to allocate thread private
  // we could do array reduction on the gradient, but it's not supported in
  // the current C-compiler for windows on cran. So do it manually.
  static double *lh, *totgrad;
  static double *spellgrad, *llspell, *dllspell, *gneuc; 
  int *nonzero; // only used temporarily in critical section, so not thread local
  if(dofisher) nonzero = (int*) R_alloc(totalpars, sizeof(int)); 

#pragma omp threadprivate(llspell,dllspell, lh, spellgrad, totgrad, gneuc)
#pragma omp parallel num_threads(nthreads)
  {  
    // Can't do R_alloc in threads
    //    if(dograd) gkahanc = new double[gradsize]();
    totgrad = new double[gradsize]();
    llspell = new double[npoints];
    lh = new double[transitions];
    if(dograd) {
      dllspell = new double[npoints*npars];
      spellgrad = new double[totalpars];
      gneuc = new double[totalpars]();
    }
  }

  // Remember not to use any R-functions (or allocate Rcpp storage) inside the parallel region.
  double neuc = 0.0;
  int memfail = 0;
#pragma omp parallel for reduction(+: LL, neuc) num_threads(nthreads) schedule(guided)
  for(int spellno = 0; spellno < nspells; spellno++) {
    if(memfail > 0) continue;
    memset(llspell, 0, npoints*sizeof(*llspell));
    if(dograd) memset(dllspell, 0, npoints*npars*sizeof(*dllspell));
    for(int i = spellidx[spellno]; i < spellidx[spellno+1]; i++) {
      // compute the log hazard for this observation for each masspoint and transition
      // loop through the mass points. For each find the hazard sum
      // fill in the loghazards in the lh-array
      const bool *riskmask = nrisks>0 ? &riskmasks[(state[i]-1)*transitions] : 0;
      for(int t = 0; t < transitions; t++) {
	lh[t] = 0.0;
	if(riskmask && !riskmask[t]) continue;  
	const double *mat = &matp[t][i*nvars[t]];
	const double *beta = betap[t];
	for(int k = 0; k < nvars[t]; k++) {
	  lh[t] += mat[k]*beta[k];
	}
	const FACTOR *fac = factors[t];
	const FACTORPAR *fpar = facpars[t];
	for(int j = 0; j < nfacs[t]; j++) {
	  const int fval = fac[j].val[i];
	  if(fval <= 0 || ISNAN(fval)) continue;  // skip NA-levels, i.e. reference
	  double x = (fac[j].x == 0) ? 1.0 : fac[j].x[i];
	  lh[t] += fpar[j].par[fval-1] * x;
	}
      }

      // update llspell with the observation log likelihood 
      obsloglik(d[i], timing, lh, duration[i], mup, npoints, transitions, riskmask, llspell);
      // update dllspell with the gradient of the observation log likelihood
      if(dograd && !onlyprobs) 
	gobsloglik(d[i], timing, lh, duration[i], i, mup, npoints, transitions, npars, nfacs,
		   faclevels, riskmask, nvars, matp, factors, dllspell, onlydist);
    }

    // We have collected the loglikelihood of a spell, one for each masspoint
    // integrate it with the probabilities
    double ll;
    if(gdiff) {
      // Use Neumaier summation to minimize roundoff errors
      ll = logsumofexp(npoints-1,llspell,logprobs);
      const double v = expm1(llspell[npoints-1] - ll);
      const double t = LL+v;
      neuc += (abs(LL) > abs(v)) ? (LL-t) + v : (v - t) + LL;
      LL = t;
    } else {
      // compute the log likelihood
      ll = logsumofexp(npoints,llspell,logprobs);
      // for many spells we should perhaps do a compensated addition, Kahan/Neumaier?
      const double t = LL+ll;
      neuc += (abs(LL) > abs(ll)) ? (LL-t)+ll : (ll-t)+LL;
      LL = t;
      
      if(dograd) {
	// compute the spell gradient, update the gradient
	updategradient(npoints, dllspell, llspell, logprobs, ll,
		       transitions, npars, nvars, faclevels, pargs, totalpars, spellgrad, onlyprobs);
	// update global gradient with spell gradient
	for(int k = 0; k < totalpars; k++) {
	  const double t = totgrad[k]+spellgrad[k];
	  gneuc[k] += (abs(totgrad[k]) > abs(spellgrad[k])) ? (totgrad[k]-t)+spellgrad[k] :
	    (spellgrad[k] - t) + totgrad[k];
	  totgrad[k] = t;
	}
	if(dofisher) {
	  // update the fisher matrix from the spellgrad
	  // we use a global fisher matrix, no omp reduction, so do it in a critical section
#pragma omp critical
	  updatefisher(&gradfill, fishblock, totalpars, gradblock, nonzero, spellgrad, fisher, &memfail);
	}
      }
    } 
  }  // end of spell loop


  LL += neuc;
  // Deallocate thread private storage
#pragma omp parallel num_threads(nthreads)
  {  
    delete [] llspell;
    delete [] lh;
    if(dograd) {
      for(int k = 0; k < gradsize; k++) {totgrad[k] += gneuc[k];}
      delete [] dllspell;
      delete [] spellgrad;
      delete [] gneuc;
    }
  }

  if(memfail > 0) {delete [] totgrad; stop("Memory allocation failed");}

  // Then set up the return value
  NumericVector ret = NumericVector::create(LL);

  if(dograd) {
    // Reduce the totgrad in a critical section
    NumericVector retgrad(totalpars);
    double *grad = REAL(retgrad);
    double *neu = new double[gradsize]();
#pragma omp parallel num_threads(nthreads) 
    {
#pragma omp critical
      for(int k = 0; k < gradsize; k++) {
	const double t = grad[k] + totgrad[k];
	neu[k] += (abs(grad[k]) > abs(totgrad[k])) ? (grad[k] - t) + totgrad[k] :
	  (totgrad[k]-t) + grad[k];
	grad[k] = t;
      }
      delete [] totgrad;
    }
    for(int k = 0; k < gradsize; k++) grad[k] += neu[k];
    delete [] neu;
    ret.attr("gradient") = retgrad;
  }


  // if anything remains in the gradient blocks, dsyrk it into the fisher matrix
  if(dofisher && gradfill > 0) {
    const double alpha = -1, beta = 1;
    const int N = totalpars, K=gradfill;
    F77_CALL(dsyrk)("U","N", &N, &K, &alpha, gradblock, &N, &beta, fisher, &N);
 
  }

  if(dofisher)  {
    // fill the lower half of the fisher matrix
    // read consecutively, write with stride, that's typically faster
    for(int i = 0; i < totalpars; i++) {
      for(int j = 0; j < i; j++) {
	fisher[j*totalpars + i] = fisher[i*totalpars + j];
      }
    }
    ret.attr("fisher") = *retfisher;
  }
  return ret;
}

// [[Rcpp::export]]
List genspell(double x1,double x2, double ve, double vp, double censor) {
  std::vector<double> x1s,x2s,alphas,ds,ts;
  bool done = false, onp = false;
  double alpha = 0, newalpha=0, newx1=x1, newx2=x2, d,t;
  int i = 0;
  double cumtime = R::runif(0,censor); // draw start time when individual enters dataset
  while(!done) {
    alpha = newalpha; x1=newx1; x2=newx2;
    double te = -log(R::runif(0,1))*exp(-(x1-x2+ve+0.2*alpha));
    double tp = onp ? DBL_MAX : -log(R::runif(0,1))*exp(-(x1+0.5*x2+vp));
    double tc = -log(R::runif(0,1))*35.0;
    if(tc < te && tc < tp) {
      newx1 = x1 + R::rnorm(0,1);
      newx2 = x2 + R::rnorm(0,1);
      d = 0;
      t = tc;
    } else if(te < tc && te < tp) {
      d = 1;
      t = te;
      done = true;
    } else {
      d = 2;
      onp = true;
      t= tp;
      newalpha = 1;
    }
    cumtime += t;
    if(cumtime > censor) {done = true; d = 0; t = censor - (cumtime-t);}
    x1s.push_back(x1); x2s.push_back(x2); alphas.push_back(alpha); ds.push_back(d); ts.push_back(t);
    i++;
  }
  
  List ret = List::create(Named("x1") = x1s,
			  Named("x2") = x2s,
			  Named("alpha") = alphas,
			  Named("d") = ds,
			  Named("duration") = ts);
  ret.attr("class") = "data.frame";
  return ret;
}
