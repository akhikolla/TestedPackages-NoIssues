#include <Rcpp.h>
using namespace Rcpp;

// Parallel computation and random numbers

// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#define thisThread omp_get_thread_num()
#else
#define thisThread 0
#endif

// [[Rcpp::depends(sitmo)]]
#include <sitmo.h>

// [[Rcpp::depends(BH)]]
#include <boost/random/normal_distribution.hpp>

#define threadRNG (parallel::rngs.r[thisThread])
namespace parallel {
    
    struct RNGS {
        uint32_t nthreads, seed;
        std::vector< sitmo::prng_engine * > r;
        
        RNGS() : nthreads(1), seed(1) { 
            r.push_back(new sitmo::prng_engine(seed));
        }

        ~RNGS() {
            for (uint32_t i = 0; i<r.size(); i++) delete r[i];
        }

    
        void setPThreads(uint32_t nt) {
#ifdef _OPENMP
            nthreads = nt;
            for (uint32_t i=r.size(); i<nthreads; i++) {
	        uint32_t thisSeed = (seed+i) % sitmo::prng_engine::max();
	        r.push_back(new sitmo::prng_engine(thisSeed));
            }
#else
            Rcpp::warning("No openmp support");
#endif
        }

        void setPSeed(double s) {
            if ((s<=0.0) || (s>=1.0)) 
                Rcpp::stop("seed must be between 0 and 1");
            seed = s * sitmo::prng_engine::max();
            for (uint32_t i=0; i<r.size(); i++) {
	        uint32_t thisSeed = (seed+i) % sitmo::prng_engine::max();
	        r[i]->seed(thisSeed);
            }
        }

    };

    RNGS rngs;

}

namespace sa {
    
    struct score {
	score() {}
	virtual ~score() {}
	virtual double operator()(double h) = 0;
    };

    // solver for P(score>target)=1-beta
    double qsolver(double h, score &s, double target, double beta,
		   double gain, int init, int iter) {
	double hpr=0.0;
	for (int i=-init+1; i <= iter; i++) {
            R_CheckUserInterrupt();
	    double sh = (s(h)>target) ? beta : beta-1; 
	    h = std::max(0.0,h-gain*sh);
	    if (i>0) hpr += (h-hpr)/i;
	}
	return hpr;
    }

    // solver for E(score)=target
    double asolver(double h, score &s, double target, double alpha,
		   double gain, double q, int init, int iter) {
	double hpr=0.0;
	for (int i=-init+1; i <= iter; i++) {
	    double sh = (s(h)-target)/target;
	    h = std::max(0.0,h-gain*sh/pow(std::max(i,1),q));
	    if (i>0) hpr += (h-hpr)/i;
	}
	return hpr;
    }

}


namespace {

    // a struct for simulating x, xbarbar0, s2barbar0)
    struct xbs {
	int m, tau;
	double sm, xdelta, sdelta;

	xbs(int m, int tau=NA_INTEGER, double xdelta=0.0, double sdelta=1.0) :
	    m(m), tau(tau), sm(1/sqrt(static_cast<double>(m))),
	    xdelta(xdelta), sdelta(sdelta) {}

	double operator()(int i) {
	    boost::random::normal_distribution<double> normal;
	    sitmo::prng_engine *rng = threadRNG;
	    double u = normal(*rng);
	    if ((tau != NA_INTEGER) && (i>=tau)) u = xdelta+sdelta*u;
	    return u;
	}


	void xs2barbar0(double &xbarbar0, double &s2barbar0) {
	    boost::random::normal_distribution<double> normal;
	    sitmo::prng_engine *rng = threadRNG;
	    xbarbar0 = sm*normal(*rng);
	    s2barbar0 = 0;
	    for (int i=1; i<m; i++) {
		double u = normal(*rng);
		s2barbar0 += u*u;
	    }
	    s2barbar0 /= m-1;
	}
    }; 


    // Self-learning limits
    struct sllimits {
	bool sl;
	int m, d;
	double Linf, Delta, A, B, dm, q, k, 
	    xbb, s2bb, muhat, s2hat, shat, lhat;

	sllimits(double *lim, double m0, double s20) :
	    sl(::R_finite(lim[2])), m(std::floor(lim[4]+0.5)),
            Linf(lim[0]), Delta(lim[1]), A(lim[2]), B(lim[3]), dm(lim[4]),
            xbb(m0), s2bb(s20) {
	    setlim();
	}

	void setlim() {
	    d = 0;
	    q = 0;
	    muhat = xbb;
	    s2hat = s2bb;
	    shat = sqrt(s2bb);
	    lhat = Linf+Delta*sqrt(std::min(1.0,dm/m));
	}

	void update(double xi) {
	    if (sl) {
		// estimates update
		double r=xi-xbb;
		m++; 
		xbb += r/m;
		s2bb = (m-2)*s2bb/(m-1)+r*r/m ;
		// update the limits?
                d++;
		q = q+(xi-muhat)*(xi-muhat)/s2hat;
		if (q<A*d-B) setlim();
	    }
	}
	
    }; // end sllimits


    struct Chart {
	bool sim;
        int lstat;
	double *limit;
	
	Chart(bool sim, int lstat, double *limit) :
	    sim(sim), lstat(lstat), limit(limit) {}

        virtual ~Chart() {}

	// used only when when sim=false
	virtual double carl(xbs &sampler, double m0, double s20) {
	    Rcpp::stop("the carl method is not defined");
	}

        virtual bool update(int i, double x, sllimits &sl, double *stat) = 0;
	
	int crl(xbs &sampler, double m0, double s20, int maxrl) {
            int i;
            double stat[8];
	    sllimits sl(limit,m0,s20);
	    for (i=1; i<maxrl; i++) {
		double xi = sampler(i);
		if (update(i, xi, sl, stat)) break;
		sl.update(xi);
	    } 
	    return i;
        }

    };


    struct ShewhartX : public Chart {
	ShewhartX(double *limit) :
	    Chart(::R_finite(limit[2]), 7, limit) {}

	double carl(xbs &sampler, double m0, double s20) {
	    double m = (sampler.xdelta-m0)/sampler.sdelta;
	    double h = (limit[0]+limit[1])*sqrt(s20)/sampler.sdelta;
	    return 1.0/(R::pnorm(-h,m,1.0,1,0)+R::pnorm(h,m,1.0,0,0));
	}

        bool update(int i, double x, sllimits &sl, double *stat) {
            stat[0] = x;
            stat[1] = sl.muhat;
            stat[2] = sl.muhat-sl.lhat*sl.shat;
            stat[3] = sl.muhat+sl.lhat*sl.shat;
            stat[4] = sl.muhat;
            stat[5] = sl.shat;
            stat[6] = sl.lhat;
            return (x<stat[2]) || (x>stat[3]);
        } 

    };


    struct EWMA : public Chart {
	double lambda, se;

	EWMA(double lambda, double *limit) :
	    Chart(true, 7, limit), lambda(lambda), se(sqrt(lambda/(2.0-lambda)))
        {}

        bool update(int i, double x, sllimits &sl, double *stat) {
            if (i==1) stat[0] = sl.muhat;
            stat[0] += lambda*(x-stat[0]);
            stat[1] = sl.muhat;
            stat[2] = sl.muhat-se*sl.lhat*sl.shat;
            stat[3] = sl.muhat+se*sl.lhat*sl.shat;
            stat[4] = sl.muhat;
            stat[5] = sl.shat;
            stat[6] = sl.lhat;
            return (stat[0]<stat[2]) || (stat[0]>stat[3]);
        } 

    };


    struct CUSUM : public Chart {
	double k;
	
	CUSUM(double k, double *limit) :
	    Chart(true, 8, limit), k(k) {}

        bool update(int i, double x, sllimits &sl, double *stat) {
            if (i==1) stat[0] = stat[1] = stat[2] = 0;
            double r = (x-sl.muhat)/sl.shat;
            stat[0] = std::min(0.0, stat[0]+r+k);
            stat[1] = std::max(0.0, stat[1]+r-k);
            stat[3] = -sl.lhat;
            stat[4] = sl.lhat;
            stat[5] = sl.muhat;
            stat[6] = sl.shat;
            stat[7] = sl.lhat;
            return (stat[0]<stat[3]) || (stat[1]>stat[4]);
        } 

    };

    
    inline Chart &getChart(List chart) {
	Chart *c;
	std::string t=as<std::string>(chart["chart"]);
	NumericVector l=as<NumericVector>(chart["limit"]);
	if (t=="X") {
	    c = new ShewhartX(l.begin());
	} else if (t=="EWMA") {
	    c = new EWMA(as<double>(chart["lambda"]),l.begin());
	} else if (t=="CUSUM") {
	    c = new CUSUM(as<double>(chart["k"]),l.begin());
	} else {
	    Rcpp::stop("Unknown chart");
	}
	return *c;
    }


    inline void simrl(Chart &c, xbs &s, double m0, double s20, 
                      int nrl, int *rl, int maxrl) {
#ifdef _OPENMP
        size_t nt = std::min(static_cast<uint32_t>(nrl/5),
                             parallel::rngs.nthreads);
#pragma omp parallel for num_threads(nt)
#endif
	for (int i=0; i<nrl; i++) {
	    rl[i] = c.crl(s,m0,s20,maxrl);
        }       
    } 
    
    
    struct aScore : public sa::score {
	Chart &c;
	xbs &s;

	aScore(Chart &c, xbs &s) : c(c), s(s) {}

	double operator()(double h) {
	    double m0, s20;
	    c.limit[1] = h;
	    s.xs2barbar0(m0, s20);
	    return c.carl(s,m0,s20);
	}
    }; 

    struct sScore : public sa::score {
	int nrl, *rl, maxrl;
	Chart &c;
	xbs &s;

	sScore(int nrl, int *rl, int maxrl, Chart &c, xbs &s) :
	    nrl(nrl), rl(rl), maxrl(maxrl), c(c), s(s) {}

       ~sScore() {}

	double operator()(double h) {
	    double m0, s20, carl;
	    c.limit[1] = h;
	    s.xs2barbar0(m0, s20);
            simrl(c, s, m0, s20, nrl, rl, maxrl);
            carl = std::accumulate(rl, rl+nrl, 0.0);
            return carl/nrl;
        }

    };
	

}


// [[Rcpp::export]]
bool hasOMP() {
    bool ans;
#ifdef _OPENMP
    ans = true;
#else
    ans = false;
#endif
    return ans;
}

// [[Rcpp::export]]
void setOMPThreads(uint32_t nthreads) {
    parallel::rngs.setPThreads(nthreads);
}

// [[Rcpp::export]]
void setSITMOSeeds(double seed) {
    parallel::rngs.setPSeed(seed);
}


// [[Rcpp::export]]
List mkChart(List chart, int m, double A, double B, double arl0, double Linf,
             double alpha=0.10, double beta=0.05, int H=200,
             int Ninit=1000, int Nfinal=30000) {
    List ans = clone(chart);
    ans["limit"] = NumericVector::create(_["Linf"]=Linf,
					 _["Delta"]=Linf/10,
					 _["A"]=A,
					 _["B"]=B,
                                         _["m"]=m);
    Chart &c=getChart(ans);
    double *limit = c.limit;
    xbs s(m); 
    if (c.sim) {
        IntegerVector rl(H);
	sScore sc(H, rl.begin(), 100*arl0, c, s);
	limit[1] = sa::qsolver(limit[1], sc, (1-alpha)*arl0, beta,
                               Linf/100, Ninit, Nfinal);
    } else {
	aScore sc(c,s);
	limit[1] = sa::qsolver(limit[1], sc, (1-alpha)*arl0, beta,
                               Linf/100, Ninit, Nfinal);
    }
    delete &c;
    return ans;
}

// [[Rcpp::export]]
NumericMatrix ruv(int n, int m) {
    if (n<0) Rcpp::stop("n cannot be negative");
    if (m<2) Rcpp::stop("m must be greater than 1");
    NumericMatrix ans(n,2);
    xbs s(m);
    double m0, s20, 
           a = sqrt(static_cast<double>(m)), b = sqrt((m-1.0)/2.0);
    for (int i=0; i<n; i++) {
        s.xs2barbar0(m0, s20);
        ans(i,0) = a*m0;
        ans(i,1) = b*(s20-1);
    }
    return ans;
}


// [[Rcpp::export]]
IntegerVector rcrl(int n, List chart, double u, double v, 
                   int tau, double delta, double omega, int maxrl=1000000) {
    if (n<0) Rcpp::stop("n cannot be negative");
    IntegerVector rl(n);
    Chart &c=getChart(chart);
    int m = std::floor(c.limit[4]+0.5);
    double m0 = u/sqrt(c.limit[4]), s20 = 1.0+v*sqrt(2.0/(c.limit[4]-1));
    xbs s(m, tau, delta, omega);
    simrl(c, s, m0, s20, n, rl.begin(), maxrl);
    delete &c;
    return rl;
}

// [[Rcpp::export]]
NumericMatrix applyChart(List chart, NumericVector x, double mu0, double s0) {
    Chart &c=getChart(chart);
    int n=x.size(), l=c.lstat;
    sllimits sl(c.limit, mu0, s0*s0);
    NumericMatrix ans(n,l);
    double stat[8];
    for (int i=0; i<n; i++) {
        c.update(i+1, x[i], sl, stat);
        sl.update(x[i]);
        for (int j=0; j<l; j++) ans(i,j) = stat[j];
    }
    delete &c;
    return ans;
}

