// The following includes 4 functions
// logsWeibloglambda
// logsWeib
// sWeibloglambda
// sWeib
// Footnote 1
#include <algorithm>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp ;

double logsWeibloglambda(double t, double alpha, double loglambda){
	if(t>0.0){
		return -exp(loglambda+alpha*log(t));
	}else{
		return 0.0;
	}
}

double logsWeib(double t, double alpha, double lambda){
	if(lambda>0.0){
			return logsWeibloglambda(t,alpha,log(lambda));
	}else{
			return 0.0;
	}
}

double sWeibloglambda(double t, double alpha, double loglambda){
	return exp(logsWeibloglambda(t,alpha,loglambda));
}

double sWeib(double t, double alpha, double lambda){
	if(t==0.0){
		return 1.0;
	}else{
			return exp(logsWeib(t,alpha,lambda));
	}

}
// The following includes 4 functions
// logdWeibloglambda
// logdWeib
// dWeibloglambda
// dWeib
// Footnote2

double logdWeibloglambda(double t, double alpha, double loglambda){
	if(t>0.0){
	return loglambda+log(alpha)+(alpha-1.0)*log(t)+logsWeibloglambda(t,alpha,loglambda);		
	}else{
		return -pow(100.0,100.0);
	}
}

double logdWeib(double t, double alpha, double lambda){
	if(lambda>0.0){
	return logdWeibloglambda(t,alpha,log(lambda));		
	}else{
		return -pow(100.0,100.0);
	}
}

double dWeibloglambda(double t, double alpha, double loglambda){
  return exp(logdWeibloglambda(t,alpha,loglambda));
}

double dWeib(double t, double alpha, double lambda){
	return exp(logdWeib(t,alpha,lambda));
}
// The following includes 2 functions
// dWeib
// logdWeib

double pWeib(double t, double alpha, double lambda){
	return 1.0-sWeib(t,alpha,lambda);
}

double logpWeib(double t, double alpha, double lambda){
	return log(pWeib(t,alpha,lambda));
}

double F1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta){
return 1.0-pow(1.0-p*pWeib(t,alpha,lambda),exp(sum(x*beta)));
	}
double F1v2(double t, double p, double alpha, double lambda, double xbeta){
return 1.0-pow(1.0-p*pWeib(t,alpha,lambda),exp(xbeta));
	}

double logf1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta){
	return sum(x*beta)+(exp(sum(x*beta))-1.0)*log(1.0-exp(log(p)+logpWeib(t,alpha,lambda)))
	+log(p)+logdWeib(t,alpha,lambda);
}
double logf1v2(double t, double p, double alpha, double lambda, double xbeta){
	return xbeta+(exp(xbeta)-1.0)*log(1.0-exp(log(p)+logpWeib(t,alpha,lambda)))
	+log(p)+logdWeib(t,alpha,lambda);
}

double f1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta){
	return exp(logf1(t,p,alpha,lambda,x,beta));
}


double F2(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta1, NumericVector beta2){
	return pow(1.0-p,exp(sum(x*beta1)))*(1.0-pow(sWeib(t,alpha,lambda),exp(sum(x*beta2))));
	}

double F2v2(double t, double p, double alpha, double lambda, double xbeta1, double xbeta2){
	//return exp(exp(xbeta1)*log(1.0-p)+log(1.0-pow(sWeib(t,alpha,lambda),exp(xbeta2))));
	return pow(1.0-p,exp(xbeta1))*(1.0-pow(sWeib(t,alpha,lambda),exp(xbeta2)));
	}
double logf2(double t, double p, double alpha, double lambda,
 NumericVector x, NumericVector beta1, NumericVector beta2){
	return exp(sum(x*beta1))*log(1.0-p)+
	(exp(sum(x*beta2))-1.0)*logsWeib(t,alpha,lambda)+
	sum(x*beta2)+logdWeib(t,alpha,lambda);
}

double logf2v2(double t, double p, double alpha, double lambda, double xbeta1, double xbeta2){
	return exp(xbeta1)*log(1.0-p)+(exp(xbeta2)-1.0)*logsWeib(t,alpha,lambda)
	+xbeta2+logdWeib(t,alpha,lambda);
}

double f2(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta1, NumericVector beta2){
	return exp(logf2(t,p,alpha,lambda,x,beta1,beta2));
}

double f1v2(double t, double p, double alpha, double lambda, double xbeta){
	return exp(logf1v2(t,p,alpha,lambda,xbeta));
}

double f2v2(double t, double p, double alpha, double lambda, double xbeta1, double xbeta2){
	return exp(logf2v2(t,p,alpha,lambda,xbeta1,xbeta2));
}

double logF1(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta){
	return log(F1(t,p,alpha,lambda,x,beta));
}

double logF2(double t, double p, double alpha, double lambda, NumericVector x, NumericVector beta1, NumericVector beta2){
	return exp(sum(x*beta1))*log(1.0-p)+log(1.0-pow(sWeib(t,alpha,lambda),exp(sum(x*beta2))));
}

double logScomp(double t, double alpha1, double lambda1,
double alpha2, double lambda2, double xbeta1, double xbeta2, double p){
		double F1=1.0-pow(1.0-p*pWeib(t,alpha1,lambda1),exp(xbeta1));
		double F2=pow(1.0-p,exp(xbeta1))*
		  (1.0-pow(sWeib(t,alpha2,lambda2),exp(xbeta2)));
		  if (F1+F2<1.0){
			  return log(1.0-F1-F2);
		  }else{
			 return -10000.0;
		  }

}

double timedWeibloglambda(double t, double ts, double alpha, double loglambda, double beta){
		if(t<ts){
			return logdWeibloglambda(t,alpha,loglambda);
		}else{
					return log(alpha)+loglambda+(alpha-1.0)*log(t)+beta
					-exp(loglambda+alpha*log(ts))-exp(loglambda+beta+alpha*log(t))
					+exp(loglambda+beta+alpha*log(ts));
		}
}

double timesWeibloglambda(double t, double ts, double alpha, double loglambda, double beta){
		if(t<ts){
			return logsWeibloglambda(t,alpha,loglambda);
		}else{
			if(ts>0.0){
				return -exp(loglambda+alpha*log(ts))
				-exp(loglambda+beta+alpha*log(t))
				+exp(loglambda+beta+alpha*log(ts));
			}else{
				return -exp(loglambda+beta+alpha*log(t));
			}
		}
}

