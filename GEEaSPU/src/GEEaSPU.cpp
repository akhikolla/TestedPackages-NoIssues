# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export()]]
arma::vec rowMinsC(arma::mat x) {
	int n = x.n_rows;
	arma::vec v(n);
	for (int i=0; i<n; i++) {
		v(i) = x.row(i).min();
	}
	return(v);
}





// [[Rcpp::export]]
double count_if(arma::uvec x) {
	double counter = 0;
    for(int i = 0; i < x.size(); i++) {
        if(x[i] == TRUE) {
            counter++;
        }
    }
    return counter;
}







// [[Rcpp::export]]
arma::vec rankC(NumericVector x) {
//   if (is_true(any(duplicated(x)))) {
//     Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
//   }
	NumericVector sorted = clone(x).sort();
	return as<arma::vec>(wrap(match(x, sorted)));
}







// [[Rcpp::export]]
NumericVector which_min(NumericVector y, LogicalVector x) {
    for(int i = 0; i < x.size(); i++) {
        if(x[i] == TRUE) {
            y[i]=0;
        }
    }
    NumericVector yy(y);
    yy = pow(y, 0.5);
    return(yy);
}












// [[Rcpp::export]]
List getEigen(arma::mat sigma) {
	
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, sigma);
	List E ;
	E["eigval"] = eigval;
	E["eigvec"] = eigvec;
	return(E) ;
}












// [[Rcpp::export()]]
arma::mat a8 (arma::mat U, int p, int k) { //p is #SNPs: k is #traits
	
	U.reshape(p, k) ;
	return U.t();  //as<NumericMatrix>(wrap(Ts));
}








// [[Rcpp::export()]]
double signC(double x) {
	
    if (x > 0) {
        return 1;
    } else if (x == 0) {
        return 0;
    } else {
        return -1;
    }
}



// [[Rcpp::export]]
arma::vec rcpp_standPow(arma::mat U1, double power) {
	
    const int k = U1.n_rows;
    arma::vec Ts(k);
    Ts.zeros();
	
    for (int t = 0; t < k; t++){
        double a = accu(pow(U1.row(t), power));
        Ts(t) = signC(a) * std::pow(std::abs(a), 1.0/power);
    }
	  
    return Ts;
}









// [[Rcpp::export()]]
arma::vec InfU(arma::mat U) {
    int k = U.n_rows;
    arma::vec Ts(k);
    Ts.zeros();
	  
    for (int t = 0; t < k; t++){
        arma::mat U0 = abs(U.row(t));
        Ts(t) = U0.max();
    }
    return Ts;
}






// [[Rcpp::export]]
arma::mat gauss_score(arma::mat invR, arma::mat G, arma::mat res, int n, int k, int p){
	
	arma::mat out(k,p);
	out.fill(0);
	 
	for (int i = 0; i < n; i++) {
    	out += (res.row(i) *invR).t() * G.row(i);
	}
	return vectorise(out, 1);
}





// [[Rcpp::export]]
List gauss_score_cov(arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp){

	arma::mat out(k,p);
	arma::mat V(tp, tp);
	out.fill(0);
	V.fill(0);
	
	for (int i = 0; i < n; i++) {
		 out += (res.row(i) *invR).t() * G.row(i);
		 arma::mat w = x.rows(k*i, k*(i+1)-1);
		 V += w.t() *invR* covres *invR * w;
	}

	arma::mat I11 = V.submat(0, 0, tp-p*k-1, tp-p*k-1) ;
	arma::mat I12 = V.submat(0, tp-p*k, tp-p*k-1, tp-1) ;
	arma::mat I21 = V.submat(tp-p*k, 0, tp-1, tp-p*k-1) ;
	arma::mat I22 = V.submat(tp-p*k, tp-p*k, tp-1,tp-1) ;

	arma::mat score_cov = I22 - I21 * pinv(I11) * I12 ;
	arma::mat score = (vectorise(out, 1) * arma::pinv(score_cov)) *  vectorise(out, 1).t();

	List S ;
	S["score_vec"] = vectorise(out, 1);
	S["score_cov"] = score_cov;
	S["score"] = score;
	return(S) ;
}







// [[Rcpp::export]]
arma::mat bin_score(arma::vec va, arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp){
	
	arma::mat out(1,tp);
	out.fill(0);
	 
	for (int i = 0; i < n; i++) {
        arma::vec v = va.subvec(k*i, k*(i+1)-1);
        arma::mat inV = diagmat(sqrt(1.0/v)) * invR* diagmat(sqrt(1.0/v));
        arma::mat w(k, x.n_cols);
        for (int j = 0; j < k; j++){
            w.row(j) = x.row(k*i+j) *v[j];
        }
        out += res.row(i) * inV *w;
    }
	return out.cols(tp- p*k, tp-1); 
}






// [[Rcpp::export]]
List bin_score_cov(arma::vec va, arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp){
	
	arma::mat out(1,tp);
	arma::mat V(tp, tp);
	out.fill(0);
	V.fill(0);

	for (int i = 0; i < n; i++) {
        arma::vec v = va.subvec(k*i, k*(i+1)-1);
        arma::mat inV = diagmat(sqrt(1.0/v)) * invR* diagmat(sqrt(1.0/v));
        arma::mat w(k,tp);
        for (int j = 0; j < k; j++){
            w.row(j) = x.row(k*i+j) *v[j];
        }
        out += res.row(i) * inV *w;
        V += w.t() * inV * covres *inV * w;
    }
	arma::mat I11 = V.submat(0, 0, tp-p*k-1, tp-p*k-1) ;
	arma::mat I12 = V.submat(0, tp-p*k, tp-p*k-1, tp-1) ;
	arma::mat I21 = V.submat(tp-p*k, 0, tp-1, tp-p*k-1) ;
	arma::mat I22 = V.submat(tp-p*k, tp-p*k, tp-1,tp-1) ;

	arma::mat score_cov = I22 - I21 * pinv(I11) * I12 ;
	arma::mat score = (out.cols(tp- p*k, tp-1) * arma::pinv(score_cov)) *  out.cols(tp- p*k, tp-1).t();

	List S ;
	S["score_vec"] = out.cols(tp- p*k, tp-1); 
	S["score_cov"] = score_cov;
	S["score"] = score;
	return(S) ;
}





// [[Rcpp::export]]
arma::mat rcpp_spuval(arma::mat U, arma::mat po, arma::mat po2, int p, int k){
    const int npo = po.n_rows;
    const int npo2 = po2.n_rows;
    arma::mat U1 = a8(U, p, k);
	arma::mat spu(npo2, npo);
	spu.fill(0);
	
    for (int i = 0; i < npo; i++) {
        arma::vec Ts(k);
        Ts.zeros();
			
        if (po(i,0) == 0) {
            Ts = InfU(U1);
        } else {
            Ts = rcpp_standPow(U1, po(i,0));
        }
        for (int j = 0; j < npo2; j++) {
            if (po2(j,0) == 0) {
                arma::vec Ts0 =abs(Ts);
                spu(j,i) = Ts0.max();
            } else {
                spu(j,i) = accu(pow(Ts, po2(j,0)));
            }
        }
    }
    return vectorise(spu.t(), 1);
}






// [[Rcpp::export]]
List perm(arma::mat invR, arma::mat G, arma::mat res, int n, int k, int p, arma::mat po, arma::mat po2, int nperm) {
    
    const int npo = po.n_rows;
    const int npo2 = po2.n_rows;
    arma::mat U = gauss_score(invR, G, res, n, k, p);
    arma::mat spu = rcpp_spuval(U, po, po2, p, k);
    arma::mat T0s(nperm, npo * npo2);
    T0s.fill(0);
	
    for (int i = 0; i < nperm; i++) {
        arma::mat U0 = gauss_score(invR, G, shuffle(res, 0), n, k, p);
        T0s.row(i) = rcpp_spuval(U0, po, po2, p, k);
    }
    arma::vec Pspu(npo * npo2); 
    arma::vec minp0(nperm);
    for (int j = 0; j < T0s.n_cols; j++) {
      Pspu(j) = (1+count_if(std::abs(spu(j)) < abs(T0s.col(j))))/(nperm+1);
      arma::vec P0s = (nperm + 1 - rankC(as<NumericVector>(wrap(abs(T0s.col(j))))))/nperm;
        if (j==0){
            minp0 = P0s;
        } else {
            minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
        }
	}
    double PaSPU = (1+count_if(minp0 < Pspu.min()))/(nperm+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    return(A) ;
}






// [[Rcpp::export]]
List perm_score(arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp, arma::mat po, arma::mat po2, int nperm) {
	
	const int npo = po.n_rows;
	const int npo2 = po2.n_rows;
	
	List S = gauss_score_cov(invR, G, x, res, covres, n, k, p, tp);
	arma::mat U = as<arma::mat>(S["score_vec"]);
	arma::mat sigma = as<arma::mat>(S["score_cov"]);
	arma::mat psigma = arma::pinv(sigma);
	double score = as<double>(S["score"]);
	double PaSPU;
	arma::mat spu_score(1, npo * npo2 + 1);
	spu_score.submat(0, 0, 0, npo * npo2 - 1)= rcpp_spuval(U, po, po2, p, k);
	spu_score(0, npo * npo2)= score;
	arma::mat T0s(nperm, npo * npo2+1);
	T0s.fill(0);
	
	for (int i = 0; i < nperm; i++) {
	    arma::mat U0 = gauss_score(invR, G, shuffle(res, 0), n, k, p);
	    T0s.submat(i, 0, i, npo * npo2 - 1) = rcpp_spuval(U0, po, po2, p, k);
	    T0s(i, npo * npo2) = as<double>(wrap(U0 * psigma *  U0.t()));
	}
	
    arma::vec Pspu(npo * npo2+1); 
    arma::vec minp0(nperm);
    for (int j = 0; j < npo * npo2+1; j++) {
        Pspu(j) = (1+count_if(std::abs(spu_score(0,j)) < abs(T0s.col(j))))/(nperm+1);
    	arma::vec P0s = (nperm +1 - rankC(as<NumericVector>(wrap(abs(T0s.col(j))))))/nperm;
        if (j==0){
            minp0 = P0s;
        } else {
            minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
            if (j==npo * npo2-1){
            	PaSPU = (1+count_if(minp0 < min(Pspu.subvec(0, npo * npo2-1))))/(nperm+1);
            }
        }
	}
    double PaSPU_score = (1+count_if(minp0 < Pspu.min()))/(nperm+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    A["P_aSPU.score"] = PaSPU_score;
    return(A) ;
}











// [[Rcpp::export]]
List permhigh(arma::mat invR, arma::mat G, arma::mat res, int n, int k, int p, arma::mat po, arma::mat po2, int nperm) {
    
  const int npo = po.n_rows;
  const int npo2 = po2.n_rows;
  arma::mat U = gauss_score(invR, G, res, n, k, p);
  arma::mat spu = rcpp_spuval(U, po, po2, p, k);
  arma::vec Pspu(npo*npo2);
  arma::vec minp0(nperm);
  minp0.fill(0);
  int temp = 0;

    for (int q = 0; q < npo; q++){
        for (int q2 = 0; q2 < npo2; q2++){	
            arma::vec T0s(nperm);
            T0s.fill(0);

            //arma::arma_rng::set_seed(1234); 
            for (int i = 0; i < nperm; i++) {
                arma::vec Ts(k);
                arma::mat U0 = gauss_score(invR, G, shuffle(res, 0), n, k, p);
                arma::mat U1 = a8(U0, p, k);
                if (po(q,0) == 0) {
                    Ts = InfU(U1);
                } else {
                    Ts = rcpp_standPow(U1, po(q,0));
                }
                if (po2(q2,0) == 0) {
                    T0s(i) = abs(Ts).max();
                } else {
                    T0s(i) = accu(pow(Ts, po2(q2,0))); 
                }
            }
            Pspu(temp) = (1+count_if(std::abs(spu(temp)) < abs(T0s)))/(nperm+1);
            arma::vec P0s = (nperm + 1 - rankC(as<NumericVector>(wrap(abs(T0s)))))/nperm;
            if (temp==0){
                minp0 = P0s;
            } else {
                minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
            }
            temp++  ;
        }
    }
    double PaSPU = (1+count_if(minp0 < Pspu.min()))/(nperm+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    return(A) ;
}



// [[Rcpp::export]]
List permhigh_score(arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp, arma::mat po, arma::mat po2, int nperm) {
	
	const int npo = po.n_rows;
	const int npo2 = po2.n_rows;
	
	List S = gauss_score_cov(invR, G, x, res, covres, n, k, p, tp);
	arma::mat U = as<arma::mat>(S["score_vec"]);
	arma::mat sigma = as<arma::mat>(S["score_cov"]);
	arma::mat psigma = arma::pinv(sigma);
	double score = as<double>(S["score"]);
	arma::mat spu_score(1, npo * npo2 + 1);
	spu_score.submat(0, 0, 0, npo * npo2 - 1)= rcpp_spuval(U, po, po2, p, k);
	spu_score(0, npo * npo2)= score;
    arma::vec Pspu(npo*npo2 + 1);
    arma::vec minp0(nperm);
    int temp = 0;

    for (int q = 0; q < npo; q++){
        double w = po(q,0);
        for (int q2 = 0; q2 < npo2; q2++){	
            double w2 = po2(q2,0);
            arma::vec T0s(nperm);
            T0s.fill(0);

            //arma::arma_rng::set_seed(1234); 
            for (int i = 0; i < nperm; i++) {
                arma::vec Ts(k);
                arma::mat U0 = gauss_score(invR, G, shuffle(res, 0), n, k, p);
                arma::mat U1 = a8(U0, p, k);
                if (w == 0) {
                    Ts = InfU(U1);
                } else {
                    Ts = rcpp_standPow(U1, w);
                }
                if (w2 == 0) {
                    T0s(i) = max(abs(Ts));
                } else {
                    T0s(i) = accu(pow(Ts, w2)); 
                }
            }
            Pspu(temp) = (1+count_if(std::abs(spu_score(temp)) < abs(T0s)))/(nperm+1);
            arma::vec P0s = (nperm + 1 - rankC(as<NumericVector>(wrap(abs(T0s)))))/nperm;
            if (temp==0){
                minp0 = P0s;
            } else {
                minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
            }
            temp++  ;  
        }
    }
    double PaSPU = (1+count_if(minp0 < min(Pspu.subvec(0, npo * npo2-1))))/(nperm+1);
    arma::vec T0s(nperm);
    T0s.fill(0);
    //arma::arma_rng::set_seed(1234); 
    for (int i = 0; i < nperm; i++) {
	    arma::mat U0 = gauss_score(invR, G, shuffle(res, 0), n, k, p);
	    T0s(i) =  as<double>(wrap(U0 * psigma *  U0.t()));
    }
	
	Pspu(temp) = (1+count_if(std::abs(spu_score(temp)) < abs(T0s)))/(nperm+1);
	arma::vec P0s = (nperm +1- rankC(as<NumericVector>(wrap(abs(T0s)))))/nperm;
	minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
	
    double PaSPU_score = (1+count_if(minp0 < Pspu.min()))/(nperm+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    A["P_aSPU.score"] = PaSPU_score;
    return(A) ;
}


// [[Rcpp::export]]
List sim(int f, arma::vec va, arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp, arma::mat po, arma::mat po2, int nsim) {

    const int npo = po.n_rows;
    const int  npo2 = po2.n_rows;
    List S;
    if (f == 0){
        S = gauss_score_cov(invR, G, x, res, covres, n, k, p, tp);
    } else {
        S = bin_score_cov(va, invR, G, x, res, covres, n, k, p, tp);
    }

	arma::mat U = as<arma::mat>(S["score_vec"]);
	arma::mat sigma = as<arma::mat>(S["score_cov"]);
	List E = getEigen(sigma);
	arma::mat eigvec = as<arma::mat>(E["eigvec"]);
	NumericVector e = as<NumericVector>(E["eigval"]);
	arma::vec eig = as<arma::vec>(which_min(e, e<0));

	arma::mat spu = rcpp_spuval(U, po, po2, p, k);
	arma::mat T0s(nsim, npo * npo2);
	T0s.fill(0);
	for (int i=0; i<nsim; i++) {
        arma::mat Y = arma::randn(1, p*k);
        arma::mat U0 = eigvec * diagmat(eig) *Y.t();
        T0s.row(i) = rcpp_spuval(U0, po, po2, p, k);
    }
	 
    arma::vec Pspu(npo * npo2); 
    arma::vec minp0(nsim);
    for (int j = 0; j < T0s.n_cols; j++) {
        Pspu(j) = (1+count_if(std::abs(spu(j)) < abs(T0s.col(j))))/(nsim+1);
    	arma::vec P0s = (nsim +1- rankC(as<NumericVector>(wrap(abs(T0s.col(j))))))/nsim;
        if (j==0){
            minp0 = P0s;
        } else {
            minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
        }
	}
    double PaSPU = (1+count_if(minp0 < Pspu.min()))/(nsim+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    return(A) ;
}



// [[Rcpp::export]]
List sim_score(int f, arma::vec va, arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp, arma::mat po, arma::mat po2, int nsim) {
    
    const int  npo = po.n_rows;
    const int  npo2 = po2.n_rows;
    List S;
    if (f == 0){
        S = gauss_score_cov(invR, G, x, res, covres, n, k, p, tp);
    } else {
        S = bin_score_cov(va, invR, G, x, res, covres, n, k, p, tp);
    }
    arma::mat U = as<arma::mat>(S["score_vec"]);
    arma::mat sigma = as<arma::mat>(S["score_cov"]);
    arma::mat psigma = arma::pinv(sigma);
    double score = as<double>(S["score"]);
    double PaSPU;
	List E = getEigen(sigma);
	arma::mat eigvec = as<arma::mat>(E["eigvec"]);
	NumericVector e = as<NumericVector>(E["eigval"]);
	arma::vec eig = as<arma::vec>(which_min(e, e<0));
	
    arma::mat spu_score(1, npo * npo2+1);
    spu_score.submat(0, 0, 0, npo * npo2 - 1)= rcpp_spuval(U, po, po2, p, k);
    spu_score(0, npo * npo2)= score;
	 
    arma::mat T0s(nsim, npo * npo2 + 1);
    T0s.fill(0);
    arma::vec mu(p*k);

    for (int i = 0; i < nsim; i++) {
        arma::mat Y = arma::randn(1, p*k);
        arma::mat U0 = eigvec * diagmat(eig) *Y.t();
        T0s.submat(i, 0, i, npo * npo2 - 1) = rcpp_spuval(U0, po, po2, p, k);
        T0s(i, npo * npo2) = arma::as_scalar(U0.t() * psigma * U0);
    }
	
    arma::vec Pspu(npo * npo2+1); 
    arma::vec minp0(nsim);
    for (int j = 0; j < npo * npo2+1; j++) {
        Pspu(j) = (1+count_if(std::abs(spu_score(0,j)) < abs(T0s.col(j))))/(nsim+1);
    	arma::vec P0s = (nsim +1- rankC(as<NumericVector>(wrap(abs(T0s.col(j))))))/nsim;
        if (j==0){
            minp0 = P0s;
        } else {
            minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
            if (j==npo * npo2-1){
            	PaSPU = (1+count_if(minp0 < min(Pspu.subvec(0, npo * npo2-1))))/(nsim+1);
            }
        }
	}
    double PaSPU_score = (1+count_if(minp0 < Pspu.min()))/(nsim+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    A["P_aSPU.score"] = PaSPU_score;
    return(A) ;
}





// [[Rcpp::export]]
List simhigh(int f, arma::vec va, arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp, arma::mat po, arma::mat po2, int nperm) {
    
    const int npo = po.n_rows;
    const int npo2 = po2.n_rows;
    List S;
    if (f == 0){
        S = gauss_score_cov(invR, G, x, res, covres, n, k, p, tp);
    } else {
        S = bin_score_cov(va, invR, G, x, res, covres, n, k, p, tp);
    }
	arma::mat U = as<arma::mat>(S["score_vec"]);
	arma::mat sigma = as<arma::mat>(S["score_cov"]);
	List E = getEigen(sigma);
	arma::mat eigvec = as<arma::mat>(E["eigvec"]);
	NumericVector e = as<NumericVector>(E["eigval"]);
	arma::vec eig = as<arma::vec>(which_min(e, e<0));
	arma::mat spu = rcpp_spuval(U, po, po2, p, k);

    arma::vec Pspu(npo*npo2);
    arma::vec minp0(nperm);
    int temp = 0;

    for (int q = 0; q < npo; q++){
        double w = po(q,0);
        for (int q2 = 0; q2 < npo2; q2++){	
            double w2 = po2(q2,0);
            arma::vec T0s(nperm);
            T0s.fill(0);

            //arma::arma_rng::set_seed(1234); 
            for (int i = 0; i < nperm; i++) {
                arma::vec Ts(k);
                arma::mat Y = arma::randn(1, p*k);
                arma::mat U0 = eigvec * diagmat(eig) *Y.t();
                arma::mat U1 = a8(U0, p, k);
                if (w == 0) {
                    Ts = InfU(U1);
                } else {
                    Ts = rcpp_standPow(U1, w);
                }
                if (w2 == 0) {
                    T0s(i) = max(abs(Ts));
                } else {
                    T0s(i) = accu(pow(Ts, w2)); 
                }
            }
            Pspu(temp) = (1+count_if(std::abs(spu(temp)) < abs(T0s)))/(nperm+1);
            arma::vec P0s = (nperm +1- rankC(as<NumericVector>(wrap(abs(T0s)))))/nperm;
            if (temp==0){
                minp0 = P0s;
            } else {
                minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
            }
            temp++  ;  
        }
    }
    double PaSPU = (1+count_if(minp0 < Pspu.min()))/(nperm+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    return(A) ;
}



// [[Rcpp::export]]
List simhigh_score(int f, arma::vec va, arma::mat invR, arma::mat G, arma::mat x, arma::mat res, arma::mat covres, int n, int k, int p, int tp, arma::mat po, arma::mat po2, int nperm) {
    
    const int  npo = po.n_rows;
    const int  npo2 = po2.n_rows;
    List S;
    if (f == 0){
        S = gauss_score_cov(invR, G, x, res, covres, n, k, p, tp);
    } else {
        S = bin_score_cov(va, invR, G, x, res, covres, n, k, p, tp);
    }
    arma::mat U = as<arma::mat>(S["score_vec"]);
    arma::mat sigma = as<arma::mat>(S["score_cov"]);
    arma::mat psigma = arma::pinv(sigma);
    double score = as<double>(S["score"]);
	List E = getEigen(sigma);
	arma::mat eigvec = as<arma::mat>(E["eigvec"]);
	NumericVector e = as<NumericVector>(E["eigval"]);
	arma::vec eig = as<arma::vec>(which_min(e, e<0));
	 
    arma::mat spu_score(1, npo * npo2+1);
    spu_score.submat(0, 0, 0, npo * npo2 - 1)= rcpp_spuval(U, po, po2, p, k);
    spu_score(0, npo * npo2)= score;
    arma::vec Pspu(npo*npo2 + 1);
    arma::vec minp0(nperm);
    int temp = 0;

    for (int q = 0; q < npo; q++){
        double w = po(q,0);
        for (int q2 = 0; q2 < npo2; q2++){	
            double w2 = po2(q2,0);
            arma::vec T0s(nperm);
            T0s.fill(0);

            //arma::arma_rng::set_seed(1234); 
  
            for (int i = 0; i < nperm; i++) {
                arma::vec Ts(k);
                arma::mat Y = arma::randn(1, p*k);
                arma::mat U0 = eigvec * diagmat(eig) *Y.t();
                arma::mat U1 = a8(U0, p, k);
                if (w == 0) {
                    Ts = InfU(U1);
                } else {
                    Ts = rcpp_standPow(U1, w);
                }
                if (w2 == 0) {
                    T0s(i) = max(abs(Ts));
                } else {
                    T0s(i) = accu(pow(Ts, w2)); 
                }
            }
            Pspu(temp) = (1+count_if(std::abs(spu_score(temp)) < abs(T0s)))/(nperm+1);
            arma::vec P0s = (nperm +1- rankC(as<NumericVector>(wrap(abs(T0s)))))/nperm;
            if (temp==0){
                minp0 = P0s;
            } else {
                minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
            }
            temp++  ;  
        }
    }
    double PaSPU = (1+count_if(minp0 < min(Pspu.subvec(0, npo * npo2-1))))/(nperm+1);
    arma::vec T0s(nperm);
    T0s.fill(0);
    //arma::arma_rng::set_seed(1234); 
    for (int i = 0; i < nperm; i++) {
        arma::mat Y = arma::randn(1, p*k);
        arma::mat U0 = eigvec * diagmat(eig) *Y.t();
	    T0s(i) =  as<double>(wrap(U0.t() * psigma *  U0));
    }

	Pspu(temp) = (1+count_if(std::abs(spu_score(temp)) < abs(T0s)))/(nperm+1);
	arma::vec P0s = (nperm +1- rankC(as<NumericVector>(wrap(abs(T0s)))))/nperm;
	minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
	
    double PaSPU_score = (1 + count_if(minp0 < Pspu.min()))/(nperm+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPU"] = PaSPU;
    A["P_aSPU.score"] = PaSPU_score;
    return(A) ;
}







// [[Rcpp::export]]
arma::mat pathval(arma::mat out, int k, arma::vec nSNPs, int nGenes, arma::mat po, arma::mat po2, arma::mat po3){
    
    const int npo = po.n_rows;
    const int npo2 = po2.n_rows;
    const int npo3 = po3.n_rows;
    int lp = npo * npo2;

	arma::mat spu(k, lp);
	spu.fill(0);

	for (int x = 0; x < k; x++) {
        arma::vec U = out.row(x).t();
        arma::mat O(npo2, npo);
        for (int i = 0; i < npo; i++){
            arma::vec Ts(nGenes);
            Ts.fill(0);
            for (int iGene = 0; iGene < nGenes; iGene++){
            	int SNPstart;
            	if (iGene==0) {
            		SNPstart = 0; 
            	} else {
            		SNPstart = accu(nSNPs.subvec(0, iGene - 1));
            	}
            	int SNPend = SNPstart + nSNPs(iGene) - 1; 
				if (po(i,0) == 0) {
            		arma::mat U01 = abs(U.subvec(SNPstart, SNPend));
        			Ts(iGene) = U01.max();
        		} else {
        			double a = accu(pow(U.subvec(SNPstart, SNPend), po(i,0)));
       				Ts(iGene) = signC(a) * std::pow(std::abs(a)/nSNPs(iGene), 1.0/po(i,0));
        		}
        	}
        	for (int j = 0; j < npo2; j++){
            	if (po2(j,0) == 0) {
            		arma::mat w = abs(Ts);
        			O(j,i) = w.max();
        		} else {
        			double a1 = accu(pow(Ts, po2(j,0)));
       				O(j,i) = signC(a1) * std::pow(std::abs(a1), 1.0/po2(j,0));
        		}
            }
        }
        spu.row(x) = vectorise(O.t(), 1);	
    }           		
    arma::mat path(npo3, lp); 
    for (int g = 0; g < npo3; g++){
    	for (int x2 = 0; x2 < lp; x2++){
			if (po3(g,0) == 0){
				path(g,x2) = max(abs(spu.col(x2)));
			} else {
				path(g,x2) = accu(pow(spu.col(x2), po3(g,0)));
			}
    	}
     }
   return(vectorise(path.t(), 1));    
}   
 

arma::mat pathscore(arma::mat invR, arma::mat res, arma::mat G, int n, int k, int p){
 	arma::mat out(k,p);
 	out.fill(0);
	for (int i = 0; i < n; i++) {
    	out += (res.row(i) *invR).t() * G.row(i);
	}
	return(out);
}

// [[Rcpp::export]]
List permpath(arma::mat invR, arma::mat G, arma::mat res, arma::vec nSNPs, int nGenes, int n, int k, int p, arma::mat po, arma::mat po2, arma::mat po3, int nperm) {
    
    const int npo = po.n_rows;
    const int npo2 = po2.n_rows;
    const int npo3 = po3.n_rows;
    arma::mat out = pathscore(invR, res, G, n,k,p);
	arma::mat spu = pathval(out, k, nSNPs, nGenes, po, po2, po3);
    arma::mat T0s(nperm, npo * npo2 * npo3);
 	T0s.fill(0);
	
    for (int i = 0; i < nperm; i++) {
        arma::mat out0 = pathscore(invR, shuffle(res, 0), G, n, k, p);
        T0s.row(i) = pathval(out0, k, nSNPs, nGenes, po, po2, po3);
    }
    arma::vec Pspu(npo * npo2* npo3); 
    arma::vec minp0(nperm);
    for (int j = 0; j < T0s.n_cols; j++) {
        Pspu(j) = (1+count_if(std::abs(spu(j)) < abs(T0s.col(j))))/(nperm+1);
    	arma::vec P0s = (nperm +1- rankC(as<NumericVector>(wrap(abs(T0s.col(j))))))/nperm;
        if (j==0){
            minp0 = P0s;
        } else {
            minp0(find(minp0 > P0s)) = P0s(find(minp0 > P0s));
        }
	}
    double PaSPUpath = (1+count_if(minp0 < Pspu.min()))/(nperm+1);
    List A ;
    A["P_SPU"] = Pspu;
    A["P_aSPUpath"] = PaSPUpath;
    return(A) ;
}



