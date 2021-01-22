// load Rcpp
#include <RcppArmadillo.h>

using namespace Rcpp;       // shorthand

// [[Rcpp::export]]
List slfm_MNN(NumericMatrix x, double a = 2.1, double b = 1.1, double gamma_a = 1, double gamma_b = 1, 
           double omega_0 = 0.01, double omega_1 = 10, int burnin = 1000, int lag = 1, int npost = 500) {
    // Number of iterations
    int ite = burnin+npost*lag;

    // X matrix dimensions
    int m = x.nrow();
    int n = x.ncol();

    // vectors of parameters for the Gibbs Sampler
    arma::colvec alpha(rep(0.0, m));
    arma::rowvec sigma2(rep(1.0, m));
    arma::rowvec z( rbinom(m,1,0.5) );
    arma::rowvec qstar(rep(0.5, m));
    arma::rowvec q(rep(0.5, m));
    arma::rowvec lambda( rnorm(n, 0, sqrt(0.1)) );

    arma::mat dinv(m,m);
    dinv.fill(0.0);
    for(int i = 0; i < m; i++){ dinv(i,i) = 1/sigma2[i]; }
  
    // matrices to store the chains
    arma::mat alpha_matrix(ite, m);
    arma::mat lambda_matrix(ite, n);
    arma::mat sigma2_matrix(ite, m);
    arma::mat z_matrix(ite, m);
    arma::mat qstar_matrix(ite, m);

    // converting the x matrix to the arma class
    arma::mat x_arma(x.begin(), m, n, false);

    // iterations loop
    for(int k = 0; k < ite; k++) {

        // j-dimensional parameter loop
        double lambda_ss = 0.0;
        for(int j = 0; j<n; j++) {
           // sampling lambda
           double aux = arma::as_scalar(alpha.t()*dinv*alpha) + 1;
           double v_lambda = 1/aux;
           aux = arma::as_scalar(alpha.t()*dinv*x_arma.col(j));
           double m_lambda = v_lambda*aux;
           lambda[j] = R::rnorm(m_lambda,sqrt(v_lambda));
           lambda_ss += pow(lambda[j], 2);
        }    

        // i-dimensional parameters loop
        for(int i = 0; i < m; i++) {
            // sampling sigma2
            double A = a + n/2;
            double B = arma::as_scalar(x_arma.row(i)*x_arma.row(i).t());
                   B = B-2*arma::as_scalar(lambda*x_arma.row(i).t())*alpha[i];
                   B = b+0.5*(B+alpha[i]*arma::as_scalar(lambda*lambda.t())*alpha[i]);
            sigma2[i] = 1/R::rgamma(A,1/B);
            
            dinv(i,i) = 1/sigma2[i];

            // sampling z
            if(R::runif(0,1) < qstar[i]) {
                z[i] = 1;
            } else {
                z[i] = 0;
            }

            // sampling q
            q[i] = R::rbeta(gamma_a + z[i], gamma_b + 1 - z[i]);

            // sampling alpha
            double lambda_x_sum = 0.0;
            for(int j = 0; j<n; j++) { lambda_x_sum += x_arma(i,j)*lambda[j]; }
            
            double v_alpha_1 = 1/(1/omega_1 + lambda_ss/sigma2[i]);
            double m_alpha_1 = v_alpha_1*(lambda_x_sum/sigma2[i]);
 
            double v_alpha_0 = 1/(1/omega_0 + lambda_ss/sigma2[i]);
            double m_alpha_0 = v_alpha_0*(lambda_x_sum/sigma2[i]);

            if(z[i]) {
                alpha[i] = R::rnorm(m_alpha_1, sqrt(v_alpha_1));
            } else {
                alpha[i] = R::rnorm(m_alpha_0, sqrt(v_alpha_0));
            }

            // qstar normalization
            A = R::dnorm(0, 0, sqrt(omega_0), true);
            A = A - R::dnorm(0, m_alpha_0, sqrt(v_alpha_0), true);
            A = A + R::dnorm(0, m_alpha_1, sqrt(v_alpha_1), true);
            A = A - R::dnorm(0, 0, sqrt(omega_1), true);
            A = q[i] + exp(A)*(1-q[i]);
            qstar[i] = q[i]/A;
 
        }

        // storing results in matrix
        alpha_matrix.row(k) = alpha.t();
        lambda_matrix.row(k) = lambda;
        sigma2_matrix.row(k) = sigma2;
        qstar_matrix.row(k) = qstar;
        z_matrix.row(k) = z;

    }
    // returning matrices
    return Rcpp::List::create(        
        Rcpp::Named("alpha", wrap(alpha_matrix)),
        Rcpp::Named("lambda", wrap(lambda_matrix)),
        Rcpp::Named("sigma2", wrap(sigma2_matrix)),
        Rcpp::Named("qstar", wrap(qstar_matrix)),
        Rcpp::Named("z", wrap(z_matrix))
        );             // Return to R
}
