// load Rcpp
#include <RcppArmadillo.h>

using namespace Rcpp;       // shorthand

// [[Rcpp::export]]
List gibbs(NumericMatrix x, int ite, double a = 2.1, double b = 1.1, double gamma_a = 1,
    double gamma_b = 1, double omega_0 = 0.01, double omega_1 = 10, bool degenerate = false) {

    int m = x.nrow();
    int n = x.ncol();

    arma::colvec alpha(rep(0.1, m));
    arma::rowvec sigma2(rep(1.0, m));
    arma::rowvec z(rep(0.0, m));
    arma::rowvec p_star(rep(0.5, m));
    arma::rowvec lambda(rnorm(n, 0, 1));

    arma::mat dinv(m, m);
    dinv.fill(0.0);

    arma::mat alpha_matrix(ite, m);
    arma::mat lambda_matrix(ite, n);
    arma::mat sigma_matrix(ite, m);
    arma::mat z_matrix(ite, m);
    arma::mat p_matrix(ite, m);

    arma::mat x_arma(x.begin(), m, n, false);

    for (int k = 0; k < ite; k++) {
        for (int i = 0; i < m; i++) {
            sigma2[i] = 1/R::rgamma(a + n/2, 1/(b + 0.5*(arma::as_scalar(x_arma.row(i)*x_arma.row(i).t()) -
                2*arma::as_scalar(lambda*x_arma.row(i).t())*alpha[i] + 
                alpha[i]*arma::as_scalar(lambda*lambda.t())*alpha[i])));

            if(R::runif(0,1) < p_star[i]) {
                z[i] = 1;
            } else {
                z[i] = 0;
            }

            double p = R::rbeta(gamma_a + z[i], gamma_b + 1 - z[i]);

            double lambda_ss = 0.0;
            for(int j = 0; j<n; j++) {
                lambda_ss += pow(lambda[j], 2);
            }
            double v_alpha_0 = 1/(1/omega_0 + lambda_ss/sigma2[i]);
            double lambda_x_sum = 0.0;
            for(int j = 0; j<n; j++) {
                lambda_x_sum += x_arma(i,j)*lambda[j];
            }
            double m_alpha_0 = v_alpha_0*(lambda_x_sum/sigma2[i]);
            
            double v_alpha_1 = 1/(1/omega_1 + lambda_ss/sigma2[i]);
            double m_alpha_1 = v_alpha_1*(lambda_x_sum/sigma2[i]);

            if(z[i]) {
                alpha[i] = R::rnorm(m_alpha_1, sqrt(v_alpha_1));
            } else {
                if(degenerate) {
                    alpha[i] = 0;
                } else {
                    alpha[i] = R::rnorm(m_alpha_0, sqrt(v_alpha_0));
                }
            }
            
            if(degenerate) {
                p_star[i] = p/(p + exp(R::dnorm(0, m_alpha_1, sqrt(v_alpha_1), true) -
                    R::dnorm(0, 0, sqrt(omega_1), true))*(1-p));
            } else {
                p_star[i] = p/(p + exp(R::dnorm(0, m_alpha_1, sqrt(v_alpha_1), true) -
                    R::dnorm(0, 0, sqrt(omega_1), true) + R::dnorm(0, 0, sqrt(omega_0), true) - 
                    R::dnorm(0, m_alpha_0, sqrt(v_alpha_0), true))*(1-p));
            }

            dinv(i,i) = 1/sigma2[i];

        }

        for(int j = 0; j<n; j++) {
            double v_lambda = 1/(arma::as_scalar(alpha.t()*dinv*alpha + 1));
            lambda[j] = R::rnorm(v_lambda*(arma::as_scalar(alpha.t()*dinv*x_arma.col(j))),
                sqrt(v_lambda));
        }

        alpha_matrix.row(k) = alpha.t();
        lambda_matrix.row(k) = lambda;
        sigma_matrix.row(k) = sigma2;
        z_matrix.row(k) = z;
        p_matrix.row(k) = p_star;

    }
    return Rcpp::List::create(        
        Rcpp::Named("alpha", wrap(alpha_matrix)),
        Rcpp::Named("lambda", wrap(lambda_matrix)),
        Rcpp::Named("sigma", wrap(sigma_matrix)),
        Rcpp::Named("p_star", wrap(p_matrix)),
        Rcpp::Named("z_matrix", wrap(z_matrix))
        );             // Return to R
}
