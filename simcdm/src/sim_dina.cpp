#include <RcppArmadillo.h>

// [[Rcpp::interfaces(r, cpp)]]

//' Simulate Binary Responses for a DINA Model
//'
//' Generate the dichotomous item matrix for a DINA Model.
//'
//' @param N     Number of Observations
//' @param J     Number of Assessment Items
//' @param CLASS Does the individual possess all the necessary attributes?
//' @param ETA   \eqn{\eta} Matrix containing indicators.
//' @param gs    A `vec` describing the probability of guessing or
//'              the probability subject correctly answers item \eqn{j} when at
//'              least one attribute is lacking.
//' @param ss    A `vec` describing the probability of slipping or
//'              the probability of an incorrect response for individuals with
//'              all of the required attributes
//'
//' @return
//' A dichotomous item matrix with dimensions \eqn{N \times J}{N x J}.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [simcdm::sim_dina_attributes()] and [simcdm::sim_dina_items()]
//'
//' @export
//' @template sim-dina-class-example
// [[Rcpp::export]]
arma::mat sim_dina_class(unsigned int N, unsigned int J, const arma::vec &CLASS,
                         const arma::mat &ETA, const arma::vec &gs,
                         const arma::vec &ss)
{
    arma::mat Y(N, J);
    for (unsigned int i = 0; i < N; ++i) {
        double class_i = CLASS(i);
        arma::vec ETA_i = ETA.col(class_i);
        for (unsigned int j = 0; j < J; ++j) {
            double u = R::runif(0, 1);
            Y(i, j) =
                1. * (gs(j) * (1. - ETA_i(j)) + (1. - ss(j)) * ETA_i(j) > u);
        }
    }
    return Y;
}

//' Simulate a DINA Model's \eqn{\eta} Matrix
//'
//' Generates a DINA model's \eqn{\eta} matrix based on alphas and
//' the \eqn{\mathbf{Q}} matrix.
//'
//' @inheritParams sim_dina_items
//'
//' @return
//' The \eqn{\eta} `matrix` with dimensions \eqn{N \times J}{N x J} under
//' the DINA model.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [simcdm::sim_dina_class()] and [simcdm::sim_dina_items()]
//'
//' @export
//' @template sim-dina-example-body
// [[Rcpp::export]]
arma::mat sim_dina_attributes(const arma::mat &alphas, const arma::mat &Q)
{
    unsigned int N = alphas.n_rows;
    unsigned int J = Q.n_rows;

    arma::mat ETA(N, J);

    for (unsigned int j = 0; j < J; ++j) {
        for (unsigned int i = 0; i < N; ++i) {
            if (arma::dot(alphas.row(i), Q.row(j)) <
                arma::dot(Q.row(j), Q.row(j))) {
                ETA(i, j) = 0.0;
            } else {
                ETA(i, j) = 1.0;
            }
        }
    }

    return ETA;
}

//' Simulation Responses from the DINA model
//'
//' Sample responses from the DINA model for given attribute profiles, Q matrix,
//' and item parmeters. Returns a `matrix` of dichotomous responses
//' generated under DINA model.
//'
//' @param alphas A \eqn{N} by \eqn{K} `matrix` of latent attributes.
//' @param Q      A \eqn{J} by \eqn{K} `matrix` indicating which skills are
//'               required for which items.
//' @param ss     A \eqn{J} `vector` of item slipping parameters. 
//' @param gs     A \eqn{J} `vector` of item guessing parameters.
//'
//' @return
//' A \eqn{N} by \eqn{J} `matrix` of responses from the DINA model.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [simcdm::sim_dina_class()] and [simcdm::sim_dina_attributes()]
//'
//' @export
//' @template sim-dina-example-body
// [[Rcpp::export]]
arma::mat sim_dina_items(const arma::mat &alphas, const arma::mat &Q,
                         const arma::vec &ss, const arma::vec &gs)
{
    unsigned int N = alphas.n_rows;
    unsigned int J = Q.n_rows;

    arma::mat Y = arma::zeros<arma::mat>(N, J);
    double eta_ij;
    double uij;

    for (unsigned int j = 0; j < J; ++j) {
        for (unsigned int i = 0; i < N; ++i) {
            uij = R::runif(0., 1.);

            if (arma::dot(alphas.row(i), Q.row(j)) <
                arma::dot(Q.row(j), Q.row(j))) {
                eta_ij = 0.0;
            } else {
                eta_ij = 1.0;
            }

            if (pow(1.0 - ss(j), eta_ij) * pow(gs(j), 1.0 - eta_ij) > uij) {
                Y(i, j) = 1;
            } else {
                Y(i, j) = 0;
            }
        }
    }

    return Y;
}
