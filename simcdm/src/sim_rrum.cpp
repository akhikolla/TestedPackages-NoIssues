#include <RcppArmadillo.h>

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat sim_rrum_main(const arma::mat &Q, const arma::mat &rstar,
                        const arma::vec &pistar, const arma::mat &alpha)
{

    unsigned int N = alpha.n_rows;
    unsigned int J = Q.n_rows;
    unsigned int K = Q.n_cols;

    arma::vec k_index = arma::linspace(0, K - 1, K);
    double kj;
    double aik;
    arma::vec pmu(J);
    arma::mat Y = arma::zeros<arma::mat>(N, J);
    for (unsigned int i = 0; i < N; i++) {
        arma::vec Yi = arma::zeros<arma::vec>(J);
        arma::vec pi = arma::ones<arma::vec>(J);
        arma::vec ui = arma::randu<arma::vec>(J);
        for (unsigned int j = 0; j < J; ++j) {
            arma::uvec task_ij = find(Q.row(j) == 1);

            for (unsigned int k = 0; k < task_ij.n_elem; ++k) {
                kj = task_ij(k);
                aik = alpha(i, kj);
                pi(j) = ((rstar(j, kj) * (1.0 - aik) + 1.0 * aik) * Q(j, kj) +
                         1.0 * (1.0 - Q(j, kj))) *
                        pi(j);
            }
            pi(j) = pistar(j) * pi(j);
        }
        pmu = pi - ui;
        Yi(arma::find(pmu > 0)).fill(1);
        Y.row(i) = Yi.t();
    }
    return Y;
}

//' Generate data from the rRUM
//'
//' Randomly generate response data according to the reduced Reparameterized
//' Unified Model (rRUM).
//'
//' @param Q      A `matrix` with \eqn{J} rows and \eqn{K} columns indicating
//'               which attributes are required to answer each of the items,
//'               where \eqn{J} represents the number of items and \eqn{K} the
//'               number of attributes.  An entry of 1 indicates attribute
//'               \eqn{k} is required to answer item \eqn{j}.  An entry of one
//'               indicates attribute \eqn{k} is not required.
//' @param rstar  A `matrix` a matrix with \eqn{J} rows and \eqn{K} columns
//'               indicating the penalties for failing to have each of the
//'               required attributes, where \eqn{J} represents the number of
//'               items and \eqn{K} the number of attributes. `rstar` and `Q`
//'               must share the same 0 entries.
//' @param pistar A `vector` of length \eqn{J} indicating the probabiliies of
//'               answering each item correctly for individuals who do not lack
//'               any required attribute, where \eqn{J} represents the number
//'               of items.
//' @param alpha  A `matrix` with \eqn{N} rows and \eqn{K} columns indicating
//'               the subjects attribute acquisition, where \eqn{K} represents
//'               the number of attributes.  An entry of 1 indicates individual
//'               \eqn{i} has attained attribute \eqn{k}. An entry of 0
//'               indicates the attribute has not been attained.
//' @return Y     A `matrix` with \eqn{N} rows and \eqn{J} columns indicating
//'               the indviduals' responses to each of the items, where \eqn{J}
//'               represents the number of items.
//' @author
//' Steven Andrew Culpepper, Aaron Hudson, and James Joseph Balamuta
//'
//' @export
//' @template rrum-example
//' @template rrum-references
// [[Rcpp::export]]
arma::mat sim_rrum_items(const arma::mat &Q, const arma::mat &rstar,
                         const arma::vec &pistar, const arma::mat &alpha)
{

    if (Q.n_rows != rstar.n_rows || Q.n_cols != rstar.n_cols) {
        Rcpp::stop("Q and rstar must have the same dimensionaltiy");
    } else if (!arma::approx_equal(find(Q == 0), find(rstar == 0), "absdiff",
                                   0)) {
        Rcpp::stop("rstar and Q must have the same 0 entries");
    } else if (pistar.n_elem != Q.n_rows) {
        Rcpp::stop("length(pistar) must be equal to nrow(Q)");
    } else if (alpha.n_cols != Q.n_cols) {
        Rcpp::stop("alpha and Q must have the same number of columns");
    }

    bool rstar_check =
        arma::all(arma::vectorise(arma::all(rstar >= 0 && rstar <= 1) == true));

    bool pistar_check = arma::all(pistar >= 0 && pistar <= 1);

    if (!(rstar_check && pistar_check)) {
        Rcpp::stop(
            "all entries of pistar and rstar must be bound between 0 and 1");
    }

    return sim_rrum_main(Q, rstar, pistar, alpha);
}
