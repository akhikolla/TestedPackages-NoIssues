#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::interfaces(r, cpp)]]

//' Constructs Unique Attribute Pattern Map
//'
//' Computes the powers of 2 from \eqn{0} up to \eqn{K - 1} for
//' \eqn{K}-dimensional attribute pattern.
//'
//' @param K  Number of Attributes.
//'
//' @return
//' A \code{vec} with length \eqn{K} detailing the power's of 2.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [simcdm::attribute_inv_bijection()]
//'
//' @export
//' @examples
//' ## Construct an attribute bijection ----
//' biject = attribute_bijection(3)
// [[Rcpp::export]]
arma::vec attribute_bijection(unsigned int K)
{
    arma::vec vv(K);
    for (unsigned int i = 0; i < K; ++i) {
        vv(i) = std::pow(2.0, static_cast<double>(K - i) - 1.0);
    }
    return vv;
}

//' Perform an Inverse Bijection of an Integer to Attribute Pattern
//'
//' Convert an integer between \eqn{0} and \eqn{2^{K-1}} to
//' \eqn{K}-dimensional attribute pattern.
//'
//' @param CL An `integer` between \eqn{0} and \eqn{2^{K-1}}
//' @inheritParams attribute_bijection
//'
//' @return
//' A \eqn{K}-dimensional vector with an attribute pattern corresponding
//' to `CL`.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [simcdm::attribute_bijection()]
//'
//' @export
//' @examples
//' ## Construct an attribute inversion bijection ----
//' inv_biject1 = attribute_inv_bijection(5, 1)
//' inv_biject2 = attribute_inv_bijection(5, 2)
// [[Rcpp::export]]
arma::vec attribute_inv_bijection(unsigned int K, double CL)
{
    arma::vec alpha(K);

    for (unsigned int k = 0; k < K; ++k) {
        double twopow = std::pow(2.0, static_cast<double>(K - k) - 1.0);
        alpha(k) = (twopow <= CL);
        CL = CL - twopow * alpha(k);
    }

    return alpha;
}

//' Generate a Random Identifiable Q Matrix
//'
//' Simulates a Q matrix containing three identity matrices after a row
//' permutation that is identifiable.
//'
//' @param J Number of Items
//' @param K Number of Attributes
//'
//' @return
//' A dichotomous \code{matrix} for Q.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [simcdm::attribute_bijection()] and [simcdm::attribute_inv_bijection()]
//'
//' @export
//' @examples
//' ## Simulate identifiable Q matrices ----
//'
//' # 7 items and 2 attributes
//' q_matrix_j7_k2 = sim_q_matrix(7, 2)
//'
//' # 10 items and 3 attributes
//' q_matrix_j10_k3 = sim_q_matrix(10, 3)
// [[Rcpp::export]]
arma::mat sim_q_matrix(unsigned int J, unsigned int K)
{
    // Impose an error check
    if (J < 3 * K - 1) {
        Rcpp::stop("J must be greater than 3*K.");
    }

    // Calculate number of classes
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    // Form a Bijection
    arma::vec vv = attribute_bijection(K);

    // Fill Q bijection vector
    arma::vec Q_biject(J);
    Q_biject(arma::span(0, K - 1)) = vv;
    Q_biject(arma::span(K, 2 * K - 1)) = vv;
    Q_biject(arma::span(2 * K, 3 * K - 1)) = vv;

    // Randomly sample from a discrete uniform
    arma::vec Jm3K =
        arma::randi<arma::vec>(J - 3 * K, arma::distr_param(1, nClass - 1));
    Q_biject(arma::span(3 * K, J - 1)) = Jm3K;

    // Randomize rows
    Q_biject = arma::shuffle(Q_biject);

    // Create Q Matrix
    arma::mat Q(J, K);
    for (unsigned int j = 0; j < J; ++j) {
        arma::vec qj = attribute_inv_bijection(K, Q_biject(j));
        Q.row(j) = qj.t();
    }

    return Q;
}

//' Generate ideal response \eqn{\eta} Matrix
//'
//' Creates the ideal response matrix for each trait
//'
//' @param K      Number of Attribute Levels
//' @param J      Number of Assessment Items
//' @param Q      Q Matrix with dimensions \eqn{K \times J}{K x J}.
//'
//' @return
//' A `mat` with dimensions \eqn{J \times 2^K}{J x 2^K}.
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [simcdm::sim_q_matrix()], [simcdm::attribute_bijection()], and
//' [simcdm::attribute_inv_bijection()]
//'
//' @export
//' @examples
//' ## Simulation Settings ----
//'
//' # Fixed Number of Assessment Items for Q
//' J = 18
//'
//' # Fixed Number of Attributes for Q
//' K = 3
//'
//' ## Pre-specified configuration ----
//'
//' # Specify Q
//' qbj = c(4, 2, 1, 4, 2, 1, 4, 2, 1, 6, 5, 3, 6, 5, 3, 7, 7, 7)
//'
//' # Fill Q Matrix
//' Q = matrix(, J, K)
//' for (j in seq_len(J)) {
//'   Q[j,] = attribute_inv_bijection(K, qbj[j])
//' }
//'
//' # Create an eta matrix
//' ETA = sim_eta_matrix(K, J, Q)
//'
//' ## Random generation of Q matrix with ETA matrix ----
//'
//' # Construct a random q matrix
//' Q_sim = sim_q_matrix(J, K)
//'
//' # Generate the eta matrix
//' ETA_gen = sim_eta_matrix(K, J, Q_sim)
// [[Rcpp::export]]
arma::mat sim_eta_matrix(unsigned int K, unsigned int J, const arma::mat &Q)
{
    // Calculate number of classes
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    arma::mat ETA(J, nClass);

    for (unsigned int cc = 0; cc < nClass; ++cc) {
        arma::vec alpha_c = attribute_inv_bijection(K, cc);

        for (unsigned int j = 0; j < J; ++j) {
            arma::rowvec qj = Q.row(j);
            // Switch to as_scalar
            double compare = arma::as_scalar(qj * alpha_c - qj * qj.t());
            ETA(j, cc) = (compare >= 0);
        }
    }

    return ETA;
}

//' Simulate all the Latent Attribute Profile \eqn{\mathbf{\alpha}_c} in
//' Matrix form
//'
//' Generate the \eqn{\mathbf{\alpha}_c = (\alpha_{c1}, \ldots, \alpha_{cK})'}
//' attribute profile matrix for members of class \eqn{c} such that 
//' \eqn{\alpha_{ck}} ' is 1 if members of class \eqn{c} possess skill \eqn{k}
//' and zero otherwise.
//'
//' @param K Number of Attributes
//'
//' @return
//' A \eqn{2^K} by \eqn{K} `matrix` of latent classes
//' corresponding to entry \eqn{c} of \eqn{pi} based upon
//' mastery and nonmastery of the \eqn{K} skills.
//'
//' @author
//' James Joseph Balamuta and Steven Andrew Culpepper
//'
//' @seealso
//' [simcdm::sim_subject_attributes()] and [simcdm::attribute_inv_bijection()]
//'
//' @export
//' @examples
//' ## Simulate Attribute Class Matrix ----
//'
//' # Define number of attributes
//' K = 3
//'
//' # Generate an Latent Attribute Profile (Alpha) Matrix
//' alphas = attribute_classes(K)
// [[Rcpp::export]]
arma::mat attribute_classes(int K)
{
    // Modified version of ETAMatrix

    // Calculate number of classes
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    // Create alpha matrix
    arma::mat alpha_matrix(nClass, K);

    // Fill alpha matrix with classes under an inverse bijection
    for (unsigned int cc = 0; cc < nClass; cc++) {
        alpha_matrix.row(cc) = attribute_inv_bijection(K, cc).t();
    }

    // Release result
    return alpha_matrix;
}

//' Simulate Subject Latent Attribute Profiles \eqn{\mathbf{\alpha}_c}
//'
//' Generate a sample from the
//'  \eqn{\mathbf{\alpha}_c = (\alpha_{c1}, \ldots, \alpha_{cK})'}
//' attribute profile matrix for members of class \eqn{c} such that
//' \eqn{\alpha_{ck}} ' is 1 if members of class \eqn{c} possess skill \eqn{k} 
//' and zero otherwise.
//'
//' @param N      Number of Observations
//' @param K      Number of Skills
//' @param probs  A `vector` of probabilities that sum to 1.
//'
//' @return
//' A \eqn{N} by \eqn{K} `matrix` of latent classes
//' corresponding to entry \eqn{c} of \eqn{pi} based upon
//' mastery and nonmastery of the \eqn{K} skills.
//'
//' @author
//' James Joseph Balamuta and Steven Andrew Culpepper
//'
//' @seealso
//' [simcdm::attribute_classes()] and [simcdm::attribute_inv_bijection()]
//'
//' @export
//' @examples
//' # Define number of subjects and attributes
//' N = 100
//' K = 3
//'
//' # Generate a sample from the Latent Attribute Profile (Alpha) Matrix
//' # By default, we sample from a uniform distribution weighting of classes.
//' alphas_builtin = sim_subject_attributes(N, K)
//'
//' # Generate a sample using custom probabilities from the
//' # Latent Attribute Profile (Alpha) Matrix
//' probs = rep(1 / (2 ^ K), 2 ^ K)
//' alphas_custom = sim_subject_attributes(N, K, probs)
// [[Rcpp::export]]
arma::mat
sim_subject_attributes(int N, int K,
                       Rcpp::Nullable<Rcpp::NumericVector> probs = R_NilValue)
{
    // Modified version of ETAMatrix

    // Calculate number of classes
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    // Nullable trick ----
    arma::vec probs_;

    // Check if probs is present
    if (probs.isNotNull()) {

        // Retrieve probabilities (costly)
        probs_ = Rcpp::as<arma::vec>(probs);

        // Verify length is okay.
        if (probs_.n_elem != nClass) {
            Rcpp::stop("`probs` must have %s elements instead of %s.", nClass,
                       probs_.size());
        }
    } else {
        probs_.set_size(nClass);
        probs_.fill(1.0 / nClass);
    }

    // --- Profile Matrix

    // Grab the attribute matrix
    arma::mat attributes = attribute_classes(K);

    // Generate indices
    arma::uvec idx = arma::linspace<arma::uvec>(0, nClass - 1, nClass);

    // Update index
    idx = Rcpp::RcppArmadillo::sample_main(idx, N, true, probs_);

    // Retrieve and return profiles
    return attributes.rows(idx);
}
