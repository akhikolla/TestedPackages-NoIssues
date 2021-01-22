#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec bijectionvector(unsigned int K)
{
    arma::vec vv(K);
    for (unsigned int k = 0; k < K; k++) {
        vv(k) = std::pow(2.0, static_cast<double>(K - k) - 1.0);
    }
    return vv;
}

// [[Rcpp::export]]
arma::vec inv_bijectionvector(unsigned int K, double CL)
{
    arma::vec alpha(K);
    for (unsigned int k = 0; k < K; k++) {
        double twopow = std::pow(2.0, static_cast<double>(K - k) - 1.0);
        alpha(k) = (twopow <= CL);
        CL = CL - twopow * alpha(k);
    }
    return alpha;
}

// [[Rcpp::export]]
arma::mat CL_invbijection_table(unsigned int K, unsigned int nClass)
{
    arma::mat CLtable(K, nClass);
    for (unsigned int cc = 0; cc < nClass; cc++) {
        CLtable.col(cc) = inv_bijectionvector(K, cc);
    }
    return CLtable;
}

//' Generate Multinomial Random Variable
//'
//' Sample a multinomial random variable for given probabilities.
//'
//' @param ps A \code{vector} for the probability of each category.
//'
//' @return
//' A \code{vector} from a multinomial with probability ps.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @noRd
// [[Rcpp::export]]
double rmultinomial(const arma::vec &ps)
{
    unsigned int C = ps.n_elem;
    double u = R::runif(0, 1);
    arma::vec cps = cumsum(ps);
    arma::vec Ips = arma::zeros<arma::vec>(C);
    Ips.elem(arma::find(cps < u)).fill(1.0);

    return sum(Ips);
}

//' Generate Dirichlet Random Variable
//'
//' Sample a Dirichlet random variable.
//'
//' @param deltas A \code{vector} of Dirichlet parameters.
//'
//' @return
//' A \code{vector} from a Dirichlet.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @noRd
// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec &deltas)
{
    unsigned int C = deltas.n_elem;
    arma::vec Xgamma(C);

    // generating gamma(deltac,1)
    for (unsigned int c = 0; c < C; c++) {
        Xgamma(c) = R::rgamma(deltas(c), 1.0);
    }
    return Xgamma / sum(Xgamma);
}

// [[Rcpp::export]]
arma::mat random_Q(unsigned int J, unsigned int K)
{
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    arma::vec vv = bijectionvector(K);
    arma::vec Q_biject(J);
    Q_biject(arma::span(0, K - 1)) = vv;
    Q_biject(arma::span(K, 2 * K - 1)) = vv;
    Q_biject(arma::span(2 * K, 3 * K - 1)) = vv;
    arma::vec Jm3K =
        arma::randi<arma::vec>(J - 3 * K, arma::distr_param(1, nClass - 1));
    Q_biject(arma::span(3 * K, J - 1)) = Jm3K;
    Q_biject = arma::shuffle(Q_biject);
    arma::mat Q(J, K);
    for (unsigned int j = 0; j < J; j++) {
        arma::vec qj = inv_bijectionvector(K, Q_biject(j));
        Q.row(j) = qj.t();
    }
    return Q;
}

// [[Rcpp::export]]
arma::mat simrRUM(unsigned int N, unsigned int J, unsigned int K,
                  const arma::mat &Q, const arma::mat &rstar,
                  const arma::vec &pistar, const arma::vec &CLASS)
{

    arma::vec k_index = arma::linspace(0, K - 1, K);
    double kj;
    double aik;
    arma::vec pmu(J);
    arma::mat Y = arma::zeros<arma::mat>(N, J);
    for (unsigned int i = 0; i < N; i++) {
        arma::vec Yi = arma::zeros<arma::vec>(J);
        arma::vec pi = arma::ones<arma::vec>(J);
        arma::vec ui = arma::randu<arma::vec>(J);
        arma::vec alpha_i = inv_bijectionvector(K, CLASS(i));
        for (unsigned int j = 0; j < J; j++) {
            arma::uvec task_ij = find(Q.row(j) == 1);

            for (unsigned int k = 0; k < task_ij.n_elem; k++) {
                kj = task_ij(k);
                aik = alpha_i(kj);
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

// rRUM log-likelihood
// [[Rcpp::export]]
double DEVIANCErRUM(unsigned int N, unsigned int J, unsigned int K,
                    unsigned int nClass, const arma::mat &Y,
                    const arma::mat &rstar, const arma::vec &pistar,
                    const arma::vec &pis)
{

    // create a table of probabilities first
    arma::mat PYeq1(J, nClass);
    for (unsigned int j = 0; j < J; j++) {
        arma::rowvec rj = rstar.row(j);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            arma::vec alpha_c = inv_bijectionvector(K, cc);
            double ptemp = 1.;
            for (unsigned int k = 0; k < K; k++) {
                double ak = alpha_c(k);
                ptemp *= (rj(k) * (1. - ak) + ak);
            }
            PYeq1(j, cc) = ptemp * pistar(j);
        }
    }
    // use table to compute marginal likelihood
    double loglike = 0.;
    for (unsigned int i = 0; i < N; i++) {
        arma::rowvec Yi = Y.row(i);
        arma::vec pitemp(nClass);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            double pctemp = 1.;
            arma::vec PYeq1cc = PYeq1.col(cc);
            for (unsigned int j = 0; j < J; j++) {
                double pyj1 = PYeq1cc(j);
                double yij = Yi(j);
                pctemp *= pyj1 * yij + (1. - pyj1) * (1. - yij);
            }
            pitemp(cc) = pctemp;
        }
        double marglik_i = arma::accu(pitemp % pis);
        loglike += log(marglik_i);
    }
    return -2. * loglike;
}

// [[Rcpp::export]]
arma::mat simgnida(unsigned int N, unsigned int J, unsigned int K,
                   const arma::mat &Q, const arma::mat &B0, const arma::mat &B1,
                   const arma::vec &CLASS)
{

    arma::vec k_index = arma::linspace(0, K - 1, K);
    double kj;
    double aik;
    arma::vec pmu(J);
    arma::mat Y = arma::zeros<arma::mat>(N, J);
    for (unsigned int i = 0; i < N; i++) {
        arma::vec Yi = arma::zeros<arma::vec>(J);
        arma::vec pi = arma::ones<arma::vec>(J);
        arma::vec ui = arma::randu<arma::vec>(J);
        arma::vec alpha_i = inv_bijectionvector(K, CLASS(i));
        for (unsigned int j = 0; j < J; j++) {
            arma::uvec task_ij = find(Q.row(j) == 1);

            for (unsigned int k = 0; k < task_ij.n_elem; k++) {
                kj = task_ij(k);
                aik = alpha_i(kj);
                pi(j) *= R::pnorm(B0(j, k) + B1(j, k) * aik, 0., 1., 1, 0);
            }
        }
        pmu = pi - ui;
        Yi(arma::find(pmu > 0)).fill(1);
        Y.row(i) = Yi.t();
    }
    return Y;
}

// [[Rcpp::export]]
double rTruncNorm(double mean, double sd, double w)
{
    double p1 = R::pnorm(mean, 0, sd, 1, 0);
    double p0 = 1 - p1;
    double uZ = R::runif(0, 1);
    double pz = w * (p0 + uZ * p1) + (1 - w) * (uZ * p0);
    double Z = R::qnorm(pz, mean, sd, 1, 0);
    return (Z);
}

// this version has no Q, but updates Bs by conditioning as before

// [[Rcpp::export]]
void parm_update_nomiss(unsigned int N, unsigned int J, unsigned int K,
                        unsigned int nClass, const arma::mat &Y,
                        arma::mat &X_biject, arma::mat &B0, arma::mat &B1,
                        arma::vec &CLASS, arma::vec &pis, arma::mat &Q,
                        const arma::mat &MU, const arma::vec &vv,
                        const arma::mat &CLtable, double v0, double v1,
                        const arma::vec &d0, arma::vec &pistar,
                        arma::mat &rstar, double a0, double b0, double a1,
                        double b1, arma::mat &Qstar, double cv1, double cv0)
{

    arma::vec k_index = arma::linspace(0, K - 1, K);
    double aik, u, c_aik_1, c_aik_0, prodXijk, Xij_biject, Yij, pi_ijk, Zijk;

    arma::mat sumZ = arma::zeros<arma::mat>(J, K);
    arma::mat sumZa = arma::zeros<arma::mat>(J, K);
    arma::vec dtilde = arma::zeros<arma::vec>(nClass);
    arma::vec suma = arma::zeros<arma::vec>(K);
    arma::vec B0k_B1k(K);
    arma::vec B1k_B1k_dividedby2(K);

    for (unsigned int k = 0; k < K; k++) {
        B0k_B1k(k) = arma::accu(B0.col(k) % B1.col(k));
        B1k_B1k_dividedby2(k) = arma::accu(B1.col(k) % B1.col(k)) / 2.;
    }
    // computing probabilities for X given alpha
    arma::cube PX_a(J, K, 2);
    for (unsigned int j = 0; j < J; j++) {
        for (unsigned int k = 0; k < K; k++) {
            PX_a(j, k, 0) = R::pnorm(B0(j, k), .0, 1., 1, 0);
            PX_a(j, k, 1) = R::pnorm(B0(j, k) + B1(j, k), .0, 1., 1, 0);
        }
    }

    for (unsigned int i = 0; i < N; i++) {
        arma::vec ai = CLtable.col(CLASS(i));
        arma::rowvec Yi = Y.row(i);
        arma::rowvec Xi_biject = X_biject.row(i);
        arma::mat Zi(J, K);

        for (unsigned int j = 0; j < J; j++) {

            Yij = Yi(j);
            Xij_biject = Xi_biject(j);
            arma::vec Xij = CLtable.col(Xij_biject);

            for (unsigned int k = 0; k < K; k++) {
                u = R::runif(0., 1.);
                aik = ai(k);
                Xij(k) = 1.;
                // prodXijk = prod(Xij(task_ij));
                prodXijk = prod(Xij);
                pi_ijk = (1. - prodXijk) * PX_a(j, k, aik);
                double compare = 1. * (pi_ijk > u);
                double Xijk = (1. - Yij) * compare + Yij;
                Xij(k) = Xijk;
                Zijk = rTruncNorm(B0(j, k) + B1(j, k) * aik, 1., Xijk);
                // update Zijk
                Zi(j, k) = Zijk;
                // update sumZ for sampling Bs later
                sumZ(j, k) += Zijk;
            }
            Xi_biject(j) = arma::accu(Xij % vv);
        }
        X_biject.row(i) = Xi_biject;
        // update alphas conditioned on Zi
        for (unsigned int k = 0; k < K; k++) {
            ai(k) = 1.0;
            c_aik_1 = arma::accu(ai % vv);
            ai(k) = 0.0;
            c_aik_0 = arma::accu(ai % vv);
            u = R::runif(0., 1.);
            double lnA = arma::accu((Zi.col(k)) % B1.col(k)) - B0k_B1k(k) -
                         B1k_B1k_dividedby2(k) + log(pis(c_aik_1));
            double lnB = log(pis(c_aik_0));
            aik = 1. * (log(1. - u) - log(u) > lnB - lnA);
            ai(k) = aik;
            suma(k) += aik;
            if (aik == 1) {
                sumZa.col(k) += Zi.col(k);
            }
        }
        double ai_biject = arma::accu(ai % vv);
        CLASS(i) = ai_biject;
        dtilde(ai_biject) += 1.;
    }
    // update pis
    pis = rDirichlet(dtilde + d0);

    // update Q and B
    double sigma02 = 1. / (double(N) + v0);
    double sigma02null = 1. / (double(N) + cv0);
    arma::vec sigma12(K);
    arma::vec sigma12null(K);
    for (unsigned int k = 0; k < K; k++) {
        sigma12(k) = 1. / (suma(k) + v1);
        sigma12null(k) = 1. / (suma(k) + cv1);
    }
    for (unsigned int j = 0; j < J; j++) {
        double pistar_temp = 1.0;

        for (unsigned int k = 0; k < K; k++) {
            // update qjk
            u = R::runif(0, 1);
            double pqjk = R::pnorm(MU(j, k), .0, 1., 1, 0);
            double B0jk = B0(j, k);
            double B1jk = B1(j, k);

            u = R::runif(0, 1);
            double pdjk = pqjk; // omega;

            double lnA = .5 * log(v1) - .5 * v1 * B1jk * B1jk + .5 * log(v0) -
                         .5 * v0 * B0jk * B0jk + log(pdjk);
            double lnB = .5 * log(cv1) - .5 * cv1 * B1jk * B1jk +
                         .5 * log(cv0) - .5 * cv0 * B0jk * B0jk +
                         log(1. - pdjk);
            double qjk = 1. * (log(1. - u) - log(u) > lnB - lnA);
            Q(j, k) = qjk;
            // update Qstar
            Qstar(j, k) = rTruncNorm(MU(j, k), 1., qjk);

            if (qjk == 1.) {
                // update B0jk
                double mu0jk = sigma02 * (sumZ(j, k) - B1(j, k) * suma(k));
                B0jk = R::rnorm(mu0jk, sqrt(sigma02));
                double mu1jk = sigma12(k) * (sumZa(j, k) - B0jk * suma(k));
                B1jk = rTruncNorm(mu1jk, sqrt(sigma12(k)), 1.);
            } else {
                double mu0jk = sigma02null * (sumZ(j, k) - B1(j, k) * suma(k));
                B0jk = R::rnorm(mu0jk, sqrt(sigma02null));
                double mu1jk = sigma12null(k) * (sumZa(j, k) - B0jk * suma(k));
                B1jk = rTruncNorm(mu1jk, sqrt(sigma12null(k)), 1.);
                // sumB1qsq+=B1jk*B1jk;
                // sumD+=1.;
            }
            double p0 = R::pnorm(B0jk, .0, 1., 1, 0);
            double p1 = R::pnorm(B0jk + B1jk, .0, 1., 1, 0);
            pistar_temp *= p1; // compute pistarj
            rstar(j, k) = p0 / p1;

            B0(j, k) = B0jk;
            B1(j, k) = B1jk;
        }
        pistar(j) = pistar_temp;
    }
}

// assume 1st predictor is intercept
// [[Rcpp::export]]
Rcpp::List update_Gamma_Delta_MVN(unsigned int N, unsigned int V,
                                  unsigned int K, unsigned int J,
                                  const arma::mat &Y, arma::mat &Gamma,
                                  arma::mat &deltas, const arma::mat &X,
                                  const arma::vec &diagXpX, double w2,
                                  double nu, double bnu)
{

    arma::vec onetoV = arma::linspace(0, V - 1, V);
    arma::vec onetoK = arma::linspace(0, K - 1, K);
    double aforw2 = .0;
    double bforw2 = .0;
    // double afornu = .0;
    // update deltas and gammas across variables
    // update intercepts first
    arma::mat Ytilde(J, K);
    if (V == 1) {
        Ytilde = Y;
    } else {
        Ytilde = Y - X.cols(1, V - 1) * (Gamma.rows(1, V - 1));
    }
    // arma::mat Ytilde = Y - X.cols(1,V-1)*(Gamma.rows(1,V-1));
    double sigv2 = 1. / (double(J) + 1.);
    arma::rowvec mu_v = sigv2 * (arma::sum(Ytilde));
    for (unsigned int k = 0; k < K; ++k) {
        Gamma(0, k) = R::rnorm(mu_v(k), sqrt(sigv2));
    }
    if (V > 1) {
        for (unsigned int vv = 1; vv < V; ++vv) {
            arma::rowvec delta_v = deltas.row(vv);
            arma::uvec notv = find(onetoV != vv);
            arma::vec Xv = X.col(vv);
            arma::mat Ytilde = Y - X.cols(notv) * (Gamma.rows(notv));
            double sigv2 = 1. / (diagXpX(vv) + w2);
            arma::vec mu_v = sigv2 * (Ytilde.t() * Xv);
            arma::vec lnb_m_lna = log(1. - nu) - .5 * pow(mu_v, 2) / sigv2 -
                                  log(nu) - .5 * (log(w2) + log(sigv2));
            for (unsigned int k = 0; k < K; ++k) {
                arma::uvec notk = find(onetoK != k);
                double uvk = R::runif(0, 1);
                double lnuvk = log(uvk);
                double ln1muvk = log(1. - uvk);
                // double dvk = (1.-(arma::accu(delta_v(notk))>0.)
                // )*(ln1muvk-lnuvk>lnb_m_lna(k));//imposes structured sparsity
                double dvk = (ln1muvk - lnuvk > lnb_m_lna(k));
                if (dvk == 1) {
                    double bvktemp = R::rnorm(mu_v(k), sqrt(sigv2));
                    bforw2 += pow(bvktemp, 2);
                    aforw2 += 1.;
                    Gamma(vv, k) = bvktemp;
                } else {
                    Gamma(vv, k) = 0.;
                }
                deltas(vv, k) = dvk;
            }
        }
    }

    // update w2
    double uw = R::runif(0, 1);
    // qgamma is parameterized as shape, scale rather than shape and rate
    double w2new =
        R::qgamma(uw, aforw2 / 2.0 + 1., 1.0 / (.5 * bforw2 + 1.), 1, 0);
    // update nu
    double u_nu = R::runif(0, 1);
    double nu_new =
        R::qbeta(u_nu, aforw2 + 1., (V - 1.) * K - aforw2 + bnu, 1, 0);

    return Rcpp::List::create(Rcpp::Named("w2new", w2new),
                              Rcpp::Named("nu_new", nu_new));
}

// [[Rcpp::export]]
arma::mat OddsRatio(unsigned int N, unsigned int J, const arma::mat &Yt)
{
    arma::mat M2_temp = arma::zeros<arma::mat>(J, J);
    for (unsigned int j1 = 0; j1 < J - 1; ++j1) {
        for (unsigned int j2 = j1 + 1; j2 < J; ++j2) {
            double n11 = arma::accu(Yt.col(j1) % Yt.col(j2));
            double n00 = arma::accu((1. - Yt.col(j1)) % (1. - Yt.col(j2)));
            double n10 = arma::accu(Yt.col(j1) % (1. - Yt.col(j2)));
            double n01 = N - n11 - n00 - n10;
            M2_temp(j1, j2) = (n11 * n00) / (n10 * n01);
        }
    }
    return M2_temp;
}

// [[Rcpp::export]]
arma::mat uppertri_matrix_logical_gt(unsigned int J, const arma::mat &A,
                                     const arma::mat &B)
{
    arma::mat AA = arma::zeros<arma::mat>(J, J);
    for (unsigned int j1 = 0; j1 < J - 1; ++j1) {
        for (unsigned int j2 = j1 + 1; j2 < J; ++j2) {
            AA(j1, j2) = 1. * (A(j1, j2) > B(j1, j2));
        }
    }
    return AA;
}

// [[Rcpp::export]]
Rcpp::List rRUM_mvnQ_Gibbs(const arma::mat &Y, unsigned int K,
                           const arma::mat &X, double v0, double v1, double cv0,
                           double cv1, double bnu, unsigned int burnin = 1000,
                           unsigned int chain_length = 10000,
                           bool verbose = false)
{

    // Parameter initialization
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    // Chain
    unsigned int chain_p_burn = chain_length + burnin;
    unsigned int tmburn;

    // Bijection
    arma::vec vv = bijectionvector(K);

    // Odds Ratio
    arma::mat Observed_ORs = OddsRatio(N, J, Y);

    // X Matrix operations
    arma::mat XpX = X.t() * X;
    arma::vec diagXpX = XpX.diag();
    unsigned int V = XpX.n_cols;

    // Q matrix construction
    arma::mat Q =
        arma::randi<arma::mat>(J, K, arma::distr_param(0, 1)); // random_Q(J,K);
    arma::mat CLtable = CL_invbijection_table(K, nClass);

    arma::mat Z = arma::randn<arma::mat>(J, K);
    arma::mat Gamma = arma::randn<arma::mat>(V, K);
    arma::mat deltas = arma::zeros<arma::mat>(
        V, K); // arma::randi<arma::vec>(V,arma::distr_param(0,K));
    arma::vec oneK = arma::ones<arma::vec>(K);
    double w2(1.), nu(1.);

    // Save output
    arma::vec deviance(chain_length);
    arma::mat PISTAR(J, chain_length);
    arma::cube RSTAR(J, K, chain_length);
    arma::mat PIs(nClass, chain_length);
    // arma::vec cv1s(chain_length);
    // arma::vec v1s(chain_length);
    // arma::vec v0s(chain_length);
    arma::cube QS(J, K, chain_length);
    arma::cube GAMMAS_OUT(V, K, chain_length);
    arma::mat Delta_tab = arma::zeros<arma::mat>(V, K);
    arma::cube B0s(J, K, chain_length);
    arma::cube B1s(J, K, chain_length);
    arma::mat OR_PPPs = arma::zeros<arma::mat>(J, J);

    // need to initialize parameters
    arma::vec CLASS =
        arma::randi<arma::vec>(N, arma::distr_param(0, nClass - 1));
    arma::mat B0 = arma::randn<arma::mat>(J, K);
    arma::mat B1 = 2. * arma::randu<arma::mat>(J, K);
    arma::mat MU = X * Gamma;
    arma::mat Qstar(J, K);
    // double v0(.5),v1(2.);//,omega(.5);
    double a0(1.), b0(1.), a1(1.), b1(1.);

    arma::vec d0 = arma::ones<arma::vec>(nClass);
    arma::vec pis = rDirichlet(d0);
    arma::mat rstar = arma::ones<arma::mat>(J, K);
    arma::vec pistar = arma::randu<arma::vec>(J);
    arma::mat X_biject = arma::ones<arma::mat>(N, J);

    // Start Markov chain
    for (unsigned int t = 0; t < chain_p_burn; ++t) {
        parm_update_nomiss(N, J, K, nClass, Y, X_biject, B0, B1, CLASS, pis, Q,
                           MU, vv, CLtable, v0, v1, d0, pistar, rstar, a0, b0,
                           a1, b1, Qstar, cv1, cv0);
        // Rcpp::List step_1
        // =parm_update_nomiss(N,J,K,nClass,Y,X_biject,B0,B1,CLASS,pis,Q,MU,vv,CLtable,v0,v1,d0,pistar,rstar,a0,b0,a1,b1,Qstar,cv1);
        // v0 = Rcpp::as< double >(step_1[0]);
        // v1 = Rcpp::as< double >(step_1[1]);
        // omega = Rcpp::as< double >(step_1[2]);

        arma::mat MU = X * Gamma;

        // Update Gamma & deltas
        Rcpp::List step_MVN = update_Gamma_Delta_MVN(
            N, V, K, J, Qstar, Gamma, deltas, X, diagXpX, w2, nu, bnu);
        w2 = Rcpp::as<double>(step_MVN[0]);
        nu = Rcpp::as<double>(step_MVN[1]);

        if (t > burnin - 1) {
            tmburn = t - burnin;
            PISTAR.col(tmburn) = pistar;
            RSTAR.slice(tmburn) = rstar;
            PIs.col(tmburn) = pis;
            QS.slice(tmburn) = Q;
            GAMMAS_OUT.slice(tmburn) = Gamma;
            Delta_tab += deltas;
            /*      B0s.slice(tmburn)     = B0;
             B1s.slice(tmburn)     = B1;
             arma::mat Yt = simrRUM(N,J,K,Q,rstar,pistar,CLASS);
             arma::mat OR_r = OddsRatio(N, J, Yt);
             arma::mat logic_r =
             uppertri_matrix_logical_gt(J,OR_r,Observed_ORs);
             OR_PPPs=1./(tmburn+1.)*logic_r+tmburn/(tmburn+1.)*OR_PPPs;
             deviance(tmburn)=DEVIANCErRUM(N,J,K,nClass,Y,rstar,pistar,pis);
             */
        }

        if (verbose && t % 100 == 0) {
            Rcpp::Rcerr << "On iteration " << t << " out of " << chain_p_burn
                        << "." << std::endl;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("PISTAR", PISTAR), Rcpp::Named("RSTAR", RSTAR),
        Rcpp::Named("PIs", PIs), Rcpp::Named("QS", QS),
        Rcpp::Named("GAMMAS", GAMMAS_OUT), Rcpp::Named("DELTAS", Delta_tab)
        // Rcpp::Named("B0s",B0s),
        // Rcpp::Named("B1s",B1s),
        // Rcpp::Named("deviance",deviance),
        // Rcpp::Named("OR_PPPs",OR_PPPs)
    );
}
