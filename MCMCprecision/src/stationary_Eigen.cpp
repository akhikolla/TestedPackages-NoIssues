// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
#include <progress.hpp>

using namespace Rcpp;

// using Eigen::Map;               	// 'maps' rather than copies
// using Eigen::MatrixXd;                  // variable size matrix, double precision
// using Eigen::VectorXd;                  // variable size vector, double precision
// using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers


Eigen::MatrixXd rdirichletPt(Eigen::MatrixXd Pt)
{
  double colsum;
  for (int j=0; j<Pt.cols(); j++)
  {
    for (int i=0; i<Pt.cols(); i++)
    {
      if (Pt(i,j) != 0)
        Pt(i,j) = R::rgamma(Pt(i,j), 1);
    }
    colsum = Pt.col(j).sum();
    if (colsum > 0)
      Pt.col(j) /= colsum;
  }
  return (Pt);
}


// Posterior distribution for stationary distribution
// N: matrix with transition frequencies
// a: prior vector for transition probabilities - Dirichlet(a[1],...,a[M])
// sample: number of (independent) posterior samples
// [[Rcpp::export]]
Eigen::MatrixXd stationaryEigen(Eigen::MatrixXd N,
                                double epsilon = 0,
                                int sample = 5000,
                                bool progress=true,
                                double digits = 8.)
{
  int M = N.cols();
  int steps = round (1000/M);
  Eigen::MatrixXd samp = -99 * Eigen::MatrixXd::Ones(sample, M);
  Eigen::MatrixXd freqt = N.transpose() + epsilon * Eigen::MatrixXd::Ones(M, M);
  Progress p(sample, progress);

  Eigen::MatrixXf::Index maxIdx;
  bool run = true;
  for (int i=0; i < sample; i++)
  {
    p.increment();   // update progress bar
    if (run && i % steps == 0)
      run = !Progress::check_abort(); // check whether to abort

    if (run)
    {
      // 1.) sample from conjugate posterior: Dirichlet
      Eigen::MatrixXd Pt = rdirichletPt(freqt);
      // 2.) get estimate for stationary distribution
      //     (normalized left eigenvector for eigenvalue = 1)
      try
      {
        Eigen::EigenSolver<Eigen::MatrixXd> es(Pt);
        Eigen::VectorXcd values = es.eigenvalues();
        values.real().maxCoeff(&maxIdx);
        // Rcout  << maxIdx << "///" << values << "\n";

        if( abs(values(maxIdx).real() - 1) < pow(10, -digits) )
        {
          Eigen::VectorXd ev = es.eigenvectors().col(maxIdx).real();
          samp.row(i) = ev / ev.sum();
        }
      }
      catch(...)
      {
        Rcout << "# RcppEigen::SelfAdjointEigenSolver unstable: \n#" <<
          "method='base' or changing the argument N.min > 1 might provide more stable results#";
      }
    }
  }
  return samp;
}
