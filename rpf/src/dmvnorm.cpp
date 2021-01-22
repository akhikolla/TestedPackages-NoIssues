#include "rpf.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>

template<typename _MatrixType, int _UpLo = Eigen::Lower>
class SimpCholesky : public Eigen::LDLT<_MatrixType, _UpLo> {
 private:
        Eigen::MatrixXd inverse;

 public:
        typedef Eigen::LDLT<_MatrixType, _UpLo> Base;

        SimpCholesky() : Base() {};
        template<typename InputType>
        explicit SimpCholesky(const Eigen::EigenBase<InputType>& matrix) : Base(matrix) {};
        template<typename InputType>
        explicit SimpCholesky(Eigen::EigenBase<InputType>& matrix) : Base(matrix) {};

        double log_determinant() const {
                typename Base::Scalar detL = Base::vectorD().array().log().sum();
                return detL * 0.5;
        }

        void refreshInverse()
        {
                inverse.setIdentity(Base::m_matrix.rows(), Base::m_matrix.rows());
                inverse = Base::solve(inverse);
        };

        const Eigen::MatrixXd &getInverse() const { return inverse; };
};

static const int ERROR_LEN = 80;

static double
_mahalanobis(char *err, int dim, double *loc, double *center, double *origCov)
{
	Eigen::VectorXd cloc(dim);
	for (int dx=0; dx < dim; dx++) {
		cloc[dx] = loc[dx] - center[dx];
	}

	Eigen::Map<Eigen::MatrixXd> covMat(origCov, dim, dim);
	SimpCholesky< Eigen::MatrixXd, Eigen::Lower > sc(covMat);
	if (sc.info() != Eigen::Success || !sc.isPositive()) {
		snprintf(err, ERROR_LEN, "Sigma is singular and cannot be inverted");
		return nan("mxThrow");
	}

	sc.refreshInverse();
	return cloc.transpose() * sc.getInverse() * cloc;
}

static double
mahalanobis(int dim, double *loc, double *center, double *origCov)
{
	char err[ERROR_LEN];
	err[0] = 0;
	double ret = _mahalanobis(err, dim, loc, center, origCov);
	if (err[0]) stop("%s", err);
	return ret;
}

static double
_dmvnorm(char *err, int dim, double *loc, double *mean, double *origSigma)
{
	Eigen::Map< Eigen::MatrixXd > Esigma(origSigma, dim, dim);
	SimpCholesky< Eigen::MatrixXd, Eigen::Upper > sc(Esigma);

	double dist = mahalanobis(dim, loc, mean, origSigma);

	double got = -(dim * M_LN_SQRT_2PI + 0.5*dist + sc.log_determinant());
	return got;
}

double
dmvnorm(int dim, double *loc, double *mean, double *sigma)
{
	char err[ERROR_LEN];
	err[0] = 0;
	double ret = _dmvnorm(err, dim, loc, mean, sigma);
	if (err[0]) stop("%s", err);
	return ret;
}
