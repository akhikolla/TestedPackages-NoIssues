#pragma once
#include "Distribution.h"
#define _USE_MATH_DEFINES
#include <math.h>

// [[Rcpp::depends(BH)]]
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>


using namespace std;
using namespace arma;
using namespace Rcpp;
class BosPredict
{

public:
	BosPredict(int kr, int kc, int m, mat pis, mat mus, int seed);
	BosPredict();
	~BosPredict();
	mat missingValuesInit(mat& x);
	mat SEstep_predict(mat W, mat x);
	cube getCubeProbs();
	cube gettabpej();
protected:
	int _m;
	mat _pis;
	mat _mus;
	int _kr;
	int _kc;
	vector<vector<int>> _miss;
	cube _tab_pejs;
	random_device _rd; // to sample distributions
	int _seed;
};

