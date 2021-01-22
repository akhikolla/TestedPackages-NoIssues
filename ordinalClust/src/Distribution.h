#pragma once
#include "LogProbs.h"
#include "TabProbsResults.h"
#include <iostream>


// [[Rcpp::depends(RcppArmadillo)]] 
//#include <armadillo>
#include <limits>
//#include <cmath>
#include <list>
#include <random>
#include <vector>

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;


class Distribution
{
public:
	explicit Distribution(mat xsep, int kr, int kc, int nbSEM, int seed);
	Distribution();
	virtual ~Distribution();
	virtual void missingValuesInit();
	virtual TabProbsResults SEstep(mat V, mat W);
	virtual mat SEstepRow(mat V, mat W);
	virtual mat SEstepCol(mat V, mat W);
	virtual LogProbs SEstep_predict(int i, int d, int k, int h, double x_id, double sumi, double sumd, vec x_id_vec);
	virtual void imputeMissingData(mat V, mat W);
	virtual void Mstep(uvec rowind, uvec colind, int k, int h, bool init);
	virtual void MstepVW(mat V, mat W, bool init);
	virtual void fillParameters(int iteration);
	virtual void getBurnedParameters(int burn);
	virtual void printResults();
	virtual List returnResults();
	virtual List returnParamsChain();
	virtual void putParamsToZero();
	virtual double computeICL(int i, int d, int k, int h);
	bool verif(mat V, mat W, int nbindmini);
	int verification(const mat& V, const mat& W, int nbindmini);
	mat colkmeans();
	double getDistance(vec &a, vec &b);
	mat returnXhat();
protected:
	rowvec getDatablockkh(uvec rowind, uvec colind);
	string _name;
	mat _xsep;
	vector<vector<int>> _miss;
	int _Nr;
	int _Jc;
	int _kr;
	int _kc;
	int _nbSEM;
	int _seed;
	
};

