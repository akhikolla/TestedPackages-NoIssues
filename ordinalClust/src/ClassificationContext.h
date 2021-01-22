#pragma once
#include "Distribution.h"
#include "Bos.h"

//#include "SparsePoisson.h"

// [[Rcpp::depends(RcppArmadillo)]] 
//#include <armadillo>
#include <limits>
#include <cmath>
#include <list>
#include <typeinfo>
#include <iostream>
#include <initializer_list> 
#include <vector>
#include <numeric>

// [[Rcpp::depends(BH)]]
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>

using namespace std;

class ClassificationContext
{
public:
	ClassificationContext(arma::mat x, arma::vec y, std::vector< arma::urowvec > dlist,
		int kr, std::vector< int > kc, std::string init, int nbSEM, int nbSEMburn, 
		int nbindmini, std::vector< int > m, std::vector<double> percentRandomB, int seed);
	ClassificationContext();
	~ClassificationContext();
	void missingValuesInit(); 
	void initialization();
	void Mstep();
	void MstepVW();
	void SEstep();
	void sampleVW();
	void sampleVWStock();
	void imputeMissingData();
	bool verif();
	vector<vector<int>> verification();
	void noColDegenerancy(vector<vector<int>> distrib_col, int iter);
	void fillParameters(int iteration);
	void fillLabels(int iteration);
	void getBurnedParameters();
	void printResults();
	void returnResults();
	void  putParamsToZero();
	double computeICL();
	// TODO:
	// find logprobas for topredict
	// find probas for topredict
	// instanciate V_topredict
	// instanciate zr_topredict
	S4 returnClassification();
	
protected:
	mat _x;
	int _Nr;
	vector<int> _Jc;
	vector<int>_m;
	vector<urowvec> _dlist;
	vector<Distribution*> _distrib_objects;
	int _number_distrib;
	int _kr;
	vector<int> _kc;
	vector<int> _zr;
	vector<vector<int>> _zc;
	vector<mat> _probaW;
	vector<mat> _logprobaW;
	mat _V;
	vector<mat> _W;
	vector<rowvec> _rho;
	rowvec _gamma;
	vector<vector<rowvec>> _allrho;
	vector<rowvec> _allgamma;
	vector<rowvec> _resrho;
	rowvec _resgamma;
	// params of co-clustering :
	string _init;
	int _nbSEM;
	int _nbSEMburn;
	int _nbindmini;
	// params for classificationContext
	vec _y;

	std::vector<mat> _zcchain;

	// to sample distributions
	random_device _rd;

	double _icl;
	int _seed;
	vector<double> _percentRandomB;

	// Utils
	rowvec getMeans(mat VorW);
	double logsum(rowvec logx);
	mat kmeansi();
	double getDistance(vec &a, vec &b);
};

