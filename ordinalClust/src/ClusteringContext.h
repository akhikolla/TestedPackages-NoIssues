#pragma once
#include "Distribution.h"
#include "Bos.h"

//#include "SparsePoisson.h"

// [[Rcpp::depends(RcppArmadillo)]] 
#include <armadillo>
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


class ClusteringContext
{
public:
	ClusteringContext(arma::mat x, std::vector< arma::urowvec > dlist, 
					int kr, std::string init, int nbSEM, int nbSEMburn, 
					int nbindmini, std::vector< int > m, 
					std::vector<double> percentRandomB, int seed);
	ClusteringContext();
	~ClusteringContext();
	void missingValuesInit();
	void initialization();
	void SEstep();
	void SEstepRow();
	void Mstep();
	void MstepVW();
	void noRowDegenerancy(vector<vector<int>> distrib_col, int iter);
	void sampleVW();
	void sampleVWStock();
	void imputeMissingData();
	bool verif();
	vector<vector<int>> verification();
	void fillParameters(int iteration);
	void fillLabels(int iteration);
	void getBurnedParameters();
	void printResults();
	List returnResults();
	void  putParamsToZero();
	S4 returnClustering();
	double computeICL();

protected:
	mat _x;
	int _Nr;
	vector<int> _Jc;
	vector<int>_m;
	vector<urowvec> _dlist;
	vector<Distribution*> _distrib_objects;
	int _number_distrib;
	int _kr;
	vector<int> _zr;
	mat _probaV;
	mat _logprobaV;
	mat _V;
	rowvec _gamma;
	vector<rowvec> _allgamma;
	rowvec _resgamma;
	// params of co-clustering :
	string _init;
	int _nbSEM;
	int _nbSEMburn;
	int _nbindmini;

	mat _zrchain;

	// to sample distributions
	random_device _rd;

	double _icl;

	int _seed;

	// complex init
	std::vector<double> _percentRandomB;

	// Utils
	double logsum(rowvec logx);
	rowvec getMeans(mat VorW);
	mat kmeansi();
	double getDistance(vec &a, vec &b);
};

