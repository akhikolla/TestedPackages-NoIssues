#pragma once
#include "Distribution.h"
#include "Bos.h"

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


class ClassificationMContext
{
public:
	ClassificationMContext(arma::mat x, arma::vec y, std::vector< arma::urowvec > dlist, int kr, 
		std::string init, int nbSEM, 
		int nbSEMburn, int nbindmini, 
		std::vector< int > m, int seed);
	ClassificationMContext();
	~ClassificationMContext();
	void missingValuesInit(); 
	void initialization();
	void Mstep();
	void MstepVW();
	void SEstep();
	void sampleVW();
	void sampleVWStock();
	void imputeMissingData();
	bool verif();
	void fillParameters(int iteration);
	void getBurnedParameters();
	void printResults();
	void returnResults();
	void  putParamsToZero();
	double computeICL();
	void predict();
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
	vector<int> _zr;
	mat _V;
	vector<mat> _W;
	rowvec _gamma;
	vector<rowvec> _allgamma;
	rowvec _resgamma;
	// params of co-clustering :
	string _init;
	int _nbSEM;
	int _nbSEMburn;
	int _nbindmini;
	// params for classificationMContext
	vec _y;

	// to sample distributions
	random_device _rd;

	double _icl;
	int _seed;

	// Utils
	rowvec getMeans(mat VorW);
	double logsum(rowvec logx);
	mat kmeansi();
	double getDistance(vec &a, vec &b);
};

