#include "Distribution.h"


Distribution::Distribution(mat xsep, int kr, int kc, int nbSEM, int seed)
{
	this->_seed = seed;
	this->_nbSEM = nbSEM;
	this->_xsep = xsep;
	this->_Nr = xsep.n_rows;
	this->_Jc = xsep.n_cols;
	vector<vector<int>> miss_tmp;
	for (int i = 0; i < this->_Nr; i++)
	{
		for (int j = 0; j < this->_Jc; j++)
		{
			if (isnan(xsep(i, j))) {
				vector<int> coordinates;
				coordinates.push_back(i);
				coordinates.push_back(j);
				miss_tmp.push_back(coordinates);
			}
		}
	}
	this->_miss = miss_tmp;
	this->_kr = kr;
	this->_kc = kc;

}

Distribution::Distribution()
{
}


Distribution::~Distribution()
{
}

void Distribution::missingValuesInit() {
	return;
}


TabProbsResults Distribution::SEstep(mat V, mat W)
{
	TabProbsResults result(_Nr, _kr, _Jc, _kc);
	return result;
}

mat Distribution::SEstepRow(mat V, mat W){
	mat result(_Nr, _kr);
	return(result);
}

mat Distribution::SEstepCol(mat V, mat W){
	mat result(_Jc, _kc);
	return(result);
}

LogProbs Distribution::SEstep_predict(int i, int d, int k, int h, double x_id, double sumi, double sumd, vec x_id_vec)
{
	LogProbs result(0, 0);
	return result;
}

void Distribution::imputeMissingData(mat V, mat W) {
	return;
}


void Distribution::Mstep(uvec rowind, uvec colind, int k, int h, bool init)
{
	return;
}

void Distribution::MstepVW(mat V, mat W, bool init)
{
	return;
}

rowvec Distribution::getDatablockkh(uvec rowind, uvec colind) {
	// function that exists so that we can get a vector of values for block kh
	// furthermore, we want to get rid of missing values:
	mat datablock_kh = _xsep.submat(rowind, colind);
	
	rowvec result = conv_to<rowvec>::from(vectorise(datablock_kh));
	// TODO: find a way for missing values:
	uvec todelete = find(result == -1);
	for (int i = 0; i < todelete.n_elem; i++) {
		result.shed_col(todelete(i));
	}
	return result;
}

void Distribution::fillParameters(int iteration) {
	return;
}

void Distribution::getBurnedParameters(int burn) {
	return;
}

void Distribution::printResults() {
	return;
}

List Distribution::returnResults() {
	List x;
	return(x);
}

List Distribution::returnParamsChain() {
	List x;
	return(x);
}

void Distribution::putParamsToZero() {
	return;
}

double Distribution::computeICL(int i, int d, int k, int h) {
	return(0);
}
bool Distribution::verif(mat V, mat W, int nbindmini) {
	bool result = true;
	for (int k = 0; k < _kr; k++) {
		for (int h = 0; h < _kc; h++) {
			uvec indicesV = arma::find(V.col(k) == 1);
			uvec indicesW = arma::find(W.col(h) == 1);
			int size = indicesV.n_elem * indicesW.n_elem;
			if (size < nbindmini) {
				return false;
			}
		}
	}

	return result;
}

int Distribution::verification(const mat& V, const mat& W, int nbindmini) {
	int result = -1;
	for (int k = 0; k < _kr; k++) {
		for (int h = 0; h < _kc; h++) {
			uvec indicesV = arma::find(V.col(k) == 1);
			uvec indicesW = arma::find(W.col(h) == 1);
			int size = indicesV.n_elem * indicesW.n_elem;
			if (size < nbindmini) {
				
				if(indicesV.n_elem<indicesW.n_elem){
					return ((-k));
				}
				else{
					return (h+1);
				}
			}
		}
	}

	return result;
}

mat Distribution::returnXhat(){
	return(this->_xsep);
}


mat Distribution::colkmeans() {
	mat result(_Jc, _kc);
	result.zeros();

	mat colmeans;

	bool status = arma::kmeans(colmeans, _xsep, _kc, random_subset, 3, false);
	if (status == false)
	{
		return result;
	}

	for (int d = 0; d < _Jc; d++) {
		int num_clust = -1;
		double dst_old = -1;
		double dst = -1;

		for (int h = 0; h < _kc; h++) {
			vec a(_Nr);
			vec b(_Nr);
			for (int ireconstruct = 0; ireconstruct < _Nr; ireconstruct++) {
				a(ireconstruct) = colmeans.col(h)(ireconstruct);
				b(ireconstruct) = _xsep.col(d)(ireconstruct);
			}
			dst = this->getDistance(a, b);
			if (dst_old < 0 || dst < dst_old) {
				dst_old = dst;
				num_clust = h;
			}
		}
		result(d, num_clust) = 1;
	}
	return result;
}

double Distribution::getDistance(vec &a, vec &b) {
	vec temp = a - b;
	return arma::norm(temp);

}
