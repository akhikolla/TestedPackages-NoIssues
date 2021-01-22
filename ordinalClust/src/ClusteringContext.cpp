#include "ClusteringContext.h"

ClusteringContext::ClusteringContext(arma::mat x, std::vector< arma::urowvec > dlist,
									int kr, std::string init, int nbSEM, int nbSEMburn, 
									int nbindmini, std::vector< int > m, 
									std::vector<double> percentRandomB, int seed) {


	this->_seed = seed;
	// attributes that are directly instanciated
	this->_x = x;
	this->_Nr = _x.n_rows;
	this->_dlist = dlist;
	this->_kr = kr;
	this->_m = m;
	this->_init = init;
	this->_nbSEM = nbSEM;
	this->_nbSEMburn = nbSEMburn;
	this->_nbindmini = nbindmini;


	//attributes to construct
	this->_number_distrib = _m.size();

	// atributes for complex init
	this->_percentRandomB = percentRandomB;

	// attributes regarding columns
	vector<int> tmp_Jc(_number_distrib);
	vector<Distribution*> tmp_distrib_objets;


	this->_zrchain = zeros(_nbSEM,_Nr);

	int im = 0;


	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		// instanciate Distribution :
		string distrib_name = "Bos";
		// define xsep
		mat xsep;
		xsep = _x.cols(_dlist.at(idistrib));
		
		// here, should be done with a map, but I can't for now

		
		if (distrib_name == "Bos") {
			unsigned int iterordiEM = 10;
			tmp_distrib_objets.push_back(new Bos(xsep, _kr, _dlist.at(idistrib).size(), _m[im], this->_nbSEM, 
			iterordiEM, this->_seed));
			im++;
		}


		tmp_Jc[idistrib] = _dlist.at(idistrib).size();
		
		
				
	}
	this->_Jc = tmp_Jc;
	this->_distrib_objects = tmp_distrib_objets;

	// attributes regarding lines
	this->_Nr = _x.n_rows;
	
	mat tmp_probaV(_Nr, _kr, fill::zeros);
	this->_probaV = tmp_probaV;
	mat tmp_logprobaV(_Nr, _kr, fill::zeros);
	this->_logprobaV = tmp_probaV;
	vector<int> tmp_zr(this->_Nr);
	std::fill(tmp_zr.begin(), tmp_zr.end(), 0);
	this->_zr = tmp_zr;


	mat tmp_V(_Nr, _kr, fill::zeros);
	this->_V = tmp_V;
	rowvec tmp_gamma(_kr);
	std::fill(tmp_gamma.begin(), tmp_gamma.end(), 0);
	this->_gamma = tmp_gamma;
	rowvec tmp_resgamma(_kr);
	std::fill(tmp_resgamma.begin(), tmp_resgamma.end(), 0);
	this->_resgamma = tmp_resgamma;

	// filling the parameters on SEM iterations
	vector<rowvec> tmp_allgamma(_nbSEM);
	for (int isem = 0; isem < _nbSEM; isem++) {
		rowvec tmp_gam(_kr);
		std::fill(tmp_gam.begin(), tmp_gam.end(), 0);
		tmp_allgamma.push_back(tmp_gam);
	}
	this->_allgamma = tmp_allgamma;
}
ClusteringContext::ClusteringContext()
{
}


ClusteringContext::~ClusteringContext()
{
}

void ClusteringContext::initialization() {
	//cout << "=============== initialization ===============" << endl;
	//arma_rng::set_seed_random();
	if (_init == "random" || _init == "randomBurnin") {
		// partitions V
		//RANDOM
		/* vector<double> vec(_kr);
		double prob = (double)1 / _kr;
		std::fill(vec.begin(), vec.end(), prob);
		discrete_distribution<> d(vec.begin(), vec.end()); */

		boost::mt19937 generator(this->_seed);
		vector<double> vec(_kr);
		double prob = (double)1 / _kr;
		std::fill(vec.begin(), vec.end(), prob);
		boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());


		for (int i = 0; i<_Nr; ++i) {
			// random!
			/* mt19937 gen(_rd());
			int sample = d(gen);*/
			int sample = distribution(generator);

			this->_V(i, sample) = 1;
		}
		// updating gamma
		this->_gamma = this->getMeans(this->_V);
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++)
		{
			for (int k = 0; k < _kr; k++) {
				for (int d = 0; d < _Jc[idistrib]; d++) {
					uvec rowind = find(this->_V.col(k) == 1);
					uvec colind;
					colind << d;
					this->_distrib_objects[idistrib]->Mstep(rowind, colind, k, d, true);
				}
			}
		}

	}
	if (_init == "kmeans") {
		this->_V = this->kmeansi();
		this->_gamma = this->getMeans(this->_V);
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
			for (int k = 0; k < _kr; k++) {
				for (int d = 0; d < _Jc[idistrib]; d++) {
					uvec rowind = find(this->_V.col(k) == 1);
					uvec colind;
					colind << d;
					this->_distrib_objects[idistrib]->Mstep(rowind, colind, k, d, true);
				}
			}
		}
	}
}

void ClusteringContext::SEstep()
{
	//cout << "=============== SE step ===============" << endl;
	// Computing the log-probabilites
	this->_logprobaV.zeros();
	for (int i = 0; i < _Nr; i++) {
		this->_logprobaV.row(i) = log(this->_gamma);
	}

	for (int idistrib = 0; idistrib < _number_distrib; idistrib++)
	{

		TabProbsResults result(_Nr, _kr, _Jc.at(idistrib), _Jc.at(idistrib));
		mat Wtmp = zeros(_Jc[idistrib], _Jc[idistrib]);
		Wtmp.eye();
		result = _distrib_objects[idistrib]->SEstep(_V, Wtmp);
		this->_logprobaV += result._tabprobaV;
	}


	// Computing the probabilites
	for (int i = 0; i < _Nr; i++) {
		for (int k = 0; k < _kr; k++) {
			this->_probaV(i, k) = exp(this->_logprobaV(i, k) - logsum(_logprobaV.row(i)));
		}
	}
}

void ClusteringContext::SEstepRow()
{
	// Computing the log-probabilites
	this->_logprobaV.zeros();
	this->_logprobaV.each_row() += log(this->_gamma);  


	for (int idistrib = 0; idistrib < _number_distrib; idistrib++)
	{
		mat Wtmp = zeros(_Jc[idistrib], _Jc[idistrib]);
		Wtmp.eye();
		mat result(_Nr, _kr);
		result = _distrib_objects[idistrib]->SEstepRow(_V, Wtmp);
		this->_logprobaV += result;
	}

	// Computing the probabilites
	for (int i = 0; i < _Nr; i++) {
		for (int k = 0; k < _kr; k++) {
			this->_probaV(i, k) = exp(this->_logprobaV(i, k) - logsum(_logprobaV.row(i)));
		}
	}
}

void ClusteringContext::missingValuesInit() {
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->missingValuesInit();
	}
	return;
}

void ClusteringContext::Mstep() {
	this->_gamma = this->getMeans(this->_V);
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int k = 0; k < _kr; k++) {
			for (int d = 0; d < _Jc[idistrib]; d++) {
				uvec rowind = find(this->_V.col(k) == 1);
				uvec colind;
				colind << d;
				this->_distrib_objects[idistrib]->Mstep(rowind, colind, k, d, false);
			}
		}
	}
}

void ClusteringContext::MstepVW() {
	this->_gamma = this->getMeans(this->_V);
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		mat Wtmp = zeros(_Jc[idistrib],_Jc[idistrib]);
		Wtmp.eye();
		this->_distrib_objects[idistrib]->MstepVW(_V, Wtmp, false);
	}
}

void ClusteringContext::noRowDegenerancy(vector<vector<int>> distrib_col, int iter){
	double percent = _percentRandomB[0]/100;

	int count = 0;
	for(int nb_degen = 0; nb_degen<distrib_col.size(); nb_degen++){
		int VorW = distrib_col.at(nb_degen)[1];
		//cout << "VorW: " << VorW << endl;
		if(VorW<0){
			count++;
			int nbToSample = ceil(percent*_Nr);
			// RANDOM
			/* std::random_device rdtest;     
			std::mt19937 rng(rdtest());    
			std::uniform_int_distribution<int> uniW(0,(int)(_Nr-1)); 
			std::uniform_int_distribution<int> unikr(0,(int)(_kr-1)); */
			boost::mt19937 generator((-VorW)+iter);
			boost::random::uniform_int_distribution<int> uniW(0,(int)(_Nr-1));
			boost::random::uniform_int_distribution<int> unikr(0,(int)(_kr-1));
			for(int i = 0; i<nbToSample; i++){
				int line = uniW(generator);
				//cout << "line: " << line << endl;
				rowvec newSample(_kr);
				newSample.zeros();
				(this->_V).row(line) = newSample;
				int cluster = unikr(generator);
				this->_V(line, cluster) = 1;

			}
		}
		if(count>0) return;
	}
	
}

void ClusteringContext::sampleVW() {
	// Sampling V and W

	this->_V.zeros();
	//std::default_random_engine gen;
	for (int i = 0; i < _Nr; i++) {
		//RANDOM
		/* rowvec vec = _probaV.row(i);
		discrete_distribution<> dis(vec.begin(), vec.end());
		mt19937 gen(_rd());
		int sample = dis(gen); */
		rowvec vec = _probaV.row(i);
		boost::mt19937 generator(this->_seed);
		boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());
		int sample = distribution(generator);
		this->_V(i, sample) = 1;
	}
	return;
}

void ClusteringContext::sampleVWStock() {
	// Sampling V and W

	mat countV = zeros(_Nr, _kr);
	
	for (int iter = 0; iter < _nbSEM; iter++) {
		this->_V.zeros();
		//std::default_random_engine gen;
		for (int i = 0; i < _Nr; i++) {
			//RANDOM
			/* rowvec vec = _probaV.row(i);
			discrete_distribution<> dis(vec.begin(), vec.end());
			mt19937 gen(_rd());
			int sample = dis(gen); */

			rowvec vec = _probaV.row(i);
			boost::mt19937 generator(this->_seed);
			boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());
			int sample = distribution(generator);
			this->_V(i, sample) = 1;
			countV(i, sample) += 1;
		}
	}

	//determinging final partitions
	this->_V.zeros();
	for (int i = 0; i < _Nr; i++) {
		int maxind = countV.row(i).index_max();
		this->_V(i, maxind) = 1;
	}
	//this->_V.print();
	return;
}

void ClusteringContext::imputeMissingData() {
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		mat diag;
		diag.eye(_Jc[idistrib], _Jc[idistrib]);
		this->_distrib_objects[idistrib]->imputeMissingData(this->_V, diag);
	}
}

bool ClusteringContext::verif() {
	bool result = true;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		mat diag;
		diag.eye(_Jc[idistrib], _Jc[idistrib]);
		result = this->_distrib_objects[idistrib]->verif(_V, diag, _nbindmini);
		if (result == false) {
			return false;
		}
	}
	return result;
}

vector<vector<int>> ClusteringContext::verification() {
	vector<vector<int>> result;
	for (int idistrib = 0; idistrib < this->_number_distrib; idistrib++) {
		mat Weye(_Jc.at(idistrib),_Jc.at(idistrib));
		Weye.eye();
		int verifD = this->_distrib_objects[idistrib]->verification(_V, Weye, _nbindmini);

		if (!(verifD==-1)) {			
			vector<int> newline(2);
			newline.at(0) = idistrib;
			newline.at(1) = verifD; 
			result.push_back(newline);
		}
	}
	return result;
}

void ClusteringContext::fillParameters(int iteration) {
	//cout << "=============== Filling parameters ===============" << endl;
	this->_allgamma.at(iteration) = this->_gamma;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->fillParameters(iteration);
	}
}

void ClusteringContext::fillLabels(int iteration) {
	//cout << "=============== Filling parameters ===============" << endl;
	for(int i = 0; i<_Nr; i++){
		uvec tmp = find(_V.row(i)==1);
		int label = tmp(0);
		_zrchain(iteration, i) = label;
	}
}

void ClusteringContext::getBurnedParameters() {
	//cout << "=============== Getting burned parameters ===============" << endl;
	// gammas:
	rowvec gamma_result = conv_to<rowvec>::from(zeros(_kr));
	for (int i = _nbSEMburn; i < _nbSEM; i++) {
		for (int k = 0; k < _kr; k++) {
			gamma_result(k) += this->_allgamma.at(i)(k);
		}
	}
	this->_resgamma = gamma_result / (_nbSEM - _nbSEMburn);

	// distributions parameters
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->getBurnedParameters(this->_nbSEMburn);
	}
}

void ClusteringContext::printResults() {
	//cout << "=============== Printing parameters ===============" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->printResults();
	}
	//cout << "mix gamma" << endl;
}

List ClusteringContext::returnResults() {
	//cout << "=============== Returning parameters ===============" << endl;
	List resultList(_number_distrib + 1);
	List mixings = List::create(Rcpp::Named("resgamma") = _resgamma);
	List Vs =  List::create(Rcpp::Named("V") = _V);

	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		resultList[idistrib] = this->_distrib_objects[idistrib]->returnResults();
	}
	resultList[_number_distrib] = mixings;

	
	return(resultList);
	
}

void ClusteringContext::putParamsToZero() {
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->putParamsToZero();
	}
}

S4 ClusteringContext::returnClustering() {
	S4 x("ResultClustOrdinal");

    // partitions:
    x.slot("V")  = _V;


    // labels:
    vec zr = zeros(_Nr);
    for(int i=0; i<_Nr; i++){
    	uvec k = find(_V.row(i)==1);
    	zr(i) = k(0)+1;
    }
    x.slot("zr")  = zr;

    //labels chain
	x.slot("zrchain") = _zrchain;


	// parameters: 
	List resultAlpha(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultAlpha[idistrib] =	this->_distrib_objects[idistrib]->returnResults();
    	
	}
	x.slot("params") = resultAlpha;

	// parameters chain: 
	List resultAlphaChain(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultAlphaChain[idistrib] =	this->_distrib_objects[idistrib]->returnParamsChain();
    	
	}
	x.slot("paramschain") = resultAlphaChain;

	//xhat
	List resultXhat(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultXhat[idistrib] =	this->_distrib_objects[idistrib]->returnXhat();
    	
	}
	x.slot("xhat") = resultXhat;


	// mixing proportions: 
	x.slot("pi")  = this->_resgamma;

	//pichain
	List resultAllGamma(_nbSEM);
    for (int i = 0; i < _nbSEM; i++) {
    	
    	resultAllGamma[i] = this->_allgamma.at(i);
    	
	}
	x.slot("pichain") = resultAllGamma;
	
	//icl
	x.slot("icl") = this->_icl;

	x.slot("name") = "Clust";
	x.slot("m") = _m;

    return(x);
	
}

double ClusteringContext::computeICL() {
	double result = 0;
	result += -(_kr - 1) / 2 * log(_Nr);
	// don't know if it should be here:
	/*for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		result += -(_Jc[idistrib] - 1) / 2 * log(_Jc[idistrib]) - _Jc[idistrib] * _kr / 2 * log(_Nr*_Jc[idistrib]);
	}*/
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int d = 0; d < _Jc[idistrib]; d++)
		{
				for (int i = 0; i < _Nr; i++)
				{
					for (int k = 0; k < _kr; k++)
					{
						result += _V(i, k)*_distrib_objects[idistrib]->computeICL(i, d, k, d);//*_Jc[idistrib];
					}
				}			
		}
	}

	for (int i = 0; i < _Nr; i++) {
		for (int k = 0; k < _kr; k++) {
			result += _V(i, k)*log(_resgamma(k));
		}
	}

	// don't know if it should be here: 
	/*for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int d = 0; d < _Jc[idistrib]; d++) {
			result += log(1/ _Jc[idistrib]);
		}
	}*/
	this->_icl = result;
	return(result);
}

/*======================================UTILS======================================*/


double ClusteringContext::logsum(rowvec logx) {
	logx.replace(datum::nan, -100000);
	if (logx.size() == 1) {
		return logx(0);
	}
	double result = 0;
	logx = sort(logx, "descend");
	double tmp = 1;
	for (int i = 1; i < logx.n_elem; i++) {
		tmp += exp(logx(i) - logx(0));
	}
	result = logx(0) + log(tmp);
	return(result);
}

rowvec ClusteringContext::getMeans(mat VorW) {
	rowvec result;
	result.zeros(VorW.n_cols);
	for (int i = 0; i < VorW.n_cols; i++)
	{
		colvec column = VorW.col(i);
		result(i) = mean(column);
	}
	return result;
}

mat ClusteringContext::kmeansi() {

	mat result(_Nr, _kr);
	result.zeros();

	mat means;
	mat xtmp = _x;
	xtmp.replace(datum::nan, 0); 
	bool status = arma::kmeans(means, xtmp.t(), _kr, random_subset, 3, false);
	if (status == false)
	{
		//cout << "clustering failed" << endl;
	}
	for (int i = 0; i < _Nr; i++) {
		//cout << i << endl;
		int num_clust = -1;
		double dst_old = -1;
		double dst = -1;
		int leng = std::accumulate(_Jc.begin(), _Jc.end(), 0);
		for (int k = 0; k < _kr; k++) {
			vec a(leng);
			vec b(leng);
			for (int ireconstruct = 0; ireconstruct < means.col(k).n_elem; ireconstruct++) {
				a(ireconstruct) = means.col(k)(ireconstruct);
				b(ireconstruct) = _x.row(i)(ireconstruct);
			}
			dst = this->getDistance(a, b);
			if (dst_old < 0 || dst < dst_old) {
				dst_old = dst;
				num_clust = k;
			}
		}

		result(i, num_clust) = 1;
	}
	return result;
}

double ClusteringContext::getDistance(vec &a, vec &b) {
	vec temp = a - b;
	return arma::norm(temp);

}

