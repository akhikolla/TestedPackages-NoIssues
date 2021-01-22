#include "ClassificationContext.h"

ClassificationContext::ClassificationContext(arma::mat x, arma::vec y, 
											std::vector< arma::urowvec > dlist, int kr, 
											std::vector< int > kc, std::string init, int nbSEM, 
											int nbSEMburn, int nbindmini, std::vector< int > m, 
											std::vector<double> percentRandomB, int seed)
{

	this->_seed = seed;
	// attributes that are directly instanciated
	this->_x = x;
	this->_Nr = _x.n_rows;
	this->_y = y;
	this->_dlist = dlist;
	this->_kr = kr;
	this->_kc = kc;
	this->_m = m;
	this->_init = init;
	this->_nbSEM = nbSEM;
	this->_nbSEMburn = nbSEMburn;
	this->_nbindmini = nbindmini;

	//attributes to construct
	this->_number_distrib = _m.size();

	// for complex init
	this->_percentRandomB = percentRandomB;

	// attributes regarding columns
	vector<int> tmp_Jc(_number_distrib);
	vector<Distribution*> tmp_distrib_objets;
	vector<vector<int>> tmp_zcvec(_number_distrib);
	vector<rowvec> tmp_rhovec(_number_distrib);
	vector<rowvec> tmp_resrhovec(_number_distrib);
	vector<mat> tmp_probaWvec(_number_distrib);
	vector<mat> tmp_logprobaWvec(_number_distrib);
	vector<mat> tmp_Wvec(_number_distrib);
	vector<mat> tmp_zcchainvec(_number_distrib);


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
			tmp_distrib_objets.push_back(new Bos(xsep, _kr, _kc.at(idistrib), _m[im], this->_nbSEM, 
			iterordiEM, this->_seed));
			im++;
		}
		tmp_Jc[idistrib] = _dlist.at(idistrib).size();
		
		vector<int> tmp_zc(tmp_Jc.at(idistrib));
		std::fill(tmp_zc.begin(), tmp_zc.end(), 0);
		tmp_zcvec[idistrib] = tmp_zc;
		rowvec tmp_rho(_kc.at(idistrib), fill::zeros);
		tmp_rhovec.push_back(tmp_rho);
		rowvec tmp_resrho(_kc.at(idistrib), fill::zeros);
		tmp_resrhovec.push_back(tmp_rho);

		// probas and logprobas tables
		mat tmp_probaW(tmp_Jc[idistrib], _kc.at(idistrib), fill::zeros);
		tmp_probaWvec[idistrib] = tmp_probaW;
		mat tmp_logprobaW(tmp_Jc[idistrib], _kc.at(idistrib), fill::zeros);
		tmp_logprobaWvec[idistrib] = tmp_probaW;
		mat tmp_W(tmp_Jc[idistrib], _kc.at(idistrib), fill::zeros);
		tmp_Wvec[idistrib] = tmp_W;

		mat tmp_zcchain(_nbSEM, tmp_Jc[idistrib], fill::zeros);
		tmp_zcchainvec[idistrib] = tmp_zcchain;
		
	}
	this->_Jc = tmp_Jc;
	this->_probaW = tmp_probaWvec;
	this->_logprobaW = tmp_logprobaWvec;
	this->_zc = tmp_zcvec;
	this->_rho = tmp_rhovec;
	this->_resrho = tmp_resrhovec;
	this->_distrib_objects = tmp_distrib_objets;
	this->_W = tmp_Wvec;

	this->_zcchain = tmp_zcchainvec;

	
	
	vector<int> tmp_zr(this->_Nr);
	std::fill(tmp_zr.begin(), tmp_zr.end(), 0);
	this->_zr = tmp_zr;

	for (int i = 0; i < _Nr; i++) {
		this->_zr[i] = y(i);
	}


	mat tmp_V(_Nr, _kr, fill::zeros);
	this->_V = tmp_V;
	for (int i = 0; i < _Nr; i++) {
		int col = y(i) - 1;
		_V(i, col) = 1;
	}


	this->_gamma = this->getMeans(this->_V);

	this->_resgamma = this->_gamma;

	// filling the parameters on SEM iterations
	vector<vector<rowvec>> tmp_allrho;
	vector<rowvec> tmp_allgamma(_nbSEM);
	for (int isem = 0; isem < _nbSEM; isem++) {
		tmp_allgamma.push_back(this->_gamma);

		vector<rowvec> tmp_rho2(_number_distrib);
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
			rowvec tmp_rh(_kc[idistrib]);
			std::fill(tmp_rh.begin(), tmp_rh.end(), 0);
			tmp_rho2.push_back(tmp_rh);
		}
		tmp_allrho.push_back(tmp_rho2);
	}
	this->_allrho = tmp_allrho;
	this->_allgamma = tmp_allgamma;

}

ClassificationContext::ClassificationContext()
{
}


ClassificationContext::~ClassificationContext()
{
}

void ClassificationContext::initialization() {
	if (_init == "random" || _init == "randomBurnin") {
		
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++)
		{
			// RANDOM
			/* vector<double> vec(_kc.at(idistrib));
			double prob = (double)1 / _kc.at(idistrib);
			std::fill(vec.begin(), vec.end(), prob);
			discrete_distribution<> d(vec.begin(), vec.end()); */
			
			boost::mt19937 generator(this->_seed);
			vector<double> vec(_kc.at(idistrib));
			double prob = (double)1 / _kc.at(idistrib);
			std::fill(vec.begin(), vec.end(), prob);
			boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());

			for (int i = 0; i<_Jc.at(idistrib); ++i) {
				// random!
				int sample = distribution(generator);
				this->_W.at(idistrib)(i, sample) = 1;
			}

			// replaced by:
			this->_distrib_objects[idistrib]->MstepVW(_V, _W.at(idistrib), false);
			// updating rho 
			this->_rho.at(idistrib) = this->getMeans(this->_W.at(idistrib));
		}

	}
	if (_init == "kmeans") {
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
			this->_W.at(idistrib) = this->_distrib_objects[idistrib]->colkmeans();
			this->_distrib_objects[idistrib]->MstepVW(_V, _W.at(idistrib), false);
			// updating rho 
			this->_rho.at(idistrib) = this->getMeans(this->_W.at(idistrib));
		}
	}
}

void ClassificationContext::missingValuesInit() {
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->missingValuesInit();
	}
	return;
}

void ClassificationContext::Mstep() {
	//cout << "=============== M step ===============" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int k = 0; k < _kr; k++) {
			for (int h = 0; h < _kc[idistrib]; h++) {
				uvec rowind = find(this->_V.col(k) == 1);
				uvec colind = find(this->_W[idistrib].col(h) == 1);
				this->_distrib_objects[idistrib]->Mstep(rowind, colind, k, h, false);
			}
		}
		this->_rho.at(idistrib) = this->getMeans(this->_W.at(idistrib));
	}
}

void ClassificationContext::MstepVW() {
	//cout << "=============== M step ===============" << endl;
	//cout << "mix gamma " << endl;
	//_gamma.print();
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->MstepVW(_V, _W.at(idistrib), false);
		this->_rho.at(idistrib) = this->getMeans(this->_W.at(idistrib));
		//cout << "mix rho " << idistrib << endl;
		//_rho.at(idistrib).print();
	}
}

/*void ClassificationContext::SEstep()
{
	//cout << "=============== SE step ===============" << endl;
	// Computing the log-probabilites

	for (int idistrib = 0; idistrib < _number_distrib; idistrib++)
	{
		this->_logprobaW.at(idistrib).zeros();
		LogProbs result(0, 0);
		for (int d = 0; d < _Jc[idistrib]; d++)
		{
			this->_logprobaW.at(idistrib).row(d) = log(this->_rho.at(idistrib));
			for (int h = 0; h < _kc[idistrib]; h++)
			{

				for (int i = 0; i < _Nr; i++)
				{

					for (int k = 0; k < _kr; k++)
					{
						result = _distrib_objects[idistrib]->SEstep(i, d, k, h);
						_logprobaW.at(idistrib)(d, h) += _V(i, k) * result._col;
					}
				}
			}
		}

	}

	// Computing the probabilites
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int d = 0; d < _Jc.at(idistrib); d++) {
			for (int h = 0; h < _kc.at(idistrib); h++) {
				_probaW.at(idistrib)(d, h) = exp(_logprobaW.at(idistrib)(d, h) - logsum(_logprobaW.at(idistrib).row(d)));
			}
		}
	}
}*/

void ClassificationContext::SEstep()
{
	//cout << "=============== SE step ===============" << endl;
	// Computing the log-probabilites
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++)
	{
		this->_logprobaW.at(idistrib).zeros();
		for (int d = 0; d < _Jc.at(idistrib); d++) {
			this->_logprobaW.at(idistrib).row(d) = log(this->_rho.at(idistrib));
		}
		TabProbsResults result(_Nr, _kr, _Jc.at(idistrib), _kc.at(idistrib));
		result = _distrib_objects[idistrib]->SEstep(_V, _W.at(idistrib));
		this->_logprobaW.at(idistrib) += result._tabprobaW;
	}

	// Computing the probabilites
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int d = 0; d < _Jc.at(idistrib); d++) {
			for (int h = 0; h < _kc.at(idistrib); h++) {
				this->_probaW.at(idistrib)(d, h) = exp(this->_logprobaW.at(idistrib)(d, h) - logsum(this->_logprobaW.at(idistrib).row(d)));
			}
		}
	}

}

void ClassificationContext::sampleVW() {
	// Sampling W

	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_W.at(idistrib).zeros();
		for (int d = 0; d < _Jc.at(idistrib); d++) {
			//RANDOM
			/* rowvec vec = _probaW.at(idistrib).row(d);
			discrete_distribution<> dis(vec.begin(), vec.end());
			mt19937 gen(_rd()); */

			boost::mt19937 generator(this->_seed);
			rowvec vec = _probaW.at(idistrib).row(d);
			boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());
			int sample = distribution(generator);

			this->_W.at(idistrib)(d, sample) = 1;
		}
	}
	return;
}

void ClassificationContext::sampleVWStock() {
	// Sampling V and W

	vector<mat> countW(_number_distrib);
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int d = 0; d < _Jc.at(idistrib); d++) {
			mat countWid(_Jc.at(idistrib), _kc.at(idistrib), fill::zeros);
			countW[idistrib] = countWid;
		}
	}
	for (int iter = 0; iter < _nbSEM; iter++) {
		
		
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
			this->_W.at(idistrib).zeros();
			for (int d = 0; d < _Jc.at(idistrib); d++) {
				//RANDOM
				/* rowvec vec = _probaW.at(idistrib).row(d);
				discrete_distribution<> dis(vec.begin(), vec.end());
				mt19937 gen(_rd());
				int sample = dis(gen); */

				boost::mt19937 generator(this->_seed);
				rowvec vec = _probaW.at(idistrib).row(d);
				boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());
				int sample = distribution(generator);

				this->_W.at(idistrib)(d, sample) = 1;
				countW.at(idistrib)(d, sample) += 1;

			}
		}
	}

	//determinging final partitions
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_W.at(idistrib).zeros();
		for (int d = 0; d < _Jc.at(idistrib); d++) {
			int maxind = countW.at(idistrib).row(d).index_max();
			this->_W.at(idistrib)(d, maxind) = 1;
		}
	}

	return;
}

void ClassificationContext::imputeMissingData() {
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->imputeMissingData(this->_V, this->_W.at(idistrib));
	}
}

bool ClassificationContext::verif() {
	bool result = true;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		result = this->_distrib_objects[idistrib]->verif(_V, _W.at(idistrib), _nbindmini);
		if (result == false) {
			return false;
		}
	}
	return result;
}

vector<vector<int>> ClassificationContext::verification() {
	vector<vector<int>> result;
	for (int idistrib = 0; idistrib < this->_number_distrib; idistrib++) {
		int verifD = this->_distrib_objects[idistrib]->verification(_V, _W.at(idistrib), _nbindmini);

		if (!(verifD==-1)) {			
			vector<int> newline(2);
			newline.at(0) = idistrib;
			newline.at(1) = verifD; 
			result.push_back(newline);
		}
		
	}
	return result;
}

void ClassificationContext::noColDegenerancy(vector<vector<int>> distrib_col, int iter){
	double percent = _percentRandomB[0]/100;

	for(int nb_degen = 0; nb_degen<distrib_col.size(); nb_degen++){
		int idistrib = distrib_col.at(nb_degen)[0];
		int VorW = distrib_col.at(nb_degen)[1];

		if(!(VorW<0)){
			int nbToSample = ceil(percent*_Jc[idistrib]);
			// RANDOM
			/* std::random_device rdtest;     // only used once to initialise (seed) engine
			std::mt19937 rng(rdtest());    // random-number engine used (Mersenne-Twister in this case)
			std::uniform_int_distribution<int> uniW(0,(int)(_Jc[idistrib]-1)); // guaranteed unbiased
			std::uniform_int_distribution<int> unikc(0,(int)(_kc[idistrib]-1)); */
			boost::mt19937 generator(VorW+iter);
			boost::random::uniform_int_distribution<int> uniW(0,(int)(_Jc[idistrib]-1)); 
			boost::random::uniform_int_distribution<int> unikc(0,(int)(_kc[idistrib]-1));
			for(int i = 0; i<nbToSample; i++){
				int column = uniW(generator);
				//cout << "column : " << column << endl;
				rowvec newSample(_kc[idistrib]);
				newSample.zeros();
				(this->_W[idistrib]).row(column) = newSample;

				int cluster = unikc(generator);
				this->_W[idistrib](column, cluster) = 1;

			}
		}

	}
	
}

void ClassificationContext::fillParameters(int iteration) {
	//cout << "=============== Filling parameters ===============" << endl;
	//this->_allgamma.at(iteration) = this->_gamma;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_allrho.at(iteration).at(idistrib) = this->_rho.at(idistrib);
		this->_distrib_objects[idistrib]->fillParameters(iteration);
	}
}

void ClassificationContext::fillLabels(int iteration) {
	//cout << "=============== Filling parameters ===============" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for(int d = 0; d<_Jc.at(idistrib); d++){
			uvec tmp = find(_W.at(idistrib).row(d)==1);
			int label = tmp(0);
			_zcchain.at(idistrib)(iteration,d) = label;
		}
	}
}

void ClassificationContext::getBurnedParameters() {
	//cout << "=============== Getting burned parameters ===============" << endl;
	// gammas: 
	/*rowvec gamma_result = conv_to<rowvec>::from(zeros(_kr));
	for (int i = _nbSEMburn; i < _nbSEM; i++) {
		for (int k = 0; k < _kr; k++) {
			gamma_result(k) += this->_allgamma.at(i)(k);
		}
	}
	this->_resgamma = gamma_result / (_nbSEM - _nbSEMburn);
	*/
	// rhos:
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {

		rowvec resrho_distrib = conv_to<rowvec>::from(zeros(_kc[idistrib]));
		for (int i = _nbSEMburn; i < _nbSEM; i++) {
			for (int h = 0; h < _kc[idistrib]; h++) {
				resrho_distrib(h) += this->_allrho.at(i).at(idistrib)(h);
			}
		}
		this->_resrho.at(idistrib) = resrho_distrib / (_nbSEM - _nbSEMburn);
	}

	// distributions parameters
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->getBurnedParameters(this->_nbSEMburn);
	}
}

void ClassificationContext::printResults() {
	//cout << "=============== Printing parameters ===============" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->printResults();
	}
	//cout << "mix gamma" << endl;
	//_resgamma.print();
	//cout << "mix rho" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		_resrho.at(idistrib).print();
	}
}

void ClassificationContext::returnResults() {
	//cout << "=============== Returning parameters ===============" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->returnResults();
	}
	//cout << "mix gamma" << endl;
	//_resgamma.print();
	//cout << "mix rho" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		_resrho.at(idistrib).print();
	}
}

void ClassificationContext::putParamsToZero(){
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->putParamsToZero();
	}
}


S4 ClassificationContext::returnClassification() {
	S4 x("ResultClassifOrdinal");

	x.slot("name") = "Classif";

    // partitions:
    x.slot("V")  = _V;
    List resultW(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		resultW[idistrib] = this->_W[idistrib];
	}
    x.slot("W") = resultW;

  
    // labels:
    vec zr = zeros(_Nr);
    for(int i=0; i<_Nr; i++){
    	uvec k = find(_V.row(i)==1);
    	zr(i) = k(0)+1;
    }
    x.slot("zr")  = zr;

  

    List resultzc(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	vec zc = zeros(_Jc[idistrib]);
    	for(int d=0; d<_Jc[idistrib]; d++){
    		uvec h = find(_W[idistrib].row(d)==1);
    		zc(d) = h(0)+1;
    	}
		resultzc[idistrib] = zc;
	}
	x.slot("zc") = resultzc;

	// parameters: 
	List resultAlpha(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultAlpha[idistrib] =	this->_distrib_objects[idistrib]->returnResults();
    	
	}
	x.slot("params") = resultAlpha;

	// parameters chain: 
	/*List resultAlphaChain(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultAlphaChain[idistrib] =	this->_distrib_objects[idistrib]->returnParamsChain();
    	
	}
	x.slot("paramschain") = resultAlphaChain;*/

	//xhat
	List resultXhat(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultXhat[idistrib] =	this->_distrib_objects[idistrib]->returnXhat();
    	
	}
	x.slot("xhat") = resultXhat;

	// mixing proportions: 
	x.slot("pi")  = this->_resgamma;
	List resultRho(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultRho[idistrib] = this->_resrho[idistrib];
    	
	}
	x.slot("rho") = resultRho;
	
    //icl
	x.slot("icl") = this->_icl;

	// ************** return to predict: *****************
	x.slot("kr") = _kr;
	x.slot("kc") = _kc;
	x.slot("number_distrib") = _number_distrib;
	x.slot("J") = _Jc;

	List resultDlist(_dlist.size());
	for(int idistrib = 0; idistrib < _dlist.size(); idistrib++){
		resultDlist[idistrib] = _dlist[idistrib];
	}
	x.slot("dlist") = resultDlist;


	x.slot("m") = _m;
	x.slot("nbSEM") = _nbSEM;

    return(x);
	
}

double ClassificationContext::computeICL() {
	double result = 0;
	result += -(_kr - 1) / 2 * log(_Nr);
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		result += -(_kc[idistrib] - 1) / 2 * log(_Jc[idistrib]) - _kc[idistrib] * _kr / 2 * log(_Nr*_Jc[idistrib]);
	}
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int d = 0; d < _Jc[idistrib]; d++)
		{
			for (int h = 0; h < _kc[idistrib]; h++)
			{
				for (int i = 0; i < _Nr; i++)
				{
					for (int k = 0; k < _kr; k++)
					{
						result += _V(i, k)*_W[idistrib](d, h)*_distrib_objects[idistrib]->computeICL(i, d, k, h);
					}
				}
			}
		}
	}

	for (int i = 0; i < _Nr; i++) {
		for (int k = 0; k < _kr; k++) {
			result += _V(i, k)*log(_resgamma(k));
		}
	}
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		for (int d = 0; d < _Jc[idistrib]; d++) {
			for (int h = 0; h < _kc[idistrib]; h++) {
				result += _W[idistrib](d, h)*log(_resrho[idistrib](h));
			}
		}
	}
	return(result);
}



/*======================================UTILS======================================*/


rowvec ClassificationContext::getMeans(mat VorW) {
	rowvec result;
	result.zeros(VorW.n_cols);
	for (int i = 0; i < VorW.n_cols; i++)
	{
		colvec column = VorW.col(i);
		result(i) = mean(column);
	}
	return result;
}

double ClassificationContext::logsum(rowvec logx) {
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
