#include "ClassificationMContext.h"

ClassificationMContext::ClassificationMContext(arma::mat x, arma::vec y, std::vector< arma::urowvec > dlist,
	int kr, std::string init, int nbSEM, int nbSEMburn, 
	int nbindmini, std::vector< int > m, int seed)
{
	this->_seed = seed;
	// attributes that are directly instanciated
	this->_x = x;
	this->_Nr = _x.n_rows;
	this->_y = y;
	this->_dlist = dlist;
	this->_kr = kr;
	this->_m = m;
	this->_init = init;
	this->_nbSEM = nbSEM;
	this->_nbSEMburn = nbSEMburn;
	this->_nbindmini = nbindmini;

	//attributes to construct
	this->_number_distrib = m.size();

	// attributes regarding columns
	vector<int> tmp_Jc(_number_distrib);
	vector<Distribution*> tmp_distrib_objets;

	vector<mat> tmp_Wvec(_number_distrib);

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
		
		mat tmp_W(tmp_Jc[idistrib], tmp_Jc[idistrib], fill::zeros);
		tmp_W.eye();
		tmp_Wvec[idistrib] = tmp_W;
		
	}
	this->_Jc = tmp_Jc;
	this->_distrib_objects = tmp_distrib_objets;
	this->_W = tmp_Wvec;


	// attributes regarding lines
	this->_Nr = _x.n_rows;
	
	
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
	vector<rowvec> tmp_allgamma(_nbSEM);
	for (int isem = 0; isem < _nbSEM; isem++) {
		tmp_allgamma.push_back(this->_gamma);
	}
	this->_allgamma = tmp_allgamma;

}

ClassificationMContext::ClassificationMContext()
{
}


ClassificationMContext::~ClassificationMContext()
{
}

void ClassificationMContext::initialization() {
	//cout << "=============== initialization ===============" << endl;
	//arma_rng::set_seed_random();
	if (_init == "random") {
		
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++)
		{
			
			mat W;
			W.eye(_Jc.at(idistrib),_Jc.at(idistrib));
			this->_distrib_objects[idistrib]->MstepVW(_V, W, false);
			
		}

	}
	if (_init == "kmeans") {
		for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {

			mat W;
			W.eye(_Jc.at(idistrib),_Jc.at(idistrib));
			this->_distrib_objects[idistrib]->MstepVW(_V, W, false);
			// updating rho 
		}
	}
}

void ClassificationMContext::missingValuesInit() {
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->missingValuesInit();
	}
	return;
}

void ClassificationMContext::Mstep() {
	//cout << "=============== M step ===============" << endl;
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

void ClassificationMContext::MstepVW() {
	//cout << "=============== M step ===============" << endl;
	//cout << "mix gamma " << endl;
	//_gamma.print();
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		mat W;
		W.eye(_Jc.at(idistrib),_Jc.at(idistrib));
		this->_distrib_objects[idistrib]->MstepVW(_V, W, false);
		//cout << "mix rho " << idistrib << endl;
		//_rho.at(idistrib).print();
	}
}



void ClassificationMContext::SEstep()
{
	//cout << "=============== SE step ===============" << endl;

	// nothing here
}

void ClassificationMContext::sampleVW() {
	// Sampling W

	// nothing here
}

void ClassificationMContext::sampleVWStock() {
	// Sampling V and W

	//nothing here
}

void ClassificationMContext::imputeMissingData() {
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		mat W;
		W.eye(_Jc.at(idistrib),_Jc.at(idistrib));
		this->_distrib_objects[idistrib]->imputeMissingData(this->_V, W);
	}
}

bool ClassificationMContext::verif() {
	bool result = true;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		mat W;
		W.eye(_Jc.at(idistrib),_Jc.at(idistrib));
		result = this->_distrib_objects[idistrib]->verif(_V, W, _nbindmini);
		if (result == false) {
			return false;
		}
	}
	return result;
}

void ClassificationMContext::fillParameters(int iteration) {
	//cout << "=============== Filling parameters ===============" << endl;
	//this->_allgamma.at(iteration) = this->_gamma;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->fillParameters(iteration);
	}
}


void ClassificationMContext::getBurnedParameters() {
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

	// distributions parameters
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->getBurnedParameters(this->_nbSEMburn);
	}
}

void ClassificationMContext::printResults() {
	//cout << "=============== Printing parameters ===============" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->printResults();
	}
	//cout << "mix gamma" << endl;
	//_resgamma.print();
}

void ClassificationMContext::returnResults() {
	//cout << "=============== Returning parameters ===============" << endl;
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->returnResults();
	}
	//cout << "mix gamma" << endl;
	//_resgamma.print();
	//cout << "mix rho" << endl;
}

void ClassificationMContext::putParamsToZero(){
	for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
		this->_distrib_objects[idistrib]->putParamsToZero();
	}
}


S4 ClassificationMContext::returnClassification() {
	S4 x("ResultClassifOrdinal");

	x.slot("name") = "ClassifM";

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

 
	// parameters: 
	List resultAlpha(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultAlpha[idistrib] =	this->_distrib_objects[idistrib]->returnResults();
    	
	}
	x.slot("params") = resultAlpha;

	//xhat
	List resultXhat(_number_distrib);
    for (int idistrib = 0; idistrib < _number_distrib; idistrib++) {
    	
    	resultXhat[idistrib] =	this->_distrib_objects[idistrib]->returnXhat();
    	
	}
	x.slot("xhat") = resultXhat;

	// mixing proportions: 
	x.slot("pi")  = this->_resgamma;

    //icl
	x.slot("icl") = this->_icl;

	// ************** return to predict: *****************
	x.slot("kr") = _kr;
	x.slot("number_distrib") = _number_distrib;
	x.slot("J") = _Jc;
	x.slot("kc") = _Jc;


	List resultDlist(_dlist.size());
	for(int idistrib = 0; idistrib < _dlist.size(); idistrib++){
		resultDlist[idistrib] = _dlist[idistrib];
	}
	x.slot("dlist") = resultDlist;

	
	x.slot("m") = _m;
	x.slot("nbSEM") = _nbSEM;

    return(x);
	
}

double ClassificationMContext::computeICL() {
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
						result += _V(i, k)*_distrib_objects[idistrib]->computeICL(i, d, k, d)*_Jc[idistrib];
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
	return(result);
}



/*======================================UTILS======================================*/


rowvec ClassificationMContext::getMeans(mat VorW) {
	rowvec result;
	result.zeros(VorW.n_cols);
	for (int i = 0; i < VorW.n_cols; i++)
	{
		colvec column = VorW.col(i);
		result(i) = mean(column);
	}
	return result;
}

double ClassificationMContext::logsum(rowvec logx) {
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
