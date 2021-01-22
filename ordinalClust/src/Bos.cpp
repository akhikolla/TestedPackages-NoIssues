#include "Bos.h"


const double inf = std::numeric_limits<double>::infinity();



Bos::Bos(mat xsep, int kr, int kc, int m, int nbSEM, int seed, unsigned int iterordiEM)
	:Distribution(xsep, kr, kc, nbSEM, seed)
{
	this->_name = "Bos";
	this->_m = m;
	this->_iterordiEM = iterordiEM;
	
	// the tmp parameters, that change at each iteration
	this->_pis = ones(_kr, _kc)*0.5;
	this->_mus = arma::conv_to<umat>::from(ones(_kr, _kc));

	// the stocked parameters for each iteration
	this->_allpis = ones(_kr, _kc, _nbSEM)*0.5;
	this->_allmus = arma::conv_to<ucube>::from(ones(_kr, _kc, _nbSEM));

	// the resulting parameters
	this->_respis = ones(_kr, _kc)*0.5;
	this->_resmus = arma::conv_to<umat>::from(ones(_kr, _kc));
	
	this->_tab_pejs = this->gettabpej();

	this->_xsepCube = zeros(_Nr, _Jc, _m);
	for(int i = 0; i < _Nr; i++){
		for(int d = 0; d < _Jc; d++){
			if(isnan(_xsep(i,d)!=1))
				this->_xsepCube(i , d, _xsep(i,d)-1) = 1;
		}
	}

	this->_cubeProbs = zeros(_Nr,_Jc,_m);
}

Bos::Bos()
{
}


Bos::~Bos()
{
}

void Bos::missingValuesInit() {
	//_xsepCube = zeros(_Nr, _Jc, _m);
	for (int imiss = 0; imiss < _miss.size(); imiss++) {
		
		// RANDOM
		/* mt19937 gen(_rd());
		double eqprob = (double)1 / _m;
		vec vecprob(_m, fill::ones);
		vecprob = vecprob*eqprob;
		discrete_distribution<> d(vecprob.begin(), vecprob.end()); // maybe a problem?
		int sample = d(gen);*/

		boost::mt19937 generator(this->_seed);
		double eqprob = (double)1 / _m;
		vec vecprob(_m, fill::ones);
		vecprob = vecprob*eqprob;
		boost::random::discrete_distribution<int> distribution (vecprob.begin(),vecprob.end());
		int sample = distribution(generator);

		_xsep(_miss.at(imiss)[0], _miss.at(imiss)[1]) = sample + 1;
		_xsepCube.tube(_miss.at(imiss)[0], _miss.at(imiss)[1]) = zeros(_m);
		_xsepCube(_miss.at(imiss)[0], _miss.at(imiss)[1], sample) = 1;
	}
	return;
}


/*LogProbs Bos::SEstep(int i, int d, int k, int h)
{
	// cout << "*****bos SE sep******" << endl;
	cube cubeProbs = this->getCubeProbs();
	LogProbs result(0, 0);
	double accum = cubeProbs(k, h, _xsep(i, d) - 1);

	result._row = log(accum);
	result._col = log(accum);

	return result;
}*/

TabProbsResults Bos::SEstep(mat V, mat W)
{
	// cout << "*****bos SE sep******" << endl;
	cube cubeProbs = this->getCubeProbs();
	TabProbsResults result = TabProbsResults(_Nr, _kr, _Jc, _kc);
	for (size_t i = 0; i < _Nr; i++)
	{

		for (size_t k = 0; k < _kr; k++)
		{

			for (size_t d = 0; d < _Jc; d++)
			{

				for (size_t h = 0; h < _kc; h++)
				{
					

					double accum = cubeProbs(k, h, _xsep(i, d) - 1);

					result._tabprobaV(i, k) = result._tabprobaV(i, k) +
						W(d, h) * log(accum);
					result._tabprobaW(d, h) = result._tabprobaW(d, h) +
						V(i, k) * log(accum);

					
				}
			}
		}
	}
	return result;
}

mat Bos::SEstepRow(mat V, mat W)
{
	// cout << "*****bos SE sep******" << endl;

	mat result(_Nr, _kr);
	result.zeros();

	for (size_t d = 0; d < _Jc; d++)
	{
		for (size_t h = 0; h < _kc; h++)
		{
			if(W(d, h)==1){
				for (size_t i = 0; i < _Nr; i++)
				{

					for (size_t k = 0; k < _kr; k++)
					{
						double accum = _cubeProbs(k, h, _xsep(i, d) - 1);
						result(i, k) +=  log(accum);	
					}
				}
			}
			
		}
	}




	return result;
}

mat Bos::SEstepCol(mat V, mat W)
{
	// cout << "*****bos SE sep******" << endl;

	mat result(_Jc, _kc);
	result.zeros();
	for (size_t i = 0; i < _Nr; i++)
	{

		for (size_t k = 0; k < _kr; k++)
		{
			if(V(i, k)==1){
				for (size_t d = 0; d < _Jc; d++)
				{

					for (size_t h = 0; h < _kc; h++)
					{
						double accum = _cubeProbs(k, h, _xsep(i, d) - 1);
						result(d, h) +=  log(accum);	
					}
				}
			}
		}
	}

	

	return result;
}

LogProbs Bos::SEstep_predict(int i, int d, int k, int h, double x_id, double sumi, double sumd, vec x_id_vec)
{
	//cout << "*****bos SE sep******" << endl;
	cube cubeProbs = this->getCubeProbs();
	LogProbs result(0, 0);
	double accum = cubeProbs(k, h, x_id - 1);

	result._row = log(accum);
	result._col = log(accum);

	return result;
}

void Bos::imputeMissingData(mat V, mat W) {
	// cout << " bos imputing missing data" << endl;

	for (int imiss = 0; imiss < _miss.size(); imiss++) {

		vector<int> coordinates = _miss.at(imiss);
		int miss_row = coordinates.at(0);
		arma::uvec k = arma::find(V.row(miss_row) == 1);
		int miss_col = coordinates.at(1);
		arma::uvec h = arma::find(W.row(miss_col) == 1);

		vector<double> proba(_m);

		
		for (int im = 0; im < _m; im++) {
			proba[im] = _cubeProbs(k(0),h(0),im);
		}

		// RANDOM
		/* mt19937 gen(_rd());
		discrete_distribution<> d(proba.begin(), proba.end());
		double sample = d(gen);*/

		boost::mt19937 generator(this->_seed);
		boost::random::discrete_distribution<int> distribution (proba.begin(),proba.end());
		int sample = distribution(generator);
		
		this->_xsep(miss_row, miss_col) = sample+1;
		
		_xsepCube.tube(_miss.at(imiss)[0], _miss.at(imiss)[1]) = zeros(_m);
		_xsepCube(_miss.at(imiss)[0], _miss.at(imiss)[1], sample) = 1;
	}


	return;
}

void Bos::Mstep(arma::uvec rowind, arma::uvec colind, int k, int h, bool init) 
{
	//cout << "BOS m step" << endl;
	colvec tabp0;
	colvec tabmu0 = linspace<colvec>(1, _m, _m);
	if (init) {
		tabp0 = linspace<colvec>(0, 1, 7);
	}
	else {
		if(_pis(k,h)>1e-3){
			tabp0 << _pis(k,h);
		}
		else{
			tabp0 = linspace<colvec>(0, 1, 7);
		}
	}

	colvec datablock_kh = arma::conv_to<colvec>::from(this->getDatablockkh(rowind, colind));
	//datablock_kh.print();
	Mupi mupi_kh = ordiemCpp(datablock_kh, tabmu0, tabp0, 1, this->_iterordiEM);
	/*cout << "mu = " << mupi_kh._mu << endl;
	cout << "pi = " << mupi_kh._pi << endl;*/
	this->_mus(k, h) = mupi_kh._mu;
	this->_pis(k, h) = mupi_kh._pi;

	this->_cubeProbs = this->getCubeProbs();


	return;
}

void Bos::MstepVW(mat V, mat W, bool init)
{
	//cout << "****** bos M step ******" << endl;
	for (int k = 0; k < _kr; k++)
	{
		for (int h = 0; h < _kc; h++) {
			mat datablock_kh;
			uvec rowind = find(V.col(k) == 1);
			uvec colind = find(W.col(h) == 1);
			datablock_kh = _xsep.submat(rowind, colind);
			colvec datablock_kh_col = vectorise(datablock_kh, 0);
			colvec tabmu0 = linspace<colvec>(1, _m, _m);
			colvec tabp0;
			if (init) {
				tabp0 = linspace<colvec>(0, 1, 7);
			}
			else {
				if(_pis(k,h)>1e-3){
					tabp0 << _pis(k,h);
				}
				else{
					tabp0 = linspace<colvec>(0, 1, 7);
				}
				
			}
			Mupi mupi_kh = ordiemCpp(datablock_kh_col, tabmu0, tabp0, 1, this->_iterordiEM);
			this->_mus(k, h) = mupi_kh._mu;
			this->_pis(k, h) = mupi_kh._pi;
		}
	}

	this->_cubeProbs = this->getCubeProbs();


	return;
}


void Bos::fillParameters(int iteration) { // TODO: check here if iteration ok (or iteration +1)
	this->_allmus.slice(iteration) = this->_mus;
	this->_allpis.slice(iteration) = this->_pis;
	return;
}

void Bos::getBurnedParameters(int burn) {
 
	// there is a problem with "median" for cube = doesn't work, so I make it myself
	for (int k = 0; k < _kr; k++) {
		for (int h = 0; h < _kc; h++) {
			vec tmp_pis_test = this->_allpis.tube(k, h);
			vec tmp_pis = tmp_pis_test.subvec(burn,_nbSEM-1);
			this->_respis(k,h) = median(tmp_pis);
			this->_pis(k,h) = this->_respis(k,h);
			// for the mode, there is no function...
			uvec tmp_mus_test = this->_allmus.tube(k, h);
			uvec tmp_mus = tmp_mus_test.subvec(burn,_nbSEM-1);
			this->_resmus(k,h) = this->getModeFromVec(tmp_mus);
			this->_mus(k,h) = this->_resmus(k,h);
		}
	}

	this->_cubeProbs = this->getCubeProbs();

	return;
}

cube Bos::getCubeProbs() {
	cube result(_kr, _kc, _m);
	result.zeros();
	for (int k = 0; k < _kr; k++) {
		for (int h = 0; h < _kc; h++) {
			for (int im = 0; im < _m; im++) {
				vec tmp_tab = this->_tab_pejs.tube(im, _mus(k, h) - 1);
				rowvec tmp_vec(_m);
				for (int im2 = 0; im2 < _m; im2++) {
					tmp_vec(im2) = std::pow(this->_pis(k, h), im2) * tmp_tab(im2);
				}
				result(k, h, im) = accu(tmp_vec);
			}
		}
	}
	return result;
}

cube Bos::gettabpej() {
	cube result(_m, _m, _m);
	result.ones();

	if (_m == 2) {
		result = result / 2;
		result(0, 1, 1) = -0.5;
		result(1, 0, 1) = -0.5;
	}
	if (_m == 3) {

		result.slice(0)  = result.slice(0) / 3;

		result.slice(1)={
		{ (double)11 / 18, (double)-5 / 18, (double) -8 / 18},
		{ (double) -1 / 6, (double)10 / 18, (double) -1 / 6},
		{ (double) -8 / 18, (double) -5 / 18, (double)11 / 18} 
		}; 


		result.slice(2)={
		{ (double)1 / 18, (double)-1 / 18, (double)1 / 9},
		{ (double) - 1 / 6, (double)1 / 9, (double)-1 / 6},
		{(double) 1 / 9, (double)-1 / 18, (double)1 / 18} 
		}; 

		/*result.slice(1) << (double)11 / 18 << (double)-5 / 18 << (double) -8 / 18 << endr
		<< (double) -1 / 6 << (double)10 / 18 << (double) -1 / 6 << endr
		<< (double) -8 / 18 << (double) -5 / 18 << (double)11 / 18 << endr;


		result.slice(2) << (double)1 / 18 << (double)-1 / 18 << (double)1 / 9 << endr
		<< (double) - 1 / 6 << (double)1 / 9 << (double)-1 / 6 << endr
		<< (double) 1 / 9 << (double)-1 / 18 << (double)1 / 18 << endr;*/
		  
	}
	if (_m == 4) {
		result.slice(0) = result.slice(0)/4;

		result.slice(1) = {
			{(double)61 / 96, (double)-5 / 32, (double)-5 / 16, (double)-19 / 48},
			{(double)-1 / 48, (double)55 / 96, (double)-5 / 48, (double)-7 / 32},
			{(double)-7 / 32, (double)-5 / 48, (double)55 / 96, (double)-1 / 48},
			{(double)-19 / 48, (double)-5 / 16, (double)-5 / 32, (double)61 / 96}
		};

		result.slice(2) = {
			{(double)1 / 9 , (double)-13 / 144 , (double)1 / 18 , (double)1 / 6 },
			{(double)-5 / 24, (double)1 / 6, (double)-19 / 144, (double)-5 / 72},
			{(double)-5 / 72, (double)-19 / 144, (double)1 / 6, (double)-5 / 24},
			{(double)1 / 6, (double)1 / 18, (double)-13 / 144, (double)1 / 9}
		};

		result.slice(3) = {
			{(double)1 / 288, (double)-1 / 288, (double)1 / 144, (double)-1 / 48},
			{(double)-1 / 48, (double)1 / 96, (double)-1 / 72, (double)11 / 288},
			{(double)11 / 288, (double)-1 / 72, (double)1 / 96, (double)-1 / 48},
			{(double)-1 / 48, (double)1 / 144, (double)-1 / 288, (double)1 / 288}
		};


		/*result.slice(1) << (double)61 / 96 << (double)-5 / 32 << (double)-5 / 16 << (double)-19 / 48 << endr
		<< (double)-1 / 48 << (double)55 / 96 << (double)-5 / 48 << (double)-7 / 32 << endr
		<< (double)-7 / 32 << (double)-5 / 48 << (double)55 / 96 << (double)-1 / 48 << endr
		<< (double)-19 / 48 << (double)-5 / 16 << (double)-5 / 32 << (double)61 / 96 << endr;

		result.slice(2) << (double)1 / 9 << (double)-13 / 144 << (double)1 / 18 << (double)1 / 6 << endr
		<< (double)-5 / 24 << (double)1 / 6 << (double)-19 / 144 << (double)-5 / 72 << endr
		<< (double)-5 / 72 << (double)-19 / 144 << (double)1 / 6 << (double)-5 / 24 << endr
		<< (double)1 / 6 << (double)1 / 18 << (double)-13 / 144 << (double)1 / 9 << endr;

		result.slice(3) << (double)1 / 288 << (double)-1 / 288 << (double)1 / 144 << (double)-1 / 48 << endr
		<< (double)-1 / 48 << (double)1 / 96 << (double)-1 / 72 << (double)11 / 288 << endr
		<< (double)11 / 288 << (double)-1 / 72 << (double)1 / 96 << (double)-1 / 48 << endr
		<< (double)-1 / 48 << (double)1 / 144 << (double)-1 / 288 << (double)1 / 288 << endr;*/

	}
	if (_m == 5) {
		result.slice(0) = result.slice(0) / 5;

		result.slice(1)={
			{(double)379 / 600 , (double)-49 / 600 , (double)-17 / 75 , (double)-23 / 75 , (double)-107 / 300},
			{(double)11 / 200, (double)17 / 30, (double)-11 / 200, (double)-33 / 200, (double)-137 / 600},
			{(double)-61 / 600, (double)-1 / 75, (double)169 / 300, (double)-1 / 75, (double)-61 / 600},
			{(double)-137 / 600, (double)-33 / 200, (double)-11 / 200, (double)17 / 30, (double)11 / 200},
			{(double)-107 / 300, (double)-23 / 75, (double)-17 / 75, (double)-49 / 600, (double)379 / 600}
		};

		result.slice(2)={
			{(double)127 / 800, (double)-263 / 2400, (double)17 / 1200, (double)47 / 400, (double)59 / 300},
			{(double)-761 / 3600, (double)379 / 1800, (double)-859 / 7200, (double)-457 / 7200, (double)-19 / 3600},
			{(double)-111 / 800, (double)-31 / 200, (double)757 / 3600, (double)-31 / 200, (double)-111 / 800},
			{(double)-19 / 3600, (double)-457 / 7200, (double)-859 / 7200, (double)379 / 1800, (double)-761 / 3600},
			{(double)59 / 300, (double)47 / 400, (double)17 / 1200, (double)-263 / 2400, (double)127 / 800}
		};

		result.slice(3)={
			{(double)17 / 1800, (double)-31 / 3600, (double)11 / 900, (double)-1 / 100, (double)-13 / 300},
			{(double)-19 / 450, (double)1 / 45, (double)-1 / 40, (double)2 / 75, (double)73 / 1800},
			{(double)8 / 225, (double)-109 / 3600, (double)23 / 900, (double)-109 / 3600, (double)8 / 225},
			{(double)73 / 1800, (double)2 / 75, (double)-1 / 40, (double)1 / 45, (double)-19 / 450},
			{(double)-13 / 300, (double)-1 / 100, (double)11 / 900, (double)-31 / 3600, (double)17 / 1800}
		};

		result.slice(4)={
			{(double)1 / 7200, (double)-1 / 7200, (double)1 / 3600, (double)-1 / 1200, (double)1 / 300},
			{(double)-1 / 720, (double)1 / 1800, (double)-1 / 1440, (double)13 / 7200, (double)-1 / 144},
			{(double)7 / 1440, (double)-1 / 720, (double)1 / 1200, (double)-1 / 720, (double)7 / 1440},
			{(double)-1 / 144, (double)13 / 7200, (double)-1 / 1440, (double)1 / 1800, (double)-1 / 720},
			{(double)1 / 300, (double)-1 / 1200, (double)1 / 3600, (double)-1 / 7200, (double)1 / 7200}
		};

		/*result.slice(1) << (double)379 / 600 << (double)-49 / 600 << (double)-17 / 75 << (double)-23 / 75 << (double)-107 / 300 << endr
			<< (double)11 / 200 << (double)17 / 30 << (double)-11 / 200 << (double)-33 / 200 << (double)-137 / 600 << endr
			<< (double)-61 / 600 << (double)-1 / 75 << (double)169 / 300 << (double)-1 / 75 << (double)-61 / 600 << endr
			<< (double)-137 / 600 << (double)-33 / 200 << (double)-11 / 200 << (double)17 / 30 << (double)11 / 200 << endr
			<< (double)-107 / 300 << (double)-23 / 75 << (double)-17 / 75 << (double)-49 / 600 << (double)379 / 600 << endr;

		result.slice(2) << (double)127 / 800 << (double)-263 / 2400 << (double)17 / 1200 << (double)47 / 400 << (double)59 / 300 << endr
			<< (double)-761 / 3600 << (double)379 / 1800 << (double)-859 / 7200 << (double)-457 / 7200 << (double)-19 / 3600 << endr
			<< (double)-111 / 800 << (double)-31 / 200 << (double)757 / 3600 << (double)-31 / 200 << (double)-111 / 800 << endr
			<< (double)-19 / 3600 << (double)-457 / 7200 << (double)-859 / 7200 << (double)379 / 1800 << (double)-761 / 3600 << endr
			<< (double)59 / 300 << (double)47 / 400 << (double)17 / 1200 << (double)-263 / 2400 << (double)127 / 800 << endr;

		result.slice(3) << (double)17 / 1800 << (double)-31 / 3600 << (double)11 / 900 << (double)-1 / 100 << (double)-13 / 300 << endr
			<< (double)-19 / 450 << (double)1 / 45 << (double)-1 / 40 << (double)2 / 75 << (double)73 / 1800 << endr
			<< (double)8 / 225 << (double)-109 / 3600 << (double)23 / 900 << (double)-109 / 3600 << (double)8 / 225 << endr
			<< (double)73 / 1800 << (double)2 / 75 << (double)-1 / 40 << (double)1 / 45 << (double)-19 / 450 << endr
			<< (double)-13 / 300 << (double)-1 / 100 << (double)11 / 900 << (double)-31 / 3600 << (double)17 / 1800 << endr;

		result.slice(4) << (double)1 / 7200 << (double)-1 / 7200 << (double)1 / 3600 << (double)-1 / 1200 << (double)1 / 300 << endr
		<< (double)-1 / 720 << (double)1 / 1800 << (double)-1 / 1440 << (double)13 / 7200 << (double)-1 / 144 << endr
		<< (double)7 / 1440 << (double)-1 / 720 << (double)1 / 1200 << (double)-1 / 720 << (double)7 / 1440 << endr
		<< (double)-1 / 144 << (double)13 / 7200 << (double)-1 / 1440 << (double)1 / 1800 << (double)-1 / 720 << endr
		<< (double)1 / 300 << (double)-1 / 1200 << (double)1 / 3600 << (double)-1 / 7200 << (double)1 / 7200 << endr;*/

	}
	if (_m == 6) {
		result.slice(0) = result.slice(0) / 6;

		result.slice(1) = {
			{(double)667 / 1080, (double)-7 / 216, (double)-361 / 2160, (double)-35 / 144, (double)-7 / 24, (double)-13 / 40},
			{(double)71 / 720, (double)397 / 720, (double)-7 / 360, (double)-133 / 1080, (double)-133 / 720, (double)-9 / 40},
			{(double)-1 / 30, (double)7 / 180, (double)1183 / 2160, (double)11 / 2160, (double)-11 / 135, (double)-287 / 2160},
			{(double)-287 / 2160, (double)-11 / 135, (double)11 / 2160, (double)1183 / 2160, (double)7 / 180, (double)-1 / 30},
			{(double)-9 / 40, (double)-133 / 720, (double)-133 / 1080, (double)-7 / 360, (double)397 / 720, (double)71 / 720},
			{(double)-13 / 40, (double)-7 / 24, (double)-35 / 144, (double)-361 / 2160, (double)-7 / 216, (double)667 / 1080}
		}; 

		result.slice(2) = {
			{(double)51413 / 259200, (double)-30983 / 259200, (double)-2011 / 129600, (double)1129 / 14400, (double)1663 / 10800, (double)461 / 2160},
			{(double)-52139 / 259200, (double)63773 / 259200, (double)-7231 / 64800, (double)-3389 / 51840, (double)-113 / 10800, (double)1607 / 43200},
			{(double)-42767 / 259200, (double)-4049 / 25920, (double)31543 / 129600, (double)-33517 / 259200, (double)-295 / 2592, (double)-21469 / 259200},
			{(double)-21469 / 259200, (double)-295 / 2592, (double)-33517 / 259200, (double)31543 / 129600, (double)-4049 / 25920, (double)-42767 / 259200},
			{(double)1607 / 43200, (double)-113 / 10800, (double)-3389 / 51840, (double)-7231 / 64800, (double)63773 / 259200, (double)-52139 / 259200},
			{(double)461 / 2160, (double)1663 / 10800, (double)1129 / 14400, (double)-2011 / 129600, (double)-30983 / 259200, (double)51413 / 259200}
		};

		result.slice(3) = {
			{(double)487 / 28800, (double)-137 / 9600, (double)73 / 4800, (double)-19 / 43200, (double)-331 / 10800, (double)-137 / 2160},
			{(double)-15617 / 259200, (double)8899 / 259200, (double)-4379 / 129600, (double)949 / 51840, (double)269 / 8100, (double)4103 / 129600},
			{(double)5771 / 259200, (double)-1481 / 32400, (double)23 / 576, (double)-2033 / 51840, (double)187 / 8100, (double)13697 / 259200},
			{(double)13697 / 259200, (double)187 / 8100, (double)-2033 / 51840, (double)23 / 576, (double)-1481 / 32400, (double)5771 / 259200},
			{(double)4103 / 129600, (double)269 / 8100, (double)949 / 51840, (double)-4379 / 129600, (double)8899 / 259200, (double)-15617 / 259200},
			{(double)-137 / 2160, (double)-331 / 10800, (double)-19 / 43200, (double)73 / 4800, (double)-137 / 9600, (double)487 / 28800}
		};

		result.slice(4) = {
			{(double)41 / 86400, (double)-13 / 28800, (double)11 / 14400, (double)-67 / 43200, (double)17 / 10800, (double)19 / 2160},
			{(double)-989 / 259200, (double)403 / 259200, (double)-59 / 32400, (double)181 / 51840, (double)-29 / 6480, (double)-1501 / 129600},
			{(double)2351 / 259200, (double)-461 / 129600, (double)11 / 4800, (double)-823 / 259200, (double)347 / 64800, (double)-763 / 259200},
			{(double)-763 / 259200, (double)347 / 64800, (double)-823 / 259200, (double)11 / 4800, (double)-461 / 129600, (double)2351 / 259200},
			{(double)-1501 / 129600, (double)-29 / 6480, (double)181 / 51840, (double)-59 / 32400, (double)403 / 259200, (double)-989 / 259200},
			{(double)19 / 2160, (double)17 / 10800, (double)-67 / 43200, (double)11 / 14400, (double)-13 / 28800, (double)41 / 86400}
		};

		result.slice(5) = {
			{(double)1 / 259200, (double)-1 / 259200, (double)1 / 129600, (double)-1 / 43200, (double)1 / 10800, (double)-1 / 2160},
			{(double)-1 / 17280, (double)1 / 51840, (double)-1 / 43200, (double)1 / 17280, (double)-7 / 32400, (double)137 / 129600},
			{(double)17 / 51840, (double)-1 / 12960, (double)1 / 25920, (double)-1 / 17280, (double)1 / 5400, (double)-1 / 1152},
			{(double)-1 / 1152, (double)1 / 5400, (double)-1 / 17280, (double)1 / 25920, (double)-1 / 12960, (double)17 / 51840},
			{(double)137 / 129600, (double)-7 / 32400, (double)1 / 17280, (double)-1 / 43200, (double)1 / 51840, (double)-1 / 17280},
			{(double)-1 / 2160, (double)1 / 10800, (double)-1 / 43200, (double)1 / 129600, (double)-1 / 259200, (double)1 / 259200}
		};
		 
		 

		/*result.slice(1) << (double)667 / 1080 << (double)-7 / 216 << (double)-361 / 2160 << (double)-35 / 144 << (double)-7 / 24 << (double)-13 / 40 << endr
			<< (double)71 / 720 << (double)397 / 720 << (double)-7 / 360 << (double)-133 / 1080 << (double)-133 / 720 << (double)-9 / 40 << endr
			<< (double)-1 / 30 << (double)7 / 180 << (double)1183 / 2160 << (double)11 / 2160 << (double)-11 / 135 << (double)-287 / 2160 << endr
			<< (double)-287 / 2160 << (double)-11 / 135 << (double)11 / 2160 << (double)1183 / 2160 << (double)7 / 180 << (double)-1 / 30 << endr
			<< (double)-9 / 40 << (double)-133 / 720 << (double)-133 / 1080 << (double)-7 / 360 << (double)397 / 720 << (double)71 / 720 << endr
			<< (double)-13 / 40 << (double)-7 / 24 << (double)-35 / 144 << (double)-361 / 2160 << (double)-7 / 216 << (double)667 / 1080 << endr;

		result.slice(2) << (double)51413 / 259200 << (double)-30983 / 259200 << (double)-2011 / 129600 << (double)1129 / 14400 << (double)1663 / 10800 << (double)461 / 2160 << endr
			<< (double)-52139 / 259200 << (double)63773 / 259200 << (double)-7231 / 64800 << (double)-3389 / 51840 << (double)-113 / 10800 << (double)1607 / 43200 << endr
			<< (double)-42767 / 259200 << (double)-4049 / 25920 << (double)31543 / 129600 << (double)-33517 / 259200 << (double)-295 / 2592 << (double)-21469 / 259200 << endr
			<< (double)-21469 / 259200 << (double)-295 / 2592 << (double)-33517 / 259200 << (double)31543 / 129600 << (double)-4049 / 25920 << (double)-42767 / 259200 << endr
			<< (double)1607 / 43200 << (double)-113 / 10800 << (double)-3389 / 51840 << (double)-7231 / 64800 << (double)63773 / 259200 << (double)-52139 / 259200 << endr
			<< (double)461 / 2160 << (double)1663 / 10800 << (double)1129 / 14400 << (double)-2011 / 129600 << (double)-30983 / 259200 << (double)51413 / 259200 << endr;

		result.slice(3) << (double)487 / 28800 << (double)-137 / 9600 << (double)73 / 4800 << (double)-19 / 43200 << (double)-331 / 10800 << (double)-137 / 2160 << endr
			<< (double)-15617 / 259200 << (double)8899 / 259200 << (double)-4379 / 129600 << (double)949 / 51840 << (double)269 / 8100 << (double)4103 / 129600 << endr
			<< (double)5771 / 259200 << (double)-1481 / 32400 << (double)23 / 576 << (double)-2033 / 51840 << (double)187 / 8100 << (double)13697 / 259200 << endr
			<< (double)13697 / 259200 << (double)187 / 8100 << (double)-2033 / 51840 << (double)23 / 576 << (double)-1481 / 32400 << (double)5771 / 259200 << endr
			<< (double)4103 / 129600 << (double)269 / 8100 << (double)949 / 51840 << (double)-4379 / 129600 << (double)8899 / 259200 << (double)-15617 / 259200 << endr
			<< (double)-137 / 2160 << (double)-331 / 10800 << (double)-19 / 43200 << (double)73 / 4800 << (double)-137 / 9600 << (double)487 / 28800 << endr;

		result.slice(4) << (double)41 / 86400 << (double)-13 / 28800 << (double)11 / 14400 << (double)-67 / 43200 << (double)17 / 10800 << (double)19 / 2160 << endr
			<< (double)-989 / 259200 << (double)403 / 259200 << (double)-59 / 32400 << (double)181 / 51840 << (double)-29 / 6480 << (double)-1501 / 129600 << endr
			<< (double)2351 / 259200 << (double)-461 / 129600 << (double)11 / 4800 << (double)-823 / 259200 << (double)347 / 64800 << (double)-763 / 259200 << endr
			<< (double)-763 / 259200 << (double)347 / 64800 << (double)-823 / 259200 << (double)11 / 4800 << (double)-461 / 129600 << (double)2351 / 259200 << endr
			<< (double)-1501 / 129600 << (double)-29 / 6480 << (double)181 / 51840 << (double)-59 / 32400 << (double)403 / 259200 << (double)-989 / 259200 << endr
			<< (double)19 / 2160 << (double)17 / 10800 << (double)-67 / 43200 << (double)11 / 14400 << (double)-13 / 28800 << (double)41 / 86400 << endr;

		result.slice(5) << (double)1 / 259200 << (double)-1 / 259200 << (double)1 / 129600 << (double)-1 / 43200 << (double)1 / 10800 << (double)-1 / 2160 << endr
			<< (double)-1 / 17280 << (double)1 / 51840 << (double)-1 / 43200 << (double)1 / 17280 << (double)-7 / 32400 << (double)137 / 129600 << endr
			<< (double)17 / 51840 << (double)-1 / 12960 << (double)1 / 25920 << (double)-1 / 17280 << (double)1 / 5400 << (double)-1 / 1152 << endr
			<< (double)-1 / 1152 << (double)1 / 5400 << (double)-1 / 17280 << (double)1 / 25920 << (double)-1 / 12960 << (double)17 / 51840 << endr
			<< (double)137 / 129600 << (double)-7 / 32400 << (double)1 / 17280 << (double)-1 / 43200 << (double)1 / 51840 << (double)-1 / 17280 << endr
			<< (double)-1 / 2160 << (double)1 / 10800 << (double)-1 / 43200 << (double)1 / 129600 << (double)-1 / 259200 << (double)1 / 259200 << endr;*/
	}
	if (_m == 7) {
		result.slice(0) = result.slice(0) / 7;

		result.slice(1)={
			{(double)529 / 882, (double)4 / 2205, (double)-437 / 3528, (double)-1151 / 5880, (double)-713 / 2940, (double)-809 / 2940, (double)-293 / 980 },
			{(double)221 / 1764, (double)2351 / 4410, (double)113 / 17640, (double)-107 / 1176, (double)-3 / 20, (double)-557 / 2940, (double)-213 / 980 },
			{(double)29 / 2940, (double)209 / 2940, (double)1166 / 2205, (double)193 / 8820, (double)-215 / 3528, (double)-281 / 2520, (double)-641 / 4410 },
			{(double)-323 / 4410, (double)-527 / 17640, (double)743 / 17640, (double)1168 / 2205, (double)743 / 17640, (double)-527 / 17640, (double)-323 / 4410 },
			{(double)-641 / 4410, (double)-281 / 2520, (double)-215 / 3528, (double)193 / 8820, (double)1166 / 2205, (double)209 / 2940, (double)29 / 2940 },
			{(double)-213 / 980, (double)-557 / 2940, (double)-3 / 20, (double)-107 / 1176, (double)113 / 17640, (double)2351 / 4410, (double)221 / 1764 },
			{(double)-293 / 980, (double)-809 / 2940, (double)-713 / 2940, (double)-1151 / 5880, (double)-437 / 3528, (double)4 / 2205, (double)529 / 882}
		};
			
		result.slice(2)={
			{(double)163103 / 705600, (double)-262159 / 2116800, (double)-38933 / 1058400, (double)16901 / 352800, (double)5237 / 44100, (double)311 / 1764, (double)3929 / 17640 },
			{(double)-65759 / 352800, (double)72697 / 264600, (double)-222569 / 2116800, (double)-48341 / 705600, (double)-6523 / 352800, (double)3223 / 117600, (double)2593 / 39200 },
			{(double)-9139 / 52920, (double)-314971 / 2116800, (double)63683 / 235200, (double)-239119 / 2116800, (double)-4387 / 43200, (double)-151427 / 2116800, (double)-87853 / 2116800 },
			{(double)-63211 / 529200, (double)-284233 / 2116800, (double)-269987 / 2116800, (double)17671 / 66150, (double)-269987 / 2116800, (double)-284233 / 2116800, (double)-63211 / 529200 },
			{(double)-87853 / 2116800, (double)-151427 / 2116800, (double)-4387 / 43200, (double)-239119 / 2116800, (double)63683 / 235200, (double)-314971 / 2116800, (double)-9139 / 52920 },
			{(double)2593 / 39200, (double)3223 / 117600, (double)-6523 / 352800, (double)-48341 / 705600, (double)-222569 / 2116800, (double)72697 / 264600, (double)-65759 / 352800 },
			{(double)3929 / 17640, (double)311 / 1764, (double)5237 / 44100, (double)16901 / 352800, (double)-38933 / 1058400, (double)-262159 / 2116800, (double)163103 / 705600}
		};
				
		result.slice(3)={
			{(double) 319847 / 12700800, (double)-50509 / 2540160, (double)104207 / 6350400, (double)991 / 141120, (double)-289 / 15120, (double)-5293 / 105840, (double)-71 / 882 },
			{(double)-37949 / 508032, (double)58843 / 1270080, (double)-19211 / 470400, (double)49337 / 4233600, (double)59933 / 2116800, (double)58621 / 2116800, (double)10361 / 529200 },
			{(double)92443 / 12700800, (double)-745111 / 12700800, (double)677017 / 12700800, (double)-590339 / 12700800, (double)33293 / 2540160, (double)534701 / 12700800, (double)12829 / 235200 },
			{(double)123481 / 2540160, (double)17551 / 1411200, (double)-216679 / 4233600, (double)176569 / 3175200, (double)-216679 / 4233600, (double)17551 / 1411200, (double)123481 / 2540160 },
			{(double)12829 / 235200, (double)534701 / 12700800, (double)33293 / 2540160, (double)-590339 / 12700800, (double)677017 / 12700800, (double)-745111 / 12700800, (double)92443 / 12700800 },
			{(double)10361 / 529200, (double)58621 / 2116800, (double)59933 / 2116800, (double)49337 / 4233600, (double)-19211 / 470400, (double)58843 / 1270080, (double)-37949 / 508032 },
			{(double)-71 / 882, (double)-5293 / 105840, (double)-289 / 15120, (double)991 / 141120, (double)104207 / 6350400, (double)-50509 / 2540160, (double)319847 / 12700800}
		};
				
		result.slice(4)={
			{(double)1433 / 1411200, (double)-29 / 31360, (double)319 / 235200, (double)-17 / 8640, (double)-17 / 105840, (double)661 / 105840, (double)3 / 196 },
			{(double)-86999 / 12700800, (double)2627 / 907200, (double)-40651 / 12700800, (double)61169 / 12700800, (double)-14227 / 6350400, (double)-58231 / 6350400, (double)-1783 / 132300 },
			{(double)21431 / 1814400, (double)-26429 / 4233600, (double)10547 / 2540160, (double)-65033 / 12700800, (double)15503 / 2540160, (double)-4559 / 4233600, (double)-13289 / 1270080 },
			{(double)33743 / 12700800, (double)105073 / 12700800, (double)-76331 / 12700800, (double)229 / 50400, (double)-76331 / 12700800, (double)105073 / 12700800, (double)33743 / 12700800 },
			{(double)-13289 / 1270080, (double)-4559 / 4233600, (double)15503 / 2540160, (double)-65033 / 12700800, (double)10547 / 2540160, (double)-26429 / 4233600, (double)21431 / 1814400 },
			{(double)-1783 / 132300, (double)-58231 / 6350400, (double)-14227 / 6350400, (double)61169 / 12700800, (double)-40651 / 12700800, (double)2627 / 907200, (double)-86999 / 12700800 },
			{(double)3 / 196, (double)661 / 105840, (double)-17 / 105840, (double)-17 / 8640, (double)319 / 235200, (double)-29 / 31360, (double)1433 / 1411200}
		};
			
		result.slice(5)={
			{(double)67 / 4233600, (double)-13 / 846720, (double)59 / 2116800, (double)-29 / 423360, (double)19 / 105840, (double)-23 / 105840, (double)-13 / 8820 },
			{(double)-2531 / 12700800, (double)17 / 254016, (double)-991 / 12700800, (double)2141 / 12700800, (double)-2767 / 6350400, (double)4097 / 6350400, (double)1259 / 529200 },
			{(double)2209 / 2540160, (double)-1007 / 4233600, (double)1591 / 12700800, (double)-2213 / 12700800, (double)1067 / 2540160, (double)-383 / 470400, (double)-193 / 907200 },
			{(double)-17509 / 12700800, (double)1039 / 1814400, (double)-607 / 2540160, (double)157 / 1058400, (double)-607 / 2540160, (double)1039 / 1814400, (double)-17509 / 12700800 },
			{(double)-193 / 907200, (double)-383 / 470400, (double)1067 / 2540160, (double)-2213 / 12700800, (double)1591 / 12700800, (double)-1007 / 4233600, (double)2209 / 2540160 },
			{(double)1259 / 529200, (double)4097 / 6350400, (double)-2767 / 6350400, (double)2141 / 12700800, (double)-991 / 12700800, (double)17 / 254016, (double)-2531 / 12700800 },
			{(double)-13 / 8820, (double)-23 / 105840, (double)19 / 105840, (double)-29 / 423360, (double)59 / 2116800, (double)-13 / 846720, (double)67 / 4233600 }
		};

		result.slice(6)={
			{(double)1 / 12700800, (double)-1 / 12700800, (double)1 / 6350400, (double)-1 / 2116800, (double)1 / 529200, (double)-1 / 105840, (double)1 / 17640 },
			{(double)-1 / 604800, (double)1 / 2116800, (double)-1 / 1814400, (double)17 / 12700800, (double)-31 / 6350400, (double)149 / 6350400, (double)-1 / 7200 },
			{(double)1 / 72576, (double)-1 / 362880, (double)1 / 846720, (double)-1 / 604800, (double)1 / 201600, (double)-281 / 12700800, (double)29 / 226800 },
			{(double)-1 / 17280, (double)19 / 1814400, (double)-1 / 362880, (double)1 / 635040, (double)-1 / 362880, (double)19 / 1814400, (double)-1 / 17280 },
			{(double)29 / 226800, (double)-281 / 12700800, (double)1 / 201600, (double)-1 / 604800, (double)1 / 846720, (double)-1 / 362880, (double)1 / 72576 },
			{(double)-1 / 7200, (double)149 / 6350400, (double)-31 / 6350400, (double)17 / 12700800, (double)-1 / 1814400, (double)1 / 2116800, (double)-1 / 604800 },
			{(double)1 / 17640, (double)-1 / 105840, (double)1 / 529200, (double)-1 / 2116800, (double)1 / 6350400, (double)-1 / 12700800, (double)1 / 12700800}
		};

		/*result.slice(1) << (double)529 / 882 << (double)4 / 2205 << (double)-437 / 3528 << (double)-1151 / 5880 << (double)-713 / 2940 << (double)-809 / 2940 << (double)-293 / 980 << endr
			<< (double)221 / 1764 << (double)2351 / 4410 << (double)113 / 17640 << (double)-107 / 1176 << (double)-3 / 20 << (double)-557 / 2940 << (double)-213 / 980 << endr
			<< (double)29 / 2940 << (double)209 / 2940 << (double)1166 / 2205 << (double)193 / 8820 << (double)-215 / 3528 << (double)-281 / 2520 << (double)-641 / 4410 << endr
			<< (double)-323 / 4410 << (double)-527 / 17640 << (double)743 / 17640 << (double)1168 / 2205 << (double)743 / 17640 << (double)-527 / 17640 << (double)-323 / 4410 << endr
			<< (double)-641 / 4410 << (double)-281 / 2520 << (double)-215 / 3528 << (double)193 / 8820 << (double)1166 / 2205 << (double)209 / 2940 << (double)29 / 2940 << endr
			<< (double)-213 / 980 << (double)-557 / 2940 << (double)-3 / 20 << (double)-107 / 1176 << (double)113 / 17640 << (double)2351 / 4410 << (double)221 / 1764 << endr
			<< (double)-293 / 980 << (double)-809 / 2940 << (double)-713 / 2940 << (double)-1151 / 5880 << (double)-437 / 3528 << (double)4 / 2205 << (double)529 / 882 << endr;
	
		result.slice(2) << (double)163103 / 705600 << (double)-262159 / 2116800 << (double)-38933 / 1058400 << (double)16901 / 352800 << (double)5237 / 44100 << (double)311 / 1764 << (double)3929 / 17640 << endr
			<< (double)-65759 / 352800 << (double)72697 / 264600 << (double)-222569 / 2116800 << (double)-48341 / 705600 << (double)-6523 / 352800 << (double)3223 / 117600 << (double)2593 / 39200 << endr
			<< (double)-9139 / 52920 << (double)-314971 / 2116800 << (double)63683 / 235200 << (double)-239119 / 2116800 << (double)-4387 / 43200 << (double)-151427 / 2116800 << (double)-87853 / 2116800 << endr
			<< (double)-63211 / 529200 << (double)-284233 / 2116800 << (double)-269987 / 2116800 << (double)17671 / 66150 << (double)-269987 / 2116800 << (double)-284233 / 2116800 << (double)-63211 / 529200 << endr
			<< (double)-87853 / 2116800 << (double)-151427 / 2116800 << (double)-4387 / 43200 << (double)-239119 / 2116800 << (double)63683 / 235200 << (double)-314971 / 2116800 << (double)-9139 / 52920 << endr
			<< (double)2593 / 39200 << (double)3223 / 117600 << (double)-6523 / 352800 << (double)-48341 / 705600 << (double)-222569 / 2116800 << (double)72697 / 264600 << (double)-65759 / 352800 << endr
			<< (double)3929 / 17640 << (double)311 / 1764 << (double)5237 / 44100 << (double)16901 / 352800 << (double)-38933 / 1058400 << (double)-262159 / 2116800 << (double)163103 / 705600 << endr;
		
		result.slice(3) << (double)319847 / 12700800 << (double)-50509 / 2540160 << (double)104207 / 6350400 << (double)991 / 141120 << (double)-289 / 15120 << (double)-5293 / 105840 << (double)-71 / 882 << endr
			<< (double)-37949 / 508032 << (double)58843 / 1270080 << (double)-19211 / 470400 << (double)49337 / 4233600 << (double)59933 / 2116800 << (double)58621 / 2116800 << (double)10361 / 529200 << endr
			<< (double)92443 / 12700800 << (double)-745111 / 12700800 << (double)677017 / 12700800 << (double)-590339 / 12700800 << (double)33293 / 2540160 << (double)534701 / 12700800 << (double)12829 / 235200 << endr
			<< (double)123481 / 2540160 << (double)17551 / 1411200 << (double)-216679 / 4233600 << (double)176569 / 3175200 << (double)-216679 / 4233600 << (double)17551 / 1411200 << (double)123481 / 2540160 << endr
			<< (double)12829 / 235200 << (double)534701 / 12700800 << (double)33293 / 2540160 << (double)-590339 / 12700800 << (double)677017 / 12700800 << (double)-745111 / 12700800 << (double)92443 / 12700800 << endr
			<< (double)10361 / 529200 << (double)58621 / 2116800 << (double)59933 / 2116800 << (double)49337 / 4233600 << (double)-19211 / 470400 << (double)58843 / 1270080 << (double)-37949 / 508032 << endr
			<< (double)-71 / 882 << (double)-5293 / 105840 << (double)-289 / 15120 << (double)991 / 141120 << (double)104207 / 6350400 << (double)-50509 / 2540160 << (double)319847 / 12700800 << endr;
		
		result.slice(4) << (double)1433 / 1411200 << (double)-29 / 31360 << (double)319 / 235200 << (double)-17 / 8640 << (double)-17 / 105840 << (double)661 / 105840 << (double)3 / 196 << endr
			<< (double)-86999 / 12700800 << (double)2627 / 907200 << (double)-40651 / 12700800 << (double)61169 / 12700800 << (double)-14227 / 6350400 << (double)-58231 / 6350400 << (double)-1783 / 132300 << endr
			<< (double)21431 / 1814400 << (double)-26429 / 4233600 << (double)10547 / 2540160 << (double)-65033 / 12700800 << (double)15503 / 2540160 << (double)-4559 / 4233600 << (double)-13289 / 1270080 << endr
			<< (double)33743 / 12700800 << (double)105073 / 12700800 << (double)-76331 / 12700800 << (double)229 / 50400 << (double)-76331 / 12700800 << (double)105073 / 12700800 << (double)33743 / 12700800 << endr
			<< (double)-13289 / 1270080 << (double)-4559 / 4233600 << (double)15503 / 2540160 << (double)-65033 / 12700800 << (double)10547 / 2540160 << (double)-26429 / 4233600 << (double)21431 / 1814400 << endr
			<< (double)-1783 / 132300 << (double)-58231 / 6350400 << (double)-14227 / 6350400 << (double)61169 / 12700800 << (double)-40651 / 12700800 << (double)2627 / 907200 << (double)-86999 / 12700800 << endr
			<< (double)3 / 196 << (double)661 / 105840 << (double)-17 / 105840 << (double)-17 / 8640 << (double)319 / 235200 << (double)-29 / 31360 << (double)1433 / 1411200 << endr;
	
		result.slice(5) << (double)67 / 4233600 << (double)-13 / 846720 << (double)59 / 2116800 << (double)-29 / 423360 << (double)19 / 105840 << (double)-23 / 105840 << (double)-13 / 8820 << endr
			<< (double)-2531 / 12700800 << (double)17 / 254016 << (double)-991 / 12700800 << (double)2141 / 12700800 << (double)-2767 / 6350400 << (double)4097 / 6350400 << (double)1259 / 529200 << endr
			<< (double)2209 / 2540160 << (double)-1007 / 4233600 << (double)1591 / 12700800 << (double)-2213 / 12700800 << (double)1067 / 2540160 << (double)-383 / 470400 << (double)-193 / 907200 << endr
			<< (double)-17509 / 12700800 << (double)1039 / 1814400 << (double)-607 / 2540160 << (double)157 / 1058400 << (double)-607 / 2540160 << (double)1039 / 1814400 << (double)-17509 / 12700800 << endr
			<< (double)-193 / 907200 << (double)-383 / 470400 << (double)1067 / 2540160 << (double)-2213 / 12700800 << (double)1591 / 12700800 << (double)-1007 / 4233600 << (double)2209 / 2540160 << endr
			<< (double)1259 / 529200 << (double)4097 / 6350400 << (double)-2767 / 6350400 << (double)2141 / 12700800 << (double)-991 / 12700800 << (double)17 / 254016 << (double)-2531 / 12700800 << endr
			<< (double)-13 / 8820 << (double)-23 / 105840 << (double)19 / 105840 << (double)-29 / 423360 << (double)59 / 2116800 << (double)-13 / 846720 << (double)67 / 4233600 << endr;

		result.slice(6) << (double)1 / 12700800 << (double)-1 / 12700800 << (double)1 / 6350400 << (double)-1 / 2116800 << (double)1 / 529200 << (double)-1 / 105840 << (double)1 / 17640 << endr
			<< (double)-1 / 604800 << (double)1 / 2116800 << (double)-1 / 1814400 << (double)17 / 12700800 << (double)-31 / 6350400 << (double)149 / 6350400 << (double)-1 / 7200 << endr
			<< (double)1 / 72576 << (double)-1 / 362880 << (double)1 / 846720 << (double)-1 / 604800 << (double)1 / 201600 << (double)-281 / 12700800 << (double)29 / 226800 << endr
			<< (double)-1 / 17280 << (double)19 / 1814400 << (double)-1 / 362880 << (double)1 / 635040 << (double)-1 / 362880 << (double)19 / 1814400 << (double)-1 / 17280 << endr
			<< (double)29 / 226800 << (double)-281 / 12700800 << (double)1 / 201600 << (double)-1 / 604800 << (double)1 / 846720 << (double)-1 / 362880 << (double)1 / 72576 << endr
			<< (double)-1 / 7200 << (double)149 / 6350400 << (double)-31 / 6350400 << (double)17 / 12700800 << (double)-1 / 1814400 << (double)1 / 2116800 << (double)-1 / 604800 << endr
			<< (double)1 / 17640 << (double)-1 / 105840 << (double)1 / 529200 << (double)-1 / 2116800 << (double)1 / 6350400 << (double)-1 / 12700800 << (double)1 / 12700800 << endr;*/

	}
	if (_m == 8) {
		result.slice(0) = result.slice(0) / 8;

		result.slice(1) = {
			{(double)46847 / 80640, (double)61 / 2304, (double)-7363 / 80640, (double)-857 / 5376, (double)-183 / 896, (double)-151 / 640, (double)-83 / 320, (double)-621 / 2240},
			{(double)2867 / 20160, (double)20737 / 40320, (double)103 / 4032, (double)-887 / 13440, (double)-3287 / 26880, (double)-719 / 4480, (double)-421 / 2240, (double)-2801 / 13440},
			{(double)223 / 5760, (double)1853 / 20160, (double)2281 / 4480, (double)143 / 4032, (double)-497 / 11520, (double)-29 / 315, (double)-5059 / 40320, (double)-1507 / 10080},
			{(double)-337 / 10080, (double)167 / 40320, (double)1331 / 20160, (double)41189 / 80640, (double)1969 / 40320, (double)-241 / 11520, (double)-1033 / 16128, (double)-7523 / 80640},
			{(double)-7523 / 80640, (double)-1033 / 16128, (double)-241 / 11520, (double)1969 / 40320, (double)41189 / 80640, (double)1331 / 20160, (double)167 / 40320, (double)-337 / 10080},
			{(double)-1507 / 10080, (double)-5059 / 40320, (double)-29 / 315, (double)-497 / 11520, (double)143 / 4032, (double)2281 / 4480, (double)1853 / 20160, (double)223 / 5760},
			{(double)-2801 / 13440, (double)-421 / 2240, (double)-719 / 4480, (double)-3287 / 26880, (double)-887 / 13440, (double)103 / 4032, (double)20737 / 40320, (double)2867 / 20160},
			{(double)-621 / 2240, (double)-83 / 320, (double)-151 / 640, (double)-183 / 896, (double)-857 / 5376, (double)-7363 / 80640, (double)61 / 2304, (double)46847 / 80640}
		};

		result.slice(2) = {
			{(double)17506193 / 67737600, (double)-8447843 / 67737600, (double)-882323 / 16934400, (double)45317 / 1881600, (double)507971 / 5644800, (double)163573 / 1128960, (double)21475 / 112896, (double)18353 / 80640},
			{(double)-11548657 / 67737600, (double)5047349 / 16934400, (double)-6712847 / 67737600, (double)-4821647 / 67737600, (double)-661 / 25088, (double)10553 / 627200, (double)12827 / 235200, (double)54197 / 627200},
			{(double)-11635229 / 67737600, (double)-233903 / 1693440, (double)473341 / 1612800, (double)-3410417 / 33868800, (double)-6383177 / 67737600, (double)-10093 / 151200, (double)-158413 / 4233600, (double)-184861 / 16934400},
			{(double)-9251699 / 67737600, (double)-9544589 / 67737600, (double)-8158961 / 67737600, (double)3891973 / 13547520, (double)-819619 / 7526400, (double)-3955831 / 33868800, (double)-1721353 / 16934400, (double)-17503 / 211680},
			{(double)-17503 / 211680, (double)-1721353 / 16934400, (double)-3955831 / 33868800, (double)-819619 / 7526400, (double)3891973 / 13547520, (double)-8158961 / 67737600, (double)-9544589 / 67737600, (double)-9251699 / 67737600},
			{(double)-184861 / 16934400, (double)-158413 / 4233600, (double)-10093 / 151200, (double)-6383177 / 67737600, (double)-3410417 / 33868800, (double)473341 / 1612800, (double)-233903 / 1693440, (double)-11635229 / 67737600},
			{(double)54197 / 627200, (double)12827 / 235200, (double)10553 / 627200, (double)-661 / 25088, (double)-4821647 / 67737600, (double)-6712847 / 67737600, (double)5047349 / 16934400, (double)-11548657 / 67737600},
			{(double)18353 / 80640, (double)21475 / 112896, (double)163573 / 1128960, (double)507971 / 5644800, (double)45317 / 1881600, (double)-882323 / 16934400, (double)-8447843 / 67737600, (double)17506193 / 67737600}
		};

		result.slice(3) = {
			{(double)366661 / 10838016, (double)-758431 / 30105600, (double)739091 / 45158400, (double)189037 / 15052800, (double)-11779 / 1254400, (double)-28351 / 752640, (double)-8399 / 125440, (double)-15289 / 161280},
			{(double)-17495893 / 203212800, (double)47135581 / 812851200, (double)-18921347 / 406425600, (double)4929857 / 812851200, (double)1664437 / 67737600, (double)3450997 / 135475200, (double)1223423 / 67737600, (double)10813 / 1505280},
			{(double)-2748883 / 406425600, (double)-28078793 / 406425600, (double)53498239 / 812851200, (double)-10629557 / 203212800, (double)1013629 / 203212800, (double)9158533 / 270950400, (double)6254159 / 135475200, (double)105533 / 2116800},
			{(double)1588373 / 40642560, (double)107129 / 203212800, (double)-12406927 / 203212800, (double)56509741 / 812851200, (double)-2845271 / 50803200, (double)388253 / 101606400, (double)79823 / 2073600, (double)46907171 / 812851200},
			{(double)46907171 / 812851200, (double)79823 / 2073600, (double)388253 / 101606400, (double)-2845271 / 50803200, (double)56509741 / 812851200, (double)-12406927 / 203212800, (double)107129 / 203212800, (double)1588373 / 40642560},
			{(double)105533 / 2116800, (double)6254159 / 135475200, (double)9158533 / 270950400, (double)1013629 / 203212800, (double)-10629557 / 203212800, (double)53498239 / 812851200, (double)-28078793 / 406425600, (double)-2748883 / 406425600},
			{(double)10813 / 1505280, (double)1223423 / 67737600, (double)3450997 / 135475200, (double)1664437 / 67737600, (double)4929857 / 812851200, (double)-18921347 / 406425600, (double)47135581 / 812851200, (double)-17495893 / 203212800},
			{(double)-15289 / 161280, (double)-8399 / 125440, (double)-28351 / 752640, (double)-11779 / 1254400, (double)189037 / 15052800, (double)739091 / 45158400, (double)-758431 / 30105600, (double)366661 / 10838016}
		};

		result.slice(4) = {
			{(double)710449 / 406425600, (double)-155623 / 101606400, (double)80483 / 40642560, (double)-47429 / 22579200, (double)-3859 / 2419200, (double)2501 / 677376, (double)685 / 56448, (double)179 / 8064},
			{(double)-1035371 / 101606400, (double)609631 / 135475200, (double)-1918663 / 406425600, (double)2353697 / 406425600, (double)-899 / 2822400, (double)-159883 / 22579200, (double)-193363 / 16934400, (double)-18841 / 1411200},
			{(double)1347851 / 101606400, (double)-3744953 / 406425600, (double)2551259 / 406425600, (double)-361507 / 50803200, (double)1361537 / 203212800, (double)49799 / 58060800, (double)-812699 / 101606400, (double)-3202991 / 203212800},
			{(double)510847 / 67737600, (double)1027969 / 101606400, (double)-405991 / 45158400, (double)576299 / 81285120, (double)-21467 / 2540160, (double)1622597 / 203212800, (double)173273 / 50803200, (double)-738287 / 135475200},
			{(double)-738287 / 135475200, (double)173273 / 50803200, (double)1622597 / 203212800, (double)-21467 / 2540160, (double)576299 / 81285120, (double)-405991 / 45158400, (double)1027969 / 101606400, (double)510847 / 67737600},
			{(double)-3202991 / 203212800, (double)-812699 / 101606400, (double)49799 / 58060800, (double)1361537 / 203212800, (double)-361507 / 50803200, (double)2551259 / 406425600, (double)-3744953 / 406425600, (double)1347851 / 101606400},
			{(double)-18841 / 1411200, (double)-193363 / 16934400, (double)-159883 / 22579200, (double)-899 / 2822400, (double)2353697 / 406425600, (double)-1918663 / 406425600, (double)609631 / 135475200, (double)-1035371 / 101606400},
			{(double)179 / 8064, (double)685 / 56448, (double)2501 / 677376, (double)-3859 / 2419200, (double)-47429 / 22579200, (double)80483 / 40642560, (double)-155623 / 101606400, (double)710449 / 406425600}
		};

		result.slice(5) = {
			{(double)1331 / 33868800, (double)-5027 / 135475200, (double)419 / 6773760, (double)-619 / 4838400, (double)979 / 4233600, (double)19 / 423360, (double)-593 / 564480, (double)-239 / 80640},
			{(double)-1087 / 2540160, (double)29791 / 203212800, (double)-67589 / 406425600, (double)64439 / 203212800, (double)-61627 / 101606400, (double)5693 / 25401600, (double)62687 / 33868800, (double)39961 / 11289600},
			{(double)43691 / 29030400, (double)-65029 / 135475200, (double)52963 / 203212800, (double)-8683 / 25401600, (double)19513 / 29030400, (double)-1789 / 2257920, (double)-78599 / 203212800, (double)82441 / 67737600},
			{(double)-297403 / 203212800, (double)109973 / 101606400, (double)-209843 / 406425600, (double)16459 / 50803200, (double)-15857 / 33868800, (double)8549 / 9676800, (double)-9529 / 8467200, (double)-293801 / 203212800},
			{(double)-293801 / 203212800, (double)-9529 / 8467200, (double)8549 / 9676800, (double)-15857 / 33868800, (double)16459 / 50803200, (double)-209843 / 406425600, (double)109973 / 101606400, (double)-297403 / 203212800},
			{(double)82441 / 67737600, (double)-78599 / 203212800, (double)-1789 / 2257920, (double)19513 / 29030400, (double)-8683 / 25401600, (double)52963 / 203212800, (double)-65029 / 135475200, (double)43691 / 29030400},
			{(double)39961 / 11289600, (double)62687 / 33868800, (double)5693 / 25401600, (double)-61627 / 101606400, (double)64439 / 203212800, (double)-67589 / 406425600, (double)29791 / 203212800, (double)-1087 / 2540160},
			{(double)-239 / 80640, (double)-593 / 564480, (double)19 / 423360, (double)979 / 4233600, (double)-619 / 4838400, (double)419 / 6773760, (double)-5027 / 135475200, (double)1331 / 33868800}
		};

		result.slice(6) = {
			{(double)17 / 45158400, (double)-1 / 2709504, (double)47 / 67737600, (double)-5 / 2709504, (double)1 / 169344, (double)-1 / 52920, (double)1 / 37632, (double)17 / 80640},
			{(double)-1387 / 203212800, (double)53 / 27095040, (double)-61 / 27095040, (double)683 / 135475200, (double)-773 / 50803200, (double)1979 / 40642560, (double)-1381 / 16934400, (double)-2209 / 5644800},
			{(double)9469 / 203212800, (double)-4099 / 406425600, (double)611 / 135475200, (double)-23 / 3763200, (double)11 / 691200, (double)-2273 / 45158400, (double)10931 / 101606400, (double)1709 / 13547520},
			{(double)-3637 / 25401600, (double)7397 / 203212800, (double)-527 / 45158400, (double)173 / 27095040, (double)-2047 / 203212800, (double)1189 / 40642560, (double)-409 / 5080320, (double)9139 / 58060800},
			{(double)9139 / 58060800, (double)-409 / 5080320, (double)1189 / 40642560, (double)-2047 / 203212800, (double)173 / 27095040, (double)-527 / 45158400, (double)7397 / 203212800, (double)-3637 / 25401600},
			{(double)1709 / 13547520, (double)10931 / 101606400, (double)-2273 / 45158400, (double)11 / 691200, (double)-23 / 3763200, (double)611 / 135475200, (double)-4099 / 406425600, (double)9469 / 203212800},
			{(double)-2209 / 5644800, (double)-1381 / 16934400, (double)1979 / 40642560, (double)-773 / 50803200, (double)683 / 135475200, (double)-61 / 27095040, (double)53 / 27095040, (double)-1387 / 203212800},
			{(double)17 / 80640, (double)1 / 37632, (double)-1 / 52920, (double)1 / 169344, (double)-5 / 2709504, (double)47 / 67737600, (double)-1 / 2709504, (double)17 / 45158400}
		};

		result.slice(7) = {
			{(double)1 / 812851200, (double)-1 / 812851200, (double)1 / 406425600, (double)-1 / 135475200, (double)1 / 33868800, (double)-1 / 6773760, (double)1 / 1128960, (double)-1 / 161280},
			{(double)-1 / 29030400, (double)1 / 116121600, (double)-1 / 101606400, (double)19 / 812851200, (double)-17 / 203212800, (double)23 / 58060800, (double)-157 / 67737600, (double)121 / 7526400},
			{(double)23 / 58060800, (double)-1 / 14515200, (double)1 / 38707200, (double)-1 / 29030400, (double)1 / 10160640, (double)-7 / 16588800, (double)961 / 406425600, (double)-67 / 4147200},
			{(double)-1 / 414720, (double)11 / 29030400, (double)-1 / 11612160, (double)1 / 23224320, (double)-1 / 14515200, (double)1 / 4147200, (double)-127 / 101606400, (double)967 / 116121600},
			{(double)967 / 116121600, (double)-127 / 101606400, (double)1 / 4147200, (double)-1 / 14515200, (double)1 / 23224320, (double)-1 / 11612160, (double)11 / 29030400, (double)-1 / 414720},
			{(double)-67 / 4147200, (double)961 / 406425600, (double)-7 / 16588800, (double)1 / 10160640, (double)-1 / 29030400, (double)1 / 38707200, (double)-1 / 14515200, (double)23 / 58060800},
			{(double)121 / 7526400, (double)-157 / 67737600, (double)23 / 58060800, (double)-17 / 203212800, (double)19 / 812851200, (double)-1 / 101606400, (double)1 / 116121600, (double)-1 / 29030400},
			{(double)-1 / 161280, (double)1 / 1128960, (double)-1 / 6773760, (double)1 / 33868800, (double)-1 / 135475200, (double)1 / 406425600, (double)-1 / 812851200, (double)1 / 812851200}
		};

		/*result.slice(1) << (double)46847 / 80640 << (double)61 / 2304 << (double)-7363 / 80640 << (double)-857 / 5376 << (double)-183 / 896 << (double)-151 / 640 << (double)-83 / 320 << (double)-621 / 2240 << endr
			<< (double)2867 / 20160 << (double)20737 / 40320 << (double)103 / 4032 << (double)-887 / 13440 << (double)-3287 / 26880 << (double)-719 / 4480 << (double)-421 / 2240 << (double)-2801 / 13440 << endr
			<< (double)223 / 5760 << (double)1853 / 20160 << (double)2281 / 4480 << (double)143 / 4032 << (double)-497 / 11520 << (double)-29 / 315 << (double)-5059 / 40320 << (double)-1507 / 10080 << endr
			<< (double)-337 / 10080 << (double)167 / 40320 << (double)1331 / 20160 << (double)41189 / 80640 << (double)1969 / 40320 << (double)-241 / 11520 << (double)-1033 / 16128 << (double)-7523 / 80640 << endr
			<< (double)-7523 / 80640 << (double)-1033 / 16128 << (double)-241 / 11520 << (double)1969 / 40320 << (double)41189 / 80640 << (double)1331 / 20160 << (double)167 / 40320 << (double)-337 / 10080 << endr
			<< (double)-1507 / 10080 << (double)-5059 / 40320 << (double)-29 / 315 << (double)-497 / 11520 << (double)143 / 4032 << (double)2281 / 4480 << (double)1853 / 20160 << (double)223 / 5760 << endr
			<< (double)-2801 / 13440 << (double)-421 / 2240 << (double)-719 / 4480 << (double)-3287 / 26880 << (double)-887 / 13440 << (double)103 / 4032 << (double)20737 / 40320 << (double)2867 / 20160 << endr
			<< (double)-621 / 2240 << (double)-83 / 320 << (double)-151 / 640 << (double)-183 / 896 << (double)-857 / 5376 << (double)-7363 / 80640 << (double)61 / 2304 << (double)46847 / 80640 << endr;

		result.slice(2) << (double)17506193 / 67737600 << (double)-8447843 / 67737600 << (double)-882323 / 16934400 << (double)45317 / 1881600 << (double)507971 / 5644800 << (double)163573 / 1128960 << (double)21475 / 112896 << (double)18353 / 80640 << endr
			<< (double)-11548657 / 67737600 << (double)5047349 / 16934400 << (double)-6712847 / 67737600 << (double)-4821647 / 67737600 << (double)-661 / 25088 << (double)10553 / 627200 << (double)12827 / 235200 << (double)54197 / 627200 << endr
			<< (double)-11635229 / 67737600 << (double)-233903 / 1693440 << (double)473341 / 1612800 << (double)-3410417 / 33868800 << (double)-6383177 / 67737600 << (double)-10093 / 151200 << (double)-158413 / 4233600 << (double)-184861 / 16934400 << endr
			<< (double)-9251699 / 67737600 << (double)-9544589 / 67737600 << (double)-8158961 / 67737600 << (double)3891973 / 13547520 << (double)-819619 / 7526400 << (double)-3955831 / 33868800 << (double)-1721353 / 16934400 << (double)-17503 / 211680 << endr
			<< (double)-17503 / 211680 << (double)-1721353 / 16934400 << (double)-3955831 / 33868800 << (double)-819619 / 7526400 << (double)3891973 / 13547520 << (double)-8158961 / 67737600 << (double)-9544589 / 67737600 << (double)-9251699 / 67737600 << endr
			<< (double)-184861 / 16934400 << (double)-158413 / 4233600 << (double)-10093 / 151200 << (double)-6383177 / 67737600 << (double)-3410417 / 33868800 << (double)473341 / 1612800 << (double)-233903 / 1693440 << (double)-11635229 / 67737600 << endr
			<< (double)54197 / 627200 << (double)12827 / 235200 << (double)10553 / 627200 << (double)-661 / 25088 << (double)-4821647 / 67737600 << (double)-6712847 / 67737600 << (double)5047349 / 16934400 << (double)-11548657 / 67737600 << endr
			<< (double)18353 / 80640 << (double)21475 / 112896 << (double)163573 / 1128960 << (double)507971 / 5644800 << (double)45317 / 1881600 << (double)-882323 / 16934400 << (double)-8447843 / 67737600 << (double)17506193 / 67737600 << endr;

		result.slice(3) << (double)366661 / 10838016 << (double)-758431 / 30105600 << (double)739091 / 45158400 << (double)189037 / 15052800 << (double)-11779 / 1254400 << (double)-28351 / 752640 << (double)-8399 / 125440 << (double)-15289 / 161280 << endr
			<< (double)-17495893 / 203212800 << (double)47135581 / 812851200 << (double)-18921347 / 406425600 << (double)4929857 / 812851200 << (double)1664437 / 67737600 << (double)3450997 / 135475200 << (double)1223423 / 67737600 << (double)10813 / 1505280 << endr
			<< (double)-2748883 / 406425600 << (double)-28078793 / 406425600 << (double)53498239 / 812851200 << (double)-10629557 / 203212800 << (double)1013629 / 203212800 << (double)9158533 / 270950400 << (double)6254159 / 135475200 << (double)105533 / 2116800 << endr
			<< (double)1588373 / 40642560 << (double)107129 / 203212800 << (double)-12406927 / 203212800 << (double)56509741 / 812851200 << (double)-2845271 / 50803200 << (double)388253 / 101606400 << (double)79823 / 2073600 << (double)46907171 / 812851200 << endr
			<< (double)46907171 / 812851200 << (double)79823 / 2073600 << (double)388253 / 101606400 << (double)-2845271 / 50803200 << (double)56509741 / 812851200 << (double)-12406927 / 203212800 << (double)107129 / 203212800 << (double)1588373 / 40642560 << endr
			<< (double)105533 / 2116800 << (double)6254159 / 135475200 << (double)9158533 / 270950400 << (double)1013629 / 203212800 << (double)-10629557 / 203212800 << (double)53498239 / 812851200 << (double)-28078793 / 406425600 << (double)-2748883 / 406425600 << endr
			<< (double)10813 / 1505280 << (double)1223423 / 67737600 << (double)3450997 / 135475200 << (double)1664437 / 67737600 << (double)4929857 / 812851200 << (double)-18921347 / 406425600 << (double)47135581 / 812851200 << (double)-17495893 / 203212800 << endr
			<< (double)-15289 / 161280 << (double)-8399 / 125440 << (double)-28351 / 752640 << (double)-11779 / 1254400 << (double)189037 / 15052800 << (double)739091 / 45158400 << (double)-758431 / 30105600 << (double)366661 / 10838016 << endr;

		result.slice(4) << (double)710449 / 406425600 << (double)-155623 / 101606400 << (double)80483 / 40642560 << (double)-47429 / 22579200 << (double)-3859 / 2419200 << (double)2501 / 677376 << (double)685 / 56448 << (double)179 / 8064 << endr
			<< (double)-1035371 / 101606400 << (double)609631 / 135475200 << (double)-1918663 / 406425600 << (double)2353697 / 406425600 << (double)-899 / 2822400 << (double)-159883 / 22579200 << (double)-193363 / 16934400 << (double)-18841 / 1411200 << endr
			<< (double)1347851 / 101606400 << (double)-3744953 / 406425600 << (double)2551259 / 406425600 << (double)-361507 / 50803200 << (double)1361537 / 203212800 << (double)49799 / 58060800 << (double)-812699 / 101606400 << (double)-3202991 / 203212800 << endr
			<< (double)510847 / 67737600 << (double)1027969 / 101606400 << (double)-405991 / 45158400 << (double)576299 / 81285120 << (double)-21467 / 2540160 << (double)1622597 / 203212800 << (double)173273 / 50803200 << (double)-738287 / 135475200 << endr
			<< (double)-738287 / 135475200 << (double)173273 / 50803200 << (double)1622597 / 203212800 << (double)-21467 / 2540160 << (double)576299 / 81285120 << (double)-405991 / 45158400 << (double)1027969 / 101606400 << (double)510847 / 67737600 << endr
			<< (double)-3202991 / 203212800 << (double)-812699 / 101606400 << (double)49799 / 58060800 << (double)1361537 / 203212800 << (double)-361507 / 50803200 << (double)2551259 / 406425600 << (double)-3744953 / 406425600 << (double)1347851 / 101606400 << endr
			<< (double)-18841 / 1411200 << (double)-193363 / 16934400 << (double)-159883 / 22579200 << (double)-899 / 2822400 << (double)2353697 / 406425600 << (double)-1918663 / 406425600 << (double)609631 / 135475200 << (double)-1035371 / 101606400 << endr
			<< (double)179 / 8064 << (double)685 / 56448 << (double)2501 / 677376 << (double)-3859 / 2419200 << (double)-47429 / 22579200 << (double)80483 / 40642560 << (double)-155623 / 101606400 << (double)710449 / 406425600 << endr;
	
		result.slice(5) << (double)1331 / 33868800 << (double)-5027 / 135475200 << (double)419 / 6773760 << (double)-619 / 4838400 << (double)979 / 4233600 << (double)19 / 423360 << (double)-593 / 564480 << (double)-239 / 80640 << endr
			<< (double)-1087 / 2540160 << (double)29791 / 203212800 << (double)-67589 / 406425600 << (double)64439 / 203212800 << (double)-61627 / 101606400 << (double)5693 / 25401600 << (double)62687 / 33868800 << (double)39961 / 11289600 << endr
			<< (double)43691 / 29030400 << (double)-65029 / 135475200 << (double)52963 / 203212800 << (double)-8683 / 25401600 << (double)19513 / 29030400 << (double)-1789 / 2257920 << (double)-78599 / 203212800 << (double)82441 / 67737600 << endr
			<< (double)-297403 / 203212800 << (double)109973 / 101606400 << (double)-209843 / 406425600 << (double)16459 / 50803200 << (double)-15857 / 33868800 << (double)8549 / 9676800 << (double)-9529 / 8467200 << (double)-293801 / 203212800 << endr
			<< (double)-293801 / 203212800 << (double)-9529 / 8467200 << (double)8549 / 9676800 << (double)-15857 / 33868800 << (double)16459 / 50803200 << (double)-209843 / 406425600 << (double)109973 / 101606400 << (double)-297403 / 203212800 << endr
			<< (double)82441 / 67737600 << (double)-78599 / 203212800 << (double)-1789 / 2257920 << (double)19513 / 29030400 << (double)-8683 / 25401600 << (double)52963 / 203212800 << (double)-65029 / 135475200 << (double)43691 / 29030400 << endr
			<< (double)39961 / 11289600 << (double)62687 / 33868800 << (double)5693 / 25401600 << (double)-61627 / 101606400 << (double)64439 / 203212800 << (double)-67589 / 406425600 << (double)29791 / 203212800 << (double)-1087 / 2540160 << endr
			<< (double)-239 / 80640 << (double)-593 / 564480 << (double)19 / 423360 << (double)979 / 4233600 << (double)-619 / 4838400 << (double)419 / 6773760 << (double)-5027 / 135475200 << (double)1331 / 33868800 << endr;

		result.slice(6) << (double)17 / 45158400 << (double)-1 / 2709504 << (double)47 / 67737600 << (double)-5 / 2709504 << (double)1 / 169344 << (double)-1 / 52920 << (double)1 / 37632 << (double)17 / 80640 << endr
			<< (double)-1387 / 203212800 << (double)53 / 27095040 << (double)-61 / 27095040 << (double)683 / 135475200 << (double)-773 / 50803200 << (double)1979 / 40642560 << (double)-1381 / 16934400 << (double)-2209 / 5644800 << endr
			<< (double)9469 / 203212800 << (double)-4099 / 406425600 << (double)611 / 135475200 << (double)-23 / 3763200 << (double)11 / 691200 << (double)-2273 / 45158400 << (double)10931 / 101606400 << (double)1709 / 13547520 << endr
			<< (double)-3637 / 25401600 << (double)7397 / 203212800 << (double)-527 / 45158400 << (double)173 / 27095040 << (double)-2047 / 203212800 << (double)1189 / 40642560 << (double)-409 / 5080320 << (double)9139 / 58060800 << endr
			<< (double)9139 / 58060800 << (double)-409 / 5080320 << (double)1189 / 40642560 << (double)-2047 / 203212800 << (double)173 / 27095040 << (double)-527 / 45158400 << (double)7397 / 203212800 << (double)-3637 / 25401600 << endr
			<< (double)1709 / 13547520 << (double)10931 / 101606400 << (double)-2273 / 45158400 << (double)11 / 691200 << (double)-23 / 3763200 << (double)611 / 135475200 << (double)-4099 / 406425600 << (double)9469 / 203212800 << endr
			<< (double)-2209 / 5644800 << (double)-1381 / 16934400 << (double)1979 / 40642560 << (double)-773 / 50803200 << (double)683 / 135475200 << (double)-61 / 27095040 << (double)53 / 27095040 << (double)-1387 / 203212800 << endr
			<< (double)17 / 80640 << (double)1 / 37632 << (double)-1 / 52920 << (double)1 / 169344 << (double)-5 / 2709504 << (double)47 / 67737600 << (double)-1 / 2709504 << (double)17 / 45158400 << endr;

		result.slice(7) << (double)1 / 812851200 << (double)-1 / 812851200 << (double)1 / 406425600 << (double)-1 / 135475200 << (double)1 / 33868800 << (double)-1 / 6773760 << (double)1 / 1128960 << (double)-1 / 161280 << endr
			<< (double)-1 / 29030400 << (double)1 / 116121600 << (double)-1 / 101606400 << (double)19 / 812851200 << (double)-17 / 203212800 << (double)23 / 58060800 << (double)-157 / 67737600 << (double)121 / 7526400 << endr
			<< (double)23 / 58060800 << (double)-1 / 14515200 << (double)1 / 38707200 << (double)-1 / 29030400 << (double)1 / 10160640 << (double)-7 / 16588800 << (double)961 / 406425600 << (double)-67 / 4147200 << endr
			<< (double)-1 / 414720 << (double)11 / 29030400 << (double)-1 / 11612160 << (double)1 / 23224320 << (double)-1 / 14515200 << (double)1 / 4147200 << (double)-127 / 101606400 << (double)967 / 116121600 << endr
			<< (double)967 / 116121600 << (double)-127 / 101606400 << (double)1 / 4147200 << (double)-1 / 14515200 << (double)1 / 23224320 << (double)-1 / 11612160 << (double)11 / 29030400 << (double)-1 / 414720 << endr
			<< (double)-67 / 4147200 << (double)961 / 406425600 << (double)-7 / 16588800 << (double)1 / 10160640 << (double)-1 / 29030400 << (double)1 / 38707200 << (double)-1 / 14515200 << (double)23 / 58060800 << endr
			<< (double)121 / 7526400 << (double)-157 / 67737600 << (double)23 / 58060800 << (double)-17 / 203212800 << (double)19 / 812851200 << (double)-1 / 101606400 << (double)1 / 116121600 << (double)-1 / 29030400 << endr
			<< (double)-1 / 161280 << (double)1 / 1128960 << (double)-1 / 6773760 << (double)1 / 33868800 << (double)-1 / 135475200 << (double)1 / 406425600 << (double)-1 / 812851200 << (double)1 / 812851200 << endr;*/
	}


	return result;
}

void Bos::printResults() {
	// cout << "------printing BOS ------" << endl;
	this->_resmus.print();
	this->_respis.print();
}

List Bos::returnResults() {

	List bosRes = List::create(Rcpp::Named("mus") = _mus,
							   Rcpp::Named("pis") = _pis);
	return(bosRes);
}

List Bos::returnParamsChain() {

	List bosRes = List::create(Rcpp::Named("mus") = _allmus,
							   Rcpp::Named("pis") = _allpis);
	return(bosRes);
}

void Bos::putParamsToZero() {
	return;
}

double Bos::computeICL(int i, int d, int k, int h) {
	double result = 0;
	if(i==0 && d==0 && k==0 && h==0){
		result = - _kc*_kr * log(_Nr*_Jc);
	}
	cube cubeProbs = this->getCubeProbs();
	double accum = cubeProbs(k, h, _xsep(i, d) - 1);
	double toadd = log(accum);
	if(isnan(toadd)) 
		toadd= log(abs(accum));
	result += toadd;
	
	return(result);
}

double Bos::logsum(rowvec logx) {
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

Mupi Bos::ordiemCpp(
	const arma::colvec& x,
	const arma::colvec& tabmu0,
	const arma::colvec& tabp0,
	double eps = 1,
	int iter_max = 100)
{

	double ml = -inf;
	double p_ml = tabp0(0);
	int mu_ml = tabmu0(0);

	int n = x.n_rows;
	arma::colvec w(n); // in code, it is as a paramter, but never used as paramteer
	w.ones();


	double ntot = arma::sum(w);

	int a = 0;

	for (unsigned int imu = 0; imu < tabmu0.n_elem; ++imu) {
		int mu = tabmu0(imu);
		double mlold = -inf;
		for (unsigned int ip = 0; ip < tabp0.n_elem; ++ip) {
			double p = tabp0(ip);
			int nostop = 1;
			int iter = 0;
			while (nostop) {
				iter++;
				// -- E step ---
				// first: compute px for each modality
				arma::colvec pallx(_m);
				pallx.zeros();

				arma::colvec px(n);
				px.zeros();

				for (int i = 0; i<_m; ++i) {
					arma::colvec tmp(_m);
					tmp = tmp.ones() * p;
					for (int deg = 0; deg < _m; ++deg) {
						tmp(deg) = std::pow(tmp(deg), deg);
					}
					//Q.subcube( first_row, first_col, first_slice, last_row, last_col, last_slice )
					arma::rowvec sub_cube =  this->_tab_pejs.subcube(i, imu, 0, i, imu, this->_tab_pejs.n_slices - 1);
					

					arma::vec vecprod = trans(sub_cube) % tmp;
					pallx(i) = arma::sum(vecprod);

					arma::uvec ids = arma::find(x == (i + 1));

					px.elem(ids).fill(std::max(1e-300, pallx(i)));
				}

				arma::vec vecprod = w % arma::log(px);

				double mlnew = arma::sum(vecprod);

				// first: compute pxz1 for each modality
				arma::mat pallxz1 = arma::mat(_m, _m - 1);
				pallxz1.zeros();


				for (int i = 0; i < _m; ++i) {
					for (int j = 0; j < (_m - 1); ++j) {
						arma::colvec z1tozmm1 = arma::colvec(_m - 1);
						z1tozmm1.zeros();
						z1tozmm1(j) = 1;
						arma::urowvec ivec; // to gain time : put it out of this loop
						ivec << (i + 1);

						double pejv = pej(ivec, _m, _m, mu, p, z1tozmm1); // TODO change that
						pallxz1(i, j) = (double)(std::max(1e-300, pejv)) / std::max(1e-300, pallx(i));
					}
				}

				//second: affect each modality value to the corresponding units
				arma::mat pxz1 = arma::mat(n, _m - 1);
				pxz1.zeros();
				for (int i = 0; i < _m; ++i) {
					arma::uvec whereisi = arma::find(x == (i + 1));
					arma::uvec all_cols_pxz1 = arma::linspace<arma::uvec>(0, pxz1.n_cols - 1, pxz1.n_cols);

					int sumwhereisi = whereisi.size();
					arma::colvec matsumwhereisi = arma::mat(sumwhereisi, 1);
					matsumwhereisi.ones();


					arma::rowvec pallxz1i = pallxz1.row(i);
					arma::mat sub_mat = matsumwhereisi * pallxz1i;
					pxz1.submat(whereisi, all_cols_pxz1) = sub_mat; // 

				}
				arma::rowvec temp1(_m-1);
				temp1.ones();
				arma::mat temp2 = w * temp1;
				pxz1 = temp2 % pxz1;

				// ---- M step ----
				double sum = arma::accu(pxz1);
				double pmean = (double)sum / (ntot*(_m - 1));
				//// cout << pmean << std::endl;
				if (!(mlnew == -inf)) {
					if ((std::abs(mlnew - mlold) / ntot < eps) || (iter>(iter_max - 1))) {
						nostop = 0;
						if (mlnew > ml) {
							ml = mlnew;
							p_ml = pmean;
							mu_ml = mu;
						}
					}
				}
				else {
					if (iter>(iter_max - 1)) {
						nostop = 0;
						if (mlnew > ml) {
							ml = mlnew;
							p_ml = pmean;
							mu_ml = mu;
						}
					}
				}
				mlold = mlnew;

			}// while
		} // p
		++a;
	}//mu

	//cout << mu_ml << "  " << p_ml << "  " << ml << endl;

	/*List result = List::create(mu_ml, p_ml, ml);
	std::string result = "ok";
	cout << mu_ml << "  " << p_ml << "  " << ml << endl;
	return(result);*/
	Mupi result(p_ml, mu_ml);
	return result;
}



int Bos::unsigned_to_signed(unsigned x)
{
	if (x <= INT_MAX)
		return static_cast<int>(x);

	if (x >= INT_MIN)
		return static_cast<int>(x - INT_MIN) + INT_MIN;

	throw x; // Or whatever you like
}

bool Bos::compare_vec(arma::urowvec vec1, arma::rowvec vec2) {
	bool result = false;
	// OLD
	//bool resold = true;
	// if (vec1.size() != vec2.size()) {
	// 	resold = false;
	// }
	// else {
	// 	for (unsigned int i = 0; i < vec1.size(); ++i) {
	// 		int signed_cast = unsigned_to_signed(vec1(i));
	// 		if (signed_cast != vec2(i)) {
	// 			resold = false;
	// 			break;
	// 		}
	// 	}
	// }
	// END OLD
	if( all(vec1 == vec2 ) == 1){
		result = true;
	}
	return(result);
}


arma::umat Bos::allej(int j, int m) { // TODO . to implement
	arma::umat result;
	if (j == 1) {
		//result << 1 << m << arma::endr;
		result = {{(unsigned int)1,(unsigned int)m}};
		return(result);
	}
	if (j > m) {
		return(result);
	}
	arma::rowvec indicesi = arma::linspace<arma::rowvec>(1, m - j + 1, m - j + 1);
	for (unsigned int i = 0; i < indicesi.size(); ++i) {
		int sizeej = indicesi(i);
		arma::rowvec indicesj = arma::linspace<arma::rowvec>(1, m - sizeej + 1, m - sizeej + 1);
		for (unsigned int j = 0; j < indicesj.size(); ++j) {
			int binf = indicesj(j);
			int bsup = binf + sizeej - 1;
			arma::umat rowvector = {{(unsigned int)binf, (unsigned int)bsup}};
			//rowvector << binf << bsup << arma::endr;
			result = arma::join_vert(result, rowvector);
		}
	}
	return(result);
}

double Bos::pejp1_yjej(arma::urowvec ejp1, int yj, arma::urowvec ej, int mu, double p) {
	double proba = 0;
	arma::rowvec ejminus;
	ejminus << ej(0) << yj - 1;
	arma::rowvec ejequal;
	ejequal << yj << yj;
	arma::rowvec ejplus;
	ejplus << yj + 1 << ej(1);

	//pejp1_yjejzj0
	double pejp1_yjejzj0 = 0;
	if (compare_vec(ejp1, ejminus) || compare_vec(ejp1, ejequal) || compare_vec(ejp1, ejplus)) {
		pejp1_yjejzj0 = (double)(ejp1(1) - ejp1(0) + 1) / (ej(1) - ej(0) + 1);
	}
	//pejp1_yjejzj1
	double dmuejminus = 0;
	if (ejminus(0) > ejminus(1)) {
		dmuejminus = inf;
	}
	else {
		arma::rowvec ejminusbis;
		ejminusbis = ejminus;
		ejminusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu;
			val = std::abs(val);
		});
		dmuejminus = arma::min(ejminusbis);
	}
	double dmuejplus = 0;
	if (ejplus(0) > ejplus(1)) {
		dmuejplus = inf;
	}
	else {
		arma::rowvec ejplusbis;
		ejplusbis = ejplus;
		ejplusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu;
			val = std::abs(val);
		});
		dmuejplus = arma::min(ejplusbis);
	}

	arma::rowvec ejequalbis;
	ejequalbis = ejequal;
	ejequalbis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; 
		val = std::abs(val);
	});
	double dmuejequal = arma::min(ejequalbis);


	//conv_to<fmat>::from(A)
	arma::rowvec ejp1bis(arma::conv_to<rowvec>::from(ejp1));
	/* 
	arma::rowvec ejp1bis(ejp1.n_elem);
	for (int in = 0; in < ejp1.n_elem; in++) {
		ejp1bis(in) = ejp1(in);
	}
	*/
	ejp1bis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; 
		val = std::abs(val);
	});
	double dmuejp1 = arma::min(ejp1bis);

	arma::rowvec all_dmu;
	all_dmu << dmuejminus << dmuejequal << dmuejplus;
	int pejp1_yjejzj1 = 0;
	if ((dmuejp1 == arma::min(all_dmu)) && ((compare_vec(ejp1, ejminus) || compare_vec(ejp1, ejequal) || compare_vec(ejp1, ejplus)))) {
		pejp1_yjejzj1 = 1;
	}
	else {
		pejp1_yjejzj1 = 0;
	}
	proba = p * pejp1_yjejzj1 + (1 - p) * pejp1_yjejzj0;
	return(proba);
}


double Bos::pejp1zj1_yjej(arma::urowvec ejp1, unsigned int yj, arma::urowvec ej, int mu, double p) {
	double proba = 0;
	arma::rowvec ejminus;
	ejminus << ej(0) << yj - 1;
	arma::rowvec ejequal;
	ejequal << yj << yj;
	arma::rowvec ejplus;
	ejplus << yj + 1 << ej(1);

	double dmuejminus = 0;
	double dmuejplus = 0;
	if (ejminus(0) > ejminus(1)) {
		dmuejminus = inf;
	}
	else {
		arma::rowvec ejminusbis;
		ejminusbis = ejminus;
		ejminusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu; 
			val = std::abs(val);
		});
		dmuejminus = arma::min(ejminusbis);
	}
	if (ejplus(0) > ejplus(1)) {
		dmuejplus = inf;
	}
	else {
		arma::rowvec ejplusbis;
		ejplusbis = ejplus;
		ejplusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu; 
			val = std::abs(val);
		});
		dmuejplus = arma::min(ejplusbis);
	}

	arma::rowvec ejequalbis;
	ejequalbis = ejequal;
	ejequalbis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; 
		val = std::abs(val);
	});
	double dmuejequal = arma::min(ejequalbis);


	/*
	arma::rowvec ejp1bis(ejp1.n_elem);
	for (int in = 0; in < ejp1.n_elem; in++) {
		ejp1bis(in) = (int)ejp1(in);
	}
	 */
	arma::rowvec ejp1bis(arma::conv_to<rowvec>::from(ejp1));

	ejp1bis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; 
		val = std::abs(val);
	});
	double dmuejp1 = arma::min(ejp1bis);

	arma::rowvec all_dmu = {dmuejminus,dmuejequal,dmuejplus};
	//all_dmu << dmuejminus << dmuejequal << dmuejplus;
	int pejp1_yjejzj1 = 0;
	if (dmuejp1 == arma::min(all_dmu) && (compare_vec(ejp1, ejminus) || compare_vec(ejp1, ejequal) || compare_vec(ejp1, ejplus))) {
		pejp1_yjejzj1 = 1;
	}
	else {
		pejp1_yjejzj1 = 0;
	}
	proba = p * pejp1_yjejzj1;
	return(proba);
}

double Bos::pejp1zj1_ej(arma::urowvec ejp1, arma::urowvec ej, int mu, double p) {
	double proba = 0;
	arma::uvec ej_bounds = linspace<uvec>(ej(0), ej(1), ej(1) - ej(0) + 1);
	for (unsigned int i = 0; i < ej_bounds.size(); ++i) {
		int yj = ej_bounds(i);
		proba += pejp1zj1_yjej(ejp1, yj, ej, mu, p);
	}
	proba = (double)proba / (ej(1) - ej(0) + 1);
	return(proba);
}


double Bos::pyj_ej(unsigned int yj, arma::urowvec ej) {
	double proba = 0;
	if (ej(0) <= yj && yj <= ej(1)) {
		proba = 1.0 / (ej(1) - ej(0) + 1);
	}
	else {
		proba = 0;
	}
	return(proba);
}

double Bos::pejp1_ej(arma::urowvec ejp1, arma::urowvec ej, int mu, double p) {
	double proba = 0;
	arma::uvec allyj;
	// we suppose ejp1 included in ej (not checked here)
	if (ejp1(1) == ejp1(0)) { // |ejp1|=1
		if (ejp1(1) < ej(1) && ejp1(0) > ej(0)) { //ejp1 doesn't touch any bound
			allyj << ejp1(0);
		}
		else {
			if (ejp1(1) < ej(1)) { // ejp1 doesn't touch right bound
				allyj << ejp1(0) << ejp1(0) + 1;
			}
			else { // ejp1 doesn't touch left bound
				allyj << ejp1(0) - 1 << ejp1(0);
			}
		}
	}
	else { // |ejp1|>1
		if (ejp1(1) < ej(1)) {// ejp1 doesn't touch right bound
			allyj << ejp1(1) + 1;
		}
		else { // ejp1 doesn't touch left bound
			allyj << ejp1(0) - 1;
		}
	}
	for (unsigned int i = 0; i < allyj.size(); ++i) {
		unsigned int yj = allyj(i);
		proba += pejp1_yjej(ejp1, yj, ej, mu, p) * pyj_ej(yj, ej);
	}
	return(proba);
}

double Bos::pej(arma::urowvec& ej,
	int j, int m, int mu, double p,
	arma::colvec& z1tozjm1)
{

	if (j == 1) {
		return(1);
	}

	if (ej.size() == 1) {
		unsigned int ejn = ej(0);
		ej = arma::urowvec(2);
		ej = ej.ones() * ejn;
	}

	arma::colvec z1tozjm2 = arma::vectorise(z1tozjm1, 0);
	double zjm1 = z1tozjm2(z1tozjm2.n_rows - 1);
	z1tozjm2.shed_row(z1tozjm2.n_rows - 1);

	double proba = 0;
	if (zjm1) { // zjm1 is known
		arma::umat tabint = allej(j - 1, m); // may have a problem here
		int nbtabint = tabint.n_rows;
		for (int i = 0; i < nbtabint; ++i) {
			arma::urowvec ejm1 = tabint.row(i);
			if ((ej(0) >= ejm1(0)) && (ej(1) <= ejm1(1))) { //to accelerate, check if ejm is included in ej
				proba += pejp1zj1_ej(ej, ejm1, mu, p) * pej(ejm1, j - 1, m, mu, p, z1tozjm2);
			}
		}
	}

	else { // zjm1 is unknown
		arma::umat tabint = allej(j - 1, m);
		int nbtabint = tabint.n_rows;
		for (int i = 0; i < nbtabint; ++i) {
			arma::urowvec ejm1 = tabint.row(i);
			if ((ej(0) >= ejm1(0)) && (ej(1) <= ejm1(1))) { //to accelerate, check if ejm is included in ej
				proba += pejp1_ej(ej, ejm1, mu, p) * pej(ejm1, j - 1, m, mu, p, z1tozjm2);
			}
		}
	}
	return(proba);
}


int Bos::getModeFromVec(uvec vector) {
	//int vecSize = vector.n_elem;
	unsigned int mode = 0;
	int count = 0;
	for (unsigned int mu = 0; mu < this->_m; mu++) {
		uvec ucount = find(vector == (mu + 1));
		int newcount = ucount.n_elem; // see if when there is none, it puts 0
		if (newcount >= count) {
			count = newcount;
			mode = mu+1;
		}
	}
	return mode;
}

void Bos::printMus() {
	// cout << " mus :" << endl;
	_mus.print();
}
void Bos::printPis() {
	// cout << " pis :" << endl;
	_pis.print();
}

