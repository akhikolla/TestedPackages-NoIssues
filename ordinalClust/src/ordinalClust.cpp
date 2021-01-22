// mixedData.cpp : définit le point d'entrée pour l'application console.
//


#include "Distribution.h"
#include "CoClusteringContext.h"
#include "ClusteringContext.h"
#include "ClassificationContext.h"
#include "ClassificationMContext.h"

#include "LogProbs.h"
#include "TabProbsResults.h"
#include <iostream>
#include <vector>
//#include "Rcpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include "BosPredict.h"



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
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
void testSeed(int seed){
	boost::mt19937 generator(seed);
	vector<double> vec(5);
	double prob = (double)1 / 5;
	std::fill(vec.begin(), vec.end(), prob);
	boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());
    //int temp = distribution(generator);
    //cout << " test: " << temp << endl;

	//int temp2 = distribution(generator);
    //cout << " test2: " << temp2 << endl;
    return;
}



//[[Rcpp::plugins(cpp11)]]


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
S4 coclust(NumericMatrix xMat, std::vector<unsigned int> myList,
	int kr, std::vector<int> kc, std::string init, int nbSEM, 
	int nbSEMburn, int nbRepeat, int nbindmini, const std::vector<int> m,
	std::vector<double> percentRandomB, int seed)
{

	
	Progress p(nbSEM, true);

	
	arma::mat x(xMat.begin(), xMat.nrow(), xMat.ncol(), false ) ;


	int dlistSize = myList.size();
	vector<urowvec> dlist(dlistSize);
	
	

	arma::urowvec tmp;
	for(int i=0; i<dlistSize;i++){
		if(i==(dlistSize-1)){
			tmp = linspace<urowvec>(myList[i], x.n_cols-1, x.n_cols - myList[i]);
			dlist[i] = tmp;
		}
		else{
			tmp = linspace<urowvec>(myList[i], myList[i+1]-1, myList[i+1] - myList[i]);
			dlist[i] = tmp;
		}

	}

	CoClusteringContext context(x, dlist, kr, kc, init, nbSEM, nbSEMburn, nbindmini, m, percentRandomB, seed);

	context.missingValuesInit();

		
	bool verif = false;
		
	int restart_init = 0;
	while (verif == false && restart_init < 15) {
		context.initialization();
		verif = context.verif();
		if (!verif) {
			restart_init++;
		}
	}
	if (!verif) {
		S4 t("ResultClustOrdinal");
		return t;
	}
	
	
	context.imputeMissingData();

	context.fillParameters(0);
	context.fillLabels(0);
	
	for (int iter = 0; iter < nbSEM; iter++) {
		p.increment();
		verif = false;

		
		for(int repeat=0; repeat<nbRepeat; repeat++){

				context.SEstepRow();
				context.sampleV();
				verif = context.verif();
				//int restart = 0;

				if(!verif){
					if( (init != "randomBurnin") || (iter > nbSEMburn)){
						S4 t("ResultCoclustOrdinal");
						return t;
					}
					else{
						vector<vector<int>> res = context.verification();
						context.noRowDegenerancy(res, iter);
						//context.MstepVW();
						verif = context.verif();
					}
				}
				else{
					context.MstepVW();
				}
			
		}

		if (!verif) {
			S4 t("ResultCoclustOrdinal");
			return t;
		}
		else {
			context.MstepVW();
			context.imputeMissingData();
		}

		for(int repeat=0; repeat<nbRepeat; repeat++){
			
				context.SEstepCol();
				context.sampleW();
				verif = context.verif();
				//int restart = 0;

				if(!verif){
					if( (init != "randomBurnin") || (iter > nbSEMburn)){
						S4 t("ResultCoclustOrdinal");
						return t;
					}
					else{
						vector<vector<int>> res = context.verification();
						context.noColDegenerancy(res, iter);
						//context.MstepVW();
						verif = context.verif();
					}
				}
				else{
					context.MstepVW();
				}
			
		}
		
		if (!verif) {
			S4 t("ResultCoclustOrdinal");
			return t;
		}
		else {
			context.MstepVW();
			context.imputeMissingData();
		}
		if(iter>0){
			context.fillParameters(iter);
			context.fillLabels(iter);
		}
			
	}
	
	
	context.getBurnedParameters();
	
	
	context.SEstepRow();
	context.SEstepCol();  // TODO : remettre // added so that probaV is actualized after burned params
	context.sampleVWStock(); // TODO : remettre
	context.computeICL();
	
    return context.returnCoclustering();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
S4 clust(NumericMatrix xMat, std::vector<unsigned int> myList, 
	int kr, std::string init, int nbSEM, int nbSEMburn, int nbindmini, 
	const std::vector<int> m, std::vector<double> percentRandomB, int seed)
{
	Progress p(nbSEM, true);

	arma::mat x(xMat.begin(), xMat.nrow(), xMat.ncol(), false ) ;


	int dlistSize = myList.size();
	vector<urowvec> dlist(dlistSize);

	arma::urowvec tmp;
	for(int i=0; i<dlistSize;i++){
		if(i==(dlistSize-1)){
			tmp = linspace<urowvec>(myList[i], x.n_cols-1, x.n_cols - myList[i]);
			dlist[i] = tmp;
		}
		else{
			tmp = linspace<urowvec>(myList[i], myList[i+1]-1, myList[i+1] - myList[i]);
			dlist[i] = tmp;
		}

	}


	ClusteringContext context(x, dlist, kr, init, nbSEM, nbSEMburn, nbindmini, m, percentRandomB, seed);

	context.missingValuesInit();


	bool verif = false;
		
	int restart_init = 0;
	while (verif == false && restart_init < 15) {
		context.initialization();
		verif = context.verif();
		if (!verif) {
			restart_init++;
		}
	}
	if (!verif) {
		S4 t("ResultClustOrdinal");
		return t;
	}
	

	context.fillParameters(0);
	context.fillLabels(0);
	
	for (int iter = 0; iter < nbSEM; iter++) {
		p.increment();
		verif = false;
		//int restart = 0;
		context.SEstepRow();
		context.sampleVW();
		verif = context.verif();
		/* while (verif == false && restart < 25) {
			verif = context.verif();
			if (!verif) {
				context.sampleVW();
				restart++;
			}
		} */


		
		if(!verif){
			if( (init != "randomBurnin") || (iter > nbSEMburn)){
				S4 t("ResultClustOrdinal");
				return t;
			}
			else{
				vector<vector<int>> res = context.verification();
				context.noRowDegenerancy(res, iter);
				//context.MstepVW();
				verif = context.verif();
			}
		}
		else{
			context.MstepVW();
		}

		if (!verif) {
			S4 t("ResultClustOrdinal");
			return t;
		}
		else {
			context.MstepVW();
			context.imputeMissingData();	
		}
		if(iter>0){
			context.fillParameters(iter);
			context.fillLabels(iter);
		}	
	}

	context.getBurnedParameters();
	
	// missing: take the median or mode for the parameters rho gamma
	/*for (int iter = 0; iter < nbSEM; iter++) {
		context.sampleVW();
	}*/
	context.SEstepRow();
	context.sampleVWStock();

	//context.printResults();

	context.computeICL();
	
    return context.returnClustering();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
S4 classif(NumericMatrix xMat, NumericVector yVec, std::vector<unsigned int> myList, 
	int kr, std::vector<int> kc, std::string init, int nbSEM, int nbSEMburn, int nbindmini, 
	const std::vector<int> m, std::vector<double> percentRandomB, int seed)
{

	Progress p(nbSEM, true);

	arma::mat x(xMat.begin(), xMat.nrow(), xMat.ncol(), false ) ;
	arma::vec y(yVec.begin(), yVec.size(), false);

	int dlistSize = myList.size();
	vector<urowvec> dlist(dlistSize);

	arma::urowvec tmp;
	for(int i=0; i<dlistSize;i++){
		if(i==(dlistSize-1)){
			tmp = linspace<urowvec>(myList[i], x.n_cols-1, x.n_cols - myList[i]);
			dlist[i] = tmp;
		}
		else{
			tmp = linspace<urowvec>(myList[i], myList[i+1]-1, myList[i+1] - myList[i]);
			dlist[i] = tmp;
		}

	}

	ClassificationContext context(x, y, dlist, kr, kc, init, nbSEM, nbSEMburn, nbindmini, m, percentRandomB, seed);



	context.missingValuesInit();

	bool verif = false;

		
	int restart_init = 0;
	while (verif == false && restart_init < 15) {
	
		context.initialization();
	
		verif = context.verif();
		if (!verif) {
			restart_init++;
		}
	}

	
	if (!verif) {
		S4 t("ResultClassifOrdinal");
		return t;
	}
	
	context.imputeMissingData();

	context.fillParameters(0);
	context.fillLabels(0);

	for (int iter = 0; iter < nbSEM; iter++) {
		p.increment();
		verif = false;
		//int restart = 0;

		context.SEstep();
		context.sampleVW();

		verif = context.verif();

		/* while (verif == false && restart < 25) {
			verif = context.verif();
			if (!verif) {
				context.sampleVW();
				restart++;
			}
		} */

		if(!verif){
			int restart = 0;
			while(!verif && restart<10){
				//cout << "restart" << restart << endl;
				if( (init != "randomBurnin") || (iter > nbSEMburn)){
					S4 t("ResultClassifOrdinal");
					return t;
				}
				else{
					vector<vector<int>> res = context.verification();
					context.noColDegenerancy(res, iter);
					//context.MstepVW();
					verif = context.verif();
				}
				restart++;
			}
			
		}
		else{
			context.MstepVW();
		}
		if (!verif) {
			S4 t("ResultClassifOrdinal");
			return t;
		}
		else {
			context.MstepVW();
			context.imputeMissingData();	
			context.putParamsToZero();
		}

		if(iter>0){
			context.fillParameters(iter);
			context.fillLabels(iter);
		}	
	}

	
	context.getBurnedParameters();
	context.SEstep(); 
	context.sampleVWStock();

	//context.printResults();
	//context.predict();

	context.computeICL();

    return context.returnClassification();

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
S4 classifM(NumericMatrix xMat, NumericVector  yVec, std::vector<unsigned int> myList, 
	int kr,  std::string init, int nbSEM, 
	int nbSEMburn, int nbindmini, const std::vector<int> m, int seed)
{
		Progress p(nbSEM, true);
		
		p.increment();

		arma::mat x(xMat.begin(), xMat.nrow(), xMat.ncol(), false ) ;
		arma::vec y(yVec.begin(), yVec.size(), false);

		int dlistSize = myList.size();
		vector<urowvec> dlist(dlistSize);

		arma::urowvec tmp;
		for(int i=0; i<dlistSize;i++){
			if(i==(dlistSize-1)){
				tmp = linspace<urowvec>(myList[i], x.n_cols-1, x.n_cols - myList[i]);
				dlist[i] = tmp;
			}
			else{
				tmp = linspace<urowvec>(myList[i], myList[i+1]-1, myList[i+1] - myList[i]);
				dlist[i] = tmp;
			}

		}


		ClassificationMContext context(x, y, dlist, kr, init, nbSEM, nbSEMburn, nbindmini, m, seed);

		p.increment();

		context.missingValuesInit();

		bool verif = false;

			
		int restart_init = 0;
		while (verif == false && restart_init < 15) {
		
			context.initialization();
		
			verif = context.verif();
			if (!verif) {
				restart_init++;
				//cout << "empty blocks--restarting number " << restart_init << endl;
			}
		}
		if (!verif) {
			//cout << "degenerancy of algorithm" << endl;
			S4 t("ResultClassifOrdinal");
			return t;
		}
		
		context.imputeMissingData();


		p.increment();

		context.MstepVW();
		context.imputeMissingData();

		//context.printResults();
		//context.predict();

		context.computeICL();

		p.increment();

	    return context.returnClassification();
}


double logsum(rowvec logx) {
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
double logsum(rowvec logx);


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
S4 prediction(S4 classif, NumericMatrix xMat_topredict, int seed)
{

	arma::mat x_topredict(xMat_topredict.begin(), xMat_topredict.nrow(), xMat_topredict.ncol(), false ) ;

	int N_topredict = 0;
	
	N_topredict = x_topredict.n_rows;
	
	int kr = classif.slot("kr");

	vec kc = classif.slot("kc");
	std::string name = classif.slot("name");
	
	// Computing the log-probabilites

	mat logprobaV_topredict(N_topredict, kr);
	logprobaV_topredict.zeros();
	for (int i = 0; i < N_topredict; i++) {
		vec gammas = classif.slot("pi");
		logprobaV_topredict.row(i) = conv_to<rowvec>::from(log(gammas));
	}

	List dlist = classif.slot("dlist");
	int number_distrib = classif.slot("number_distrib");
	vec J = classif.slot("J");

	List W = classif.slot("W");

	List params = classif.slot("params");

	mat V = classif.slot("V");

	int nbSEM = classif.slot("nbSEM");

	vec m = classif.slot("m");

	
	int iterm = 0;


	for (int idistrib = 0; idistrib < number_distrib; idistrib++)
	{

		//TabProbsResults result(0, 0);
		
		List paramsIdistrib = params[idistrib];
		
		mat tmpW = W[idistrib];

								
		//if(name == "Bos") {
			uvec indexes = dlist[idistrib];
			mat tmp_x = x_topredict.cols(indexes);
			mat pis = paramsIdistrib["pis"];
			mat mus = paramsIdistrib["mus"];
			BosPredict object(kr, kc(idistrib), m(iterm), pis, mus, seed);
			// for now we just initialize at random TODO: impute Missing data in a better way
			tmp_x = object.missingValuesInit(tmp_x);
			iterm ++;
			//tmp_x.print();
			logprobaV_topredict += object.SEstep_predict(tmpW, tmp_x);
		//}
		

	} // end for logprobas
	
	mat probaV_topredict(N_topredict, kr);
	probaV_topredict.zeros();

	mat V_topredict(N_topredict, kr);
	V_topredict.zeros();

	// Computing the probabilites
	for (int i = 0; i < N_topredict; i++) {
		for (int k = 0; k < kr; k++) {
			probaV_topredict(i, k) = exp(logprobaV_topredict(i, k) - logsum(logprobaV_topredict.row(i)));
		}
	}

	
	// sampling _V_topredict
	mat countV_topredict = zeros(N_topredict, kr);

	for (int iter = 0; iter < nbSEM; iter++) {
		V_topredict.zeros();
		//std::default_random_engine gen;
		for (int i = 0; i < N_topredict; i++) {
			
			//RANDOM
			/* rowvec vec = probaV_topredict.row(i);			
			discrete_distribution<> dis(vec.begin(), vec.end());			
			random_device _rd;
			mt19937 gen(_rd());
			int sample = dis(gen); */
			rowvec vec = probaV_topredict.row(i);
			boost::mt19937 generator(seed);
			boost::random::discrete_distribution<int> distribution (vec.begin(),vec.end());
			int sample = distribution(generator);


			V_topredict(i, sample) = 1;
			countV_topredict(i, sample) += 1;
		}
	}

	//determinging final partitions
	V_topredict.zeros();
	for (int i = 0; i < N_topredict; i++) {
		int maxind = countV_topredict.row(i).index_max();
		V_topredict(i, maxind) = 1;
	}
	vector<int> zr_topredict(N_topredict);
	// labels
	for (int i = 0; i < N_topredict; i++) {
		uvec col1 = find(V_topredict.row(i) == 1);
		zr_topredict[i] = col1(0) + 1;
	}
	//context.predict();
    //return context.returnClassification();


    S4 x("ResultPredictionOrdinal");

    x.slot("zr_topredict") = zr_topredict;
    x.slot("V_topredict") = V_topredict;

    return x;
}


int unsigned_to_signed(unsigned x)
{
	if (x <= INT_MAX)
		return static_cast<int>(x);

	if (x >= INT_MIN)
		return static_cast<int>(x - INT_MIN) + INT_MIN;

	throw x; // Or whatever you like
}


bool compare_vec(arma::urowvec vec1, arma::rowvec vec2) {
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix allej(int j, int m) { // TODO . to implement
	arma::umat result;
	if (j == 1) {
		//result << 1 << m << arma::endr;
		result = {{(unsigned int)1,(unsigned int)m}};
		return(wrap(result));
	}
	if (j > m) {
		return(wrap(result));
	}
	arma::rowvec indicesi = arma::linspace<arma::rowvec>(1, m - j + 1, m - j + 1);
	for (unsigned int i = 0; i < indicesi.size(); ++i) {
		int sizeej = indicesi(i);
		arma::rowvec indicesj = arma::linspace<arma::rowvec>(1, m - sizeej + 1, m - sizeej + 1);
		for (unsigned int j = 0; j < indicesj.size(); ++j) {
			int binf = indicesj(j);
			int bsup = binf + sizeej - 1;
			arma::umat rowvector = {{(unsigned int)binf, (unsigned int)bsup}};
			//arma::umat rowvector;
			//rowvector << binf << bsup << arma::endr;
			result = arma::join_vert(result, rowvector);
		}
	}

	return(wrap(result));
}


double pejp1_yjej(arma::urowvec ejp1, int yj, arma::urowvec ej, int mu, double p) {

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

	arma::rowvec ejp1bis(ejp1.n_elem);
	for (int in = 0; in < ejp1.n_elem; in++) {
		ejp1bis(in) = ejp1(in);
	}
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


double pejp1zj1_yjej(arma::urowvec ejp1, unsigned int yj, arma::urowvec ej, int mu, double p) {
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

	/*arma::urowvec ejp1bis;
	ejp1bis = ejp1;
	ejp1bis.for_each([mu](arma::urowvec::elem_type& val) {
		val -= mu; 
		val = val;
	});*/
	arma::rowvec ejp1bis(ejp1.n_elem);
	for (int in = 0; in < ejp1.n_elem; in++) {
		ejp1bis(in) = (int)ejp1(in);
	}
	ejp1bis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; 
		val = std::abs(val);
	});
	double dmuejp1 = arma::min(ejp1bis);

	arma::rowvec all_dmu;
	all_dmu << dmuejminus << dmuejequal << dmuejplus;
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double pejp1zj1_ej(NumericVector ejp1Vec , NumericVector ejVec, int mu, double p) {

	arma::vec ejp1R(ejp1Vec.begin(), ejp1Vec.size(), false);
	urowvec ejp1 = conv_to<urowvec>::from(ejp1R);

	arma::vec ejR(ejVec.begin(), ejVec.size(), false);
	urowvec ej = conv_to<urowvec>::from(ejR);

	double proba = 0;
	arma::uvec ej_bounds = linspace<uvec>(ej(0), ej(1), ej(1) - ej(0) + 1);
	for (unsigned int i = 0; i < ej_bounds.size(); ++i) {
		int yj = ej_bounds(i);
		proba += pejp1zj1_yjej(ejp1, yj, ej, mu, p);
	}
	proba = (double)proba / (ej(1) - ej(0) + 1);
	return(proba);
}


double pyj_ej(unsigned int yj, arma::urowvec ej) {
	double proba = 0;
	if (ej(0) <= yj && yj <= ej(1)) {
		proba = 1.0 / (ej(1) - ej(0) + 1);
	}
	else {
		proba = 0;
	}
	return(proba);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double pejp1_ej(NumericVector ejp1Vec, NumericVector ejVec, int mu, double p) {


	arma::vec ejp1R(ejp1Vec.begin(), ejp1Vec.size(), false);
	urowvec ejp1 = conv_to<urowvec>::from(ejp1R);

	arma::vec ejR(ejVec.begin(), ejVec.size(), false);
	urowvec ej = conv_to<urowvec>::from(ejR);

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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double pej(NumericVector ejVec,
	int j, int m, int mu, double p,
	NumericVector z1tozjm1Vec)
{


	arma::vec ejR(ejVec.begin(), ejVec.size(), false);
	urowvec ej = conv_to<urowvec>::from(ejR);

	arma::vec z1tozjm1R(z1tozjm1Vec.begin(), z1tozjm1Vec.size(), false);
	colvec z1tozjm1 = conv_to<colvec>::from(z1tozjm1R);

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
		NumericMatrix tabintMat = allej(j - 1, m); // may have a problem here

		arma::mat tabintR(tabintMat.begin(), tabintMat.nrow(), tabintMat.ncol(), false ) ;
		arma::umat tabint = conv_to<umat>::from(tabintR);

		int nbtabint = tabint.n_rows;
		for (int i = 0; i < nbtabint; ++i) {
			arma::urowvec ejm1 = tabint.row(i);
			if ((ej(0) >= ejm1(0)) && (ej(1) <= ejm1(1))) { //to accelerate, check if ejm is included in ej
				//f(NumericVector(a.begin(),a.end()))
				proba += pejp1zj1_ej(NumericVector(ej.begin(),ej.end()), 
					NumericVector(ejm1.begin(),ejm1.end()), mu, p) * 
				pej(NumericVector(ejm1.begin(),ejm1.end()), j - 1, m, mu, p, 
					NumericVector(z1tozjm2.begin(),z1tozjm2.end()));
			}
		}
	}

	else { // zjm1 is unknown
		NumericMatrix tabintMat = allej(j - 1, m); // may have a problem here

		arma::mat tabintR(tabintMat.begin(), tabintMat.nrow(), tabintMat.ncol(), false ) ;
		arma::umat tabint = conv_to<umat>::from(tabintR);
		int nbtabint = tabint.n_rows;
		for (int i = 0; i < nbtabint; ++i) {
			arma::urowvec ejm1 = tabint.row(i);
			if ((ej(0) >= ejm1(0)) && (ej(1) <= ejm1(1))) { //to accelerate, check if ejm is included in ej
				proba += pejp1_ej(NumericVector(ej.begin(),ej.end()), NumericVector(ejm1.begin(),ejm1.end()), mu, p) * 
				pej(NumericVector(ejm1.begin(),ejm1.end()), j - 1, m, mu, p, 
					NumericVector(z1tozjm2.begin(),z1tozjm2.end()));
			}
		}
	}
	return(proba);
}





