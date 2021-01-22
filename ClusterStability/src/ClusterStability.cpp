#include <Rcpp.h>
#include <algorithm>
#include <map>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace Rcpp;

inline bool isNan(double x) {
    return x != x;
}
 
std::vector<std::vector<int> > combinations;
 std::vector<int> combination;
   
   // [[Rcpp::export]]   
 IntegerVector Reorder(NumericVector data) { 	
	 std::map<int, int> d;	
	std::vector<int> results= Rcpp::as<std::vector<int> > (data);	
	int startgroup=1; 
	for (int j=0; j<(int)results.size();j++) {
		if(d.count(results[j])==0) {
			d[results[j]]=startgroup;
			startgroup++;
		}
	}
	 for (int j=0; j<(int)results.size(); j++) {    
		results[j]=d[results[j]];
	 }	 
	return Rcpp::wrap(results);
 }
   
 List k6combination(NumericVector data, int selector=0) {
 	combination.clear();
 	combinations.clear();
 	for (int i=0; i< data.size();i++) {
 		for (int j=i+1;j<data.size();j++) {
 			for (int k=j+1;k<data.size();k++) {
 				for (int l=k+1;l<data.size();l++) {
 					for (int m=l+1;m<data.size();m++) {
 						for (int n=m+1;n<data.size();n++) {
							if (data[i]==selector||data[j]==selector||data[k]==selector||data[l]==selector||data[m]==selector||data[n]==selector||selector==0) {
								combination.clear();
								combination.push_back(data[i]);
								combination.push_back(data[j]);
								combination.push_back(data[k]);
								combination.push_back(data[l]);
								combination.push_back(data[m]);
								combination.push_back(data[n]);
								combinations.push_back(combination);
							}
						}
					}			
				}
 			}  
 		}
 	}
 	return Rcpp::wrap( combinations );
 }
 
 List k5combination(NumericVector data, int selector=0) {
 	combination.clear();
 	combinations.clear();
 	for (int i=0; i< data.size();i++) {
 		for (int j=i+1;j<data.size();j++) {
 			for (int k=j+1;k<data.size();k++) {
 				for (int l=k+1;l<data.size();l++) {
 					for (int m=l+1;m<data.size();m++) {
						if (data[i]==selector||data[j]==selector||data[k]==selector||data[l]==selector||data[m]==selector||selector==0) {
							combination.clear();
							combination.push_back(data[i]);
							combination.push_back(data[j]);
							combination.push_back(data[k]);
							combination.push_back(data[l]);
							combination.push_back(data[m]);
							combinations.push_back(combination);
						}
					}			
				}
 			}  
 		}
 	}
 	return Rcpp::wrap( combinations );
 }
 
 List k4combination(NumericVector data, int selector=0) {
 	combination.clear();
 	combinations.clear();
 	for (int i=0; i< data.size();i++) {
 		for (int j=i+1;j<data.size();j++) {
 			for (int k=j+1;k<data.size();k++) {
 				for (int l=k+1;l<data.size();l++) {
 					if (data[i]==selector||data[j]==selector||data[k]==selector||data[l]==selector||selector==0) {
						combination.clear();
						combination.push_back(data[i]);
						combination.push_back(data[j]);
						combination.push_back(data[k]);
						combination.push_back(data[l]);
						combinations.push_back(combination);
					}
				}
 			}  
 		}
 	}
 	return Rcpp::wrap( combinations );
 }
 
 
 List k3combination(NumericVector data, int selector=0) {
 	combination.clear();
 	combinations.clear();
 	for (int i=0; i< data.size();i++) {
 		for (int j=i+1;j<data.size();j++) {
 			for (int k=j+1;k<data.size();k++) {
 				if (data[i]==selector||data[j]==selector||data[k]==selector||selector==0) {
 					combination.clear();
 					combination.push_back(data[i]);
 					combination.push_back(data[j]);
 					combination.push_back(data[k]);
 					combinations.push_back(combination);
 				}
 			}  
 		}
 	}
 	return Rcpp::wrap( combinations );
 }
 
 List k2combination(NumericVector data, int selector=0) {	
 	combination.clear();
 	combinations.clear();
 	for (int i=0; i< data.size();i++) {
 		for (int j=i+1;j<data.size();j++) {
 			if (data[i]==selector||data[j]==selector||selector==0) {
 				combination.clear();
 				combination.push_back(data[i]);
 				combination.push_back(data[j]);
 			 	combinations.push_back(combination); 
 			}
 		}
 	}
 	return Rcpp::wrap( combinations );
 }
 
 
  // [[Rcpp::export]]
 List Kcombination(NumericVector data, int k, int selector=0) {
 	combinations.clear();
 	combination.clear();
 	if (k==2) return k2combination(data, selector);
 	if (k==3) return k3combination(data, selector);
 	if (k==4) return k4combination(data, selector);
 	if (k==5) return k5combination(data, selector);
 	if (k==6) return k6combination(data, selector);
 	if (k>6) {
 		Rcpp::Rcout <<"Current function is limited to k between 2-6"<< std::endl;
 	}
 	return Rcpp::wrap( combinations );
 }

 
 // [[Rcpp::export]]
 NumericVector calculate_indices(Rcpp::List r_combinations,Rcpp::List r_partition, NumericVector r_indice, double total_indice) {	 
	 std::vector<double>  total_indices;
	 std::vector<std::vector<int> > partitions;
	 combinations.clear();		  
		  for (int i=0; i<r_combinations.size();i++) {		 			
			 std::vector<int> tmp = Rcpp::as< std::vector<int> >(r_combinations[i]);		 
			 combinations.push_back(tmp);		
		  }	  
	   for (int i=0; i<r_partition.size();i++) {	  		  
		  std::vector<int> tmp = Rcpp::as<std::vector<int> >(r_partition[i]);		 
		  partitions.push_back(tmp);	  
	   }	 
	   std::vector<double> indice= Rcpp::as<std::vector<double> >(r_indice);	  	 	 
	 //int len_combination=combinations[0].size();	 
	  for (int i=0; i<(int)combinations.size();i++) {		
	   total_indices.push_back(0);
		for (int j=0; j<(int)partitions.size(); j++) {
			if (partitions[j][combinations[i][0]-1]==partitions[j][combinations[i][1]-1]) 	total_indices[i]+=indice[j];					
		}		
		total_indices[i]/=total_indice;		
	  }
	return Rcpp::wrap( total_indices );
}


 // [[Rcpp::export]]
 NumericVector calculate_individual_PSG(double k, Rcpp::List r_combinations, NumericVector total_indices, NumericVector indices) {
	//--This is the new function (old version Janvier 2015)	
	std::vector<double>  PSG;
	double len=(double)indices.size();
	for (int i=0;i<len;i++) PSG.push_back(0.0);	
	double kl1=(k/(k-1.0));
	double kl2=((k-1.0)/k);
	for (int i=0; i<r_combinations.size();i++) {		 					
		std::vector<int> combination = Rcpp::as< std::vector<int> >(r_combinations[i]);		 
		int len_combination=combination.size();
		//The max is done here instead of in _calculate_indices()	
		total_indices[i]=std::max(kl1*(total_indices[i]-(1.0/k)), k*( 1.0-total_indices[i]-kl2));
		
		for (int j=0; j<len_combination;j++) {												
			PSG[combination[j]-1]+=total_indices[i];
		}			
	}	  
	
		 for (int i=0;i<len;i++) {				
			 PSG[i]=((PSG[i])/(len-1.0));
		 }	
	return Rcpp::wrap(PSG);
}
 
 // [[Rcpp::export]]
 NumericVector calculate_individual_PSG_exact(double k, Rcpp::List r_combinations, NumericVector total_singletons,NumericVector total_indices, NumericVector indices, double pnk, double ptildenk) {
	//--This is the new function (Mars 2015 following Matthieu Willems discussion)
	//--Note: for speed we expect the precalculated pnk p(n,k) and ptildenk(n,k)
	std::vector<double>  PSG;
	
	// Precalculate variables
	double n=(double)indices.size();
	double A=(1.0/(1.0-pnk));
	double B=(1.0/pnk);
	double C=(1.0/(1.0-ptildenk));
	double D=(1.0/ptildenk);
	
	//Handle limit cases
	if ((ptildenk-1.0)<0.000000001) C=1.0;
	if ((ptildenk-1.0)<0.000000001) D=0.0;
	if (std::isnan(C)) C=1.0;
	if (std::isnan(D)) D=0.0;
	
	double len=(double)indices.size();
	for (int i=0;i<len;i++) PSG.push_back(0.0);	
	
	for (int i=0; i<r_combinations.size();i++) {		 					
		std::vector<int> combination = Rcpp::as< std::vector<int> >(r_combinations[i]);		 
		int len_combination=combination.size();
		
		total_indices[i]=std::max(A*(total_indices[i]-pnk), B*(pnk-total_indices[i]));
		
		for (int j=0; j<len_combination;j++) {												
			PSG[combination[j]-1]+=total_indices[i];
		}			
	}	  
		
		 for (int i=0;i<len;i++) {				
			 PSG[i]=((PSG[i])/(len));			
		 }	

		for (int i=0;i<len;i++)  {
			 //PSG[i]= PSG[i]+(std::max(C*(total_singletons[i]-ptildenk),0.0)*(1.0/n));
			 // old 10 mars 2015  
			 PSG[i]= PSG[i]+(std::max(C*(total_singletons[i]-ptildenk),D*(ptildenk-total_singletons[i]))*(1.0/n));
		}
		
	return Rcpp::wrap(PSG);
}
 
  // [[Rcpp::export]]
 NumericVector calculate_individual_PSG_approximative(double k, Rcpp::List r_combinations, NumericVector total_singletons,NumericVector total_indices, NumericVector indices) {
	//--This is the new function (Mars 2015 following Matthieu Willems discussion)
	//printf("Calculating PSG");
	std::vector<double>  PSG;
	
	// Precalculate variables
	double n=(double)indices.size();
	double A=(k/(k-1.0));
	double B=(1.0/k);
	double C=pow(k,n-1.0)/(pow(k, n-1.0)-pow(k-1.0,n-1.0));
	double D=pow(k-1.0,n-1.0)/pow(k, n-1.0);
	double E=pow(k,n-1.0)/pow(k-1.0, n-1.0);
	
	if (std::isnan(C)) C=1.0;
	if (std::isnan(D)) D=0.0;
	if (std::isnan(D)) E=1.0;
	double len=(double)indices.size();
	for (int i=0;i<len;i++) PSG.push_back(0.0);	
	
	for (int i=0; i<r_combinations.size();i++) {		 					
		std::vector<int> combination = Rcpp::as< std::vector<int> >(r_combinations[i]);		 
		int len_combination=combination.size();
		//This is the first max		
		total_indices[i]=std::max(A*(total_indices[i]-B), k*(B-total_indices[i]));
		
		for (int j=0; j<len_combination;j++) {												
			PSG[combination[j]-1]+=total_indices[i];
		}			
	}	  
	
	 for (int i=0;i<len;i++) {				
			 PSG[i]=((PSG[i])/(len));			
	 }	
		
	for (int i=0;i<len;i++)  {
			 PSG[i]= PSG[i]+(std::max(C*(total_singletons[i]-D),E*(D-total_singletons[i]))*(1.0/n));
	}
	return Rcpp::wrap(PSG);
}
