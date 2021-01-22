#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector cga_generate_chromosome(NumericVector prob_vec) {
   int len_prob_vec = prob_vec.length();
   NumericVector vect(len_prob_vec);
   NumericVector unifs = runif(len_prob_vec);
   int i;
   for (i=0;i<len_prob_vec;i++){
		if ( unifs[i] < prob_vec[i] ) {
			vect[i] = 1;
		}
	}
	return(vect);
}


// [[Rcpp::export]]
NumericVector cga (int chsize, int popsize, Function evalFunc){
    NumericVector prob_vec = rep(0.5, chsize);
    NumericVector chromosome1, chromosome2, winner, loser;
    NumericVector cost1, cost2;
    int i, t;
	while (1){
		chromosome1 = cga_generate_chromosome(prob_vec);
		chromosome2 = cga_generate_chromosome(prob_vec);
		cost1 = evalFunc(chromosome1);
		cost2 = evalFunc(chromosome2);
		winner = chromosome1;
		loser = chromosome2;
		if( cost2[0] < cost1[0] ) {
			winner = chromosome2;
			loser = chromosome1;
		}
	
		for (i=0;i < chsize;i++){	
			if (winner[i] != loser[i]){
				if(winner[i]==1){
					prob_vec[i] = prob_vec[i] + (1.0/popsize);
				}else{
					prob_vec[i] = prob_vec[i] - (1.0/popsize);
				}
			}
		}

		t = 0;
		for (i=0; i<chsize; i++){
			if(prob_vec[i] <= 0.001 || prob_vec[i] >= 0.999) {
				t++;
			}		
		}
		if (t == chsize) {
			break;
		}
	}
	return(prob_vec);
}	
