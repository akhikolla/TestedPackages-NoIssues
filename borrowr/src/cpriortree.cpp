/*
 *  borrowr: estimate population average treatment effects with borrowing between data sources.
 *  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP getRunif() {
  RNGScope scope;
  NumericVector x = runif(1);
  return x;
}

RcppExport SEXP getRnorm() {
  RNGScope scope;
  NumericVector x = rnorm(1, 0, 1);
  return x;
}

// [[Rcpp::export]]
List priortree(NumericVector cut_lens, double alpha = 0.95, double beta = 2) {
  // at minimum, function needs to return:
  //  - nodes
  //  - var
  //  - cut
  //  - leaves

  List res;
  double p = cut_lens.size();
  std::vector<int> node(1); node[0] = 1;
  std::vector<double> leave(0);
  std::vector<int> var(0);
  std::vector<int> cut(0);
  // std::vector<int> depth(1);
  // std::vector<double> prob_terminal(1);
  //nodes[0] = 1;

  // terminate loop when I'm on the the last element, and it's terminal
  bool done = false;
  int i = 0;
  while(!done) {
    //int nd = node[i];
    double nd = node[i];
    // determine whether terminal or non-terminal
    double dp = std::floor(log2(nd));
    //depth.push_back(dp);`
    double pt = 1.0 - alpha * pow(1.0 + dp, - beta);
    NumericVector r1 = getRunif();
    //prob_terminal.push_back(pt);
    // if non-terminal, generate children.
    if (r1[0] > pt) { // if not a terminal node...
      int nn = nd * 2;
      leave.push_back(0.0);
      node.push_back(nn);
      node.push_back(nn + 1);
      // choose splitting variable:
      NumericVector r2 = getRunif();
      NumericVector r3 = getRunif();
      int wvar = std::floor(r2[0] * p);
      double cl = cut_lens[wvar];
      int wcut = std::floor(r3[0] * cl);
      var.push_back(wvar);
      cut.push_back(wcut);
    } else { // then it's a terminal node
      NumericVector lv = getRnorm();
      leave.push_back(lv[0]);
      var.push_back(0);
      cut.push_back(0);
    }
    // if terminal and no more children, exit
    //bool on_last = FALSE;
    //if (i == 15) {
      //done = true;
    //}
    // node.push_back(i);
    int s = node.size();
    if(i + 1 == s) {
      done = true;
    }
    i++;
  }
  // int p = node.size();
  res["node"] = node;
  res["leave"] = leave;
  res["var"] = var;
  res["cut"] = cut;
  //res["depth"] = depth;
  //res["prob_terminal"] = prob_terminal;
  //res["done"] = done;
  // res["p"] = p;
  return res;
}


/*** R
#set.seed(21)
#(out <- priortree(c(8001, 3, 2), beta = 2))

# library(microbenchmark)
# microbenchmark(
#   priortree(c(8001, 3, 2), beta = 1),
#   p2(c(8001, 3, 2), beta = 1),
#   times = 1e3)

*/




// List priortree(NumericVector cut_lens, NumericVector cum_sum_var_probs, double sigma_mu, int tiers = 6, double alpha = 0.95, double beta = 2) {
// 	List res;
//
//   int p = cut_lens.size();
//
//   //std::vector<double> nodes(pow(2, tiers)  * 2 - 1);
//   std::vector<int> nodes(pow(2, tiers)  * 2 - 1);
//   // NumericVector nodes(pow(2, tiers)  * 2 - 1);
//   // std::generate (nodes.begin(), nodes.end(), gen);
//
//   //nn    <- length(nodes)
//   int n_nodes = nodes.size();
//   //tier  <- floor(log2(nodes))
//   std::vector<int> tier(n_nodes);
//   //NumericVector tier(n_nodes);
//   NumericVector terminal_prob(n_nodes);
//   //NumericVector terminal_prob = 1 - alpha * pow(1 + tier, - beta);
//
//   NumericVector ran1 = runif(n_nodes);
//   NumericVector ran2 = runif(n_nodes);
//   NumericVector ran3 = runif(n_nodes);
//   NumericVector leaves = rnorm(n_nodes, 0, sigma_mu);
//
//   std::vector<int> var(n_nodes);
//   std::vector<int> cut(n_nodes);
//
//   LogicalVector terminal(n_nodes); // = ran1 < terminal_prob;
//   LogicalVector killed(n_nodes);
//   // LogicalVector branch(n_nodes);
//
//   for(int i = 0; i < n_nodes; ++i){
//     //NumericVector child_nodes(2);
//     nodes[i] = i + 1;
//     tier[i] = floor(log2(nodes[i]));
//     terminal_prob[i] = 1 - alpha * pow(1 + tier[i], - beta);
//     terminal[i] = ran1[i] < terminal_prob[i];
//
//     if(!terminal[i]) {
//       // var[i] = floor(ran2[i] * p);
//       for(int jj = 0; jj < p; ++jj) {
//         if(ran2[i] < cum_sum_var_probs[jj]) {
//           var[i] = jj;
//           break;
//         }
//       }
//       cut[i] = floor(ran3[i] * cut_lens[var[i]]);
//       leaves[i] = 0;
//     }
//     //printf("foo\n
//     // if terminal and not at max tiers, kill children nodes
//     //printf("tier = %i, tiers = %d, node = %i, i = %d; child 1 = %d, child 2 = %d\n",
//       //tier[i], tiers, nodes[i], i, 2 * i + 1, 2 * i + 2);
//     if((tier[i] < tiers) & (terminal[i] | killed[i])) {
//       //printf("if\n");
//       //printf("tier = %f, tiers = %d, i = %d; child 1 = %d, child 2 = %d\n", tier[i], tiers, i, 2 * i, 2 * i + 1);
//       killed[2 * i + 1] = true;
//       killed[2 * i + 2] = true;
//     }
//     //child_nodes[0] = 2 * nodes[i];
//     //child_nodes[1] = 2 * nodes[i] + 1;
//   }
//
//
//   // for debugging
// 	//res["seed"]     = seed;
// 	//res["p"]             = p;
// 	//res["cut_lens"]      = cut_lens;
// 	//res["sigma_mu"]      = sigma_mu;
// 	//res["tiers"]         = tiers;
// 	//res["alpha"]         = alpha;
// 	//res["beta"]          = beta;
// 	//res["nodes"]         = nodes;
// 	//res["n_nodes"]       = n_nodes;
// 	//res["tier"]          = tier;
// 	//res["terminal_prob"] = terminal_prob;
// 	//res["ran1"]          = ran1;
// 	//res["ran2"]          = ran2;
// 	//res["ran3"]          = ran3;
// 	//res["var"]           = var;
// 	//res["cut"]           = cut;
// 	//res["leaves"]        = leaves;
// 	//res["terminal"]      = terminal;
// 	//res["killed"]        = killed;
//
// 	res["p"]             = p;
//   res["cum_sum_var_probs"] = cum_sum_var_probs;
// 	//res["cut_lens"]      = cut_lens;
// 	//res["sigma_mu"]      = sigma_mu;
// 	//res["tiers"]         = tiers;
// 	//res["alpha"]         = alpha;
// 	//res["beta"]          = beta;
// 	res["nodes"]         = nodes;
// 	//res["n_nodes"]       = n_nodes;
// 	//res["tier"]          = tier;
// 	//res["terminal_prob"] = terminal_prob;
// 	//res["ran1"]          = ran1;
// 	//res["ran2"]          = ran2;
// 	//res["ran3"]          = ran3;
// 	res["var"]           = var;
// 	res["cut"]           = cut;
// 	res["leaves"]        = leaves;
// 	//res["terminal"]      = terminal;
// 	res["killed"]        = killed;
// 	return res;
// }
