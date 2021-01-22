/*
 *  borrowr: estimate population average treatment effects with borrowing between data sources.
 *  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List matchesToCor(NumericMatrix x) {
  List res;

  int n = x.ncol();
  int m = x.nrow();

  NumericMatrix R(n, n);

  for(int i = 0; i < n; i++) {
    R(i, i) = 1;
  }

  for(int i = 0; i < (n - 1); i++) {
    for(int j = i + 1; j < n; j++) {
      double numMatches = 0;
      for(int k = 0; k < m; k ++) {
        if(x(k, i) == x(k, j)) {
          numMatches += 1;
        }
      }
      R(i, j) = numMatches / m;
      R(j, i) = numMatches / m;
    }
  }

  //res["res"] = x;
  res["n"]   = n;
  res["m"]   = m;
  res["R"]   = R;

  return res;
}

/*** R
set.seed(123)
foo <- matrix(sample(1:3, 3 * 5, replace = TRUE), 3, 5)
out <- matchesToCor(foo)
out$R[1:8, 1:8]
*/
