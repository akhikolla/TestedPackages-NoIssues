#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
/*
 * c_GiniMd.cpp
 * 
 * Copyright 2020 Joern <joern@joern-daheim>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */
// [[Rcpp::export]]
double c_gmd(Rcpp::NumericVector j) {
  unsigned long i;
  double sumX = 0, avg = 0, gm = 0, gmd;
  
  std::size_t n = j.size() ;
  std::sort(j.begin(), j.end());sumX = sum(j);
  avg = sumX / n;
  for (i = 0; i < n; i++) {
    gm += (i + 1) * j[i];
  }
  gmd = 4 * gm / (n * (n - 1)) - 2 * avg * (n + 1) / (n - 1);
  
  return (gmd);
} 

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
