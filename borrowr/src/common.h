/*
 *  borrowr: estimate population average treatment effects with borrowing between data sources.
 *  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
 *  This code is a modified version from the BART R package from April 2019.
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

#ifndef GUARD_common_h
#define GUARD_common_h

#ifdef MATHLIB_STANDALONE
#define NoRcpp
#else
#define RNG_Rcpp
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef NoRcpp

#include <stdio.h> // for printf

using std::cout;

#define PI 3.1415926535897931

#else // YesRcpp

#include <Rcpp.h>

#define printf Rprintf
#define cout Rcpp::Rcout

#endif

// log(2*pi)
#define LTPI 1.83787706640934536

#include "rn.h"

#endif
