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

#ifndef GUARD_treefuns_h
#define GUARD_treefuns_h

#include "tree.h"

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi);
//--------------------------------------------------
//evaluate tree tr on grid xi, write to os
void grm(tree& tr, xinfo& xi, std::ostream& os);
//--------------------------------------------------
//fit tree at matrix of x, matrix is stacked columns x[i,j] is *(x+p*i+j)
void fit(tree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv);
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);

#endif
