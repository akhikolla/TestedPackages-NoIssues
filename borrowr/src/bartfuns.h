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

#ifndef GUARD_bartfuns_h
#define GUARD_bartfuns_h

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include <algorithm>

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, int* nc);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots);
//--------------------------------------------------
//compute n and \sum y_i for left and right give bot and v,c
void getsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);
//--------------------------------------------------
//lh, replacement for lil that only depends on sum y.
double lh(size_t n, double sy, double sigma, double gamma);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else a/(1+d)^b
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi);
//--------------------------------------------------
//compute n and \sum y_i for left and right bots
void getsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<size_t>& nv, std::vector<double>& syv);
//--------------------------------------------------
// draw all the bottom node mu's
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);
//--------------------------------------------------
//birth proposal
void bprop(tree& x, xinfo& xi, pinfo& pi, tree::npv& goodbots, double& PBx, tree::tree_p& nx, size_t& v, size_t& c, double& pr, std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen);
//--------------------------------------------------
// death proposal
void dprop(tree& x, xinfo& xi, pinfo& pi, tree::npv& goodbots, double& PBx, tree::tree_p& nx, double& pr, rn& gen);
//--------------------------------------------------
//draw one mu from post
double drawnodemu(size_t n, double sy, double gamma, double sigma, rn& gen);
//--------------------------------------------------
//draw variable splitting probabilities from Dirichlet (Linero, 2018)
void draw_s(std::vector<size_t>& nv, std::vector<double>& lpv, double& theta, rn& gen);
//--------------------------------------------------
//draw Dirichlet sparsity parameter from posterior using grid
void draw_theta0(bool const_theta, double& theta, std::vector<double>& lpv,
		 double a, double b, double rho, rn& gen);
#endif
