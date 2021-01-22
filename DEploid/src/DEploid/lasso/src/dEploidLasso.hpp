/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample. DEploid-Lasso-lib is a submodule for
 * choosing the appropriate reference panel using the LASSO algorithm.
 *
 * Copyright (C) 2018 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of DEploid-Lasso-lib.
 *
 * dEploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "lasso.hpp"

#ifndef DEPLOIDLASSO
#define DEPLOIDLASSO

class DEploidLASSO : public Lasso {
#ifdef UNITTEST
  friend class TestDEploidLASSO;
#endif
  friend class DEploidIO;

 public:
    // DEploidLASSO();
    DEploidLASSO(const vector < vector <double> > &x,  // nObs x nVariable
                 const vector < double > &y,
                 size_t nLambda = 100);
    ~DEploidLASSO() {}

 private:
    size_t choiceAt_;
    void setChoiceAt(const double setTo) {this->choiceAt_ = setTo;}
    double choiceAt() const { return this->choiceAt_; }
    vector <size_t> choiceIdx;
    vector <double> choiceBeta;

    // selections
    void determineTheCutOff();
    void shrinkThePanel(const vector < vector <double> > &x);
    vector < vector <double> > reducedPanel;
};


#endif
