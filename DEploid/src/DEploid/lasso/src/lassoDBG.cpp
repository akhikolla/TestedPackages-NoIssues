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

#include <assert.h>     // assert
#include <string>
#include <fstream>      // std::ifstream
#include <limits>       // std::numeric_limits<double>::infinity();
#include <iomanip>      // std::setw
#include "dEploidLasso.hpp"


bool Lasso::print_normalized_struff() {
    dout <<"First 5 elements of normalized x 5 variables" << endl;
    for (size_t i = 0; i < min((size_t)5, nVars_); i++) {
        vector<double> xx = this->standardized_x_transposed[i];
        for (size_t j = 0; j < min((size_t)5, xx.size()); j++) {
            dout << std::setw(9) << xx[j] << ", ";
        }
    }
    dout << endl;

    dout << "First 5 elements of normalized y" << endl;
    for (size_t i = 0; i < min((size_t)5, nObs_) ; i++) {
        dout << std::setw(8) << this->standardized_y[i] << ", ";
    }
    dout << endl;

    return(true);
}


bool Lasso::print_initial_gk() {
    dout << "Initial gk:" << endl;
    for (size_t i = 0; i < min((size_t)5, nVars_); i++) {
        dout << g[i] << ", ";
    }
    dout << endl;

    return(true);
}


bool Lasso::print_homogeneous_input() {
    // DEBUG MESSAGE
    dout << "Variables: ";
    for (size_t i = 0; i < this->nVars_; i++) {
        if (ju[i] == 0) {dout << i <<", ";}
    }
    dout << " are homogeneous vectors." << endl;

    return(true);
}
