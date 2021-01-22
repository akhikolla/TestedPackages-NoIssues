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

#include <iomanip>      // std::setw
#include "lasso.hpp"

using std::setw;

void Lasso::printResults() {
    cout << setw(15) << "TABLE"
         << setw(15) << "df"
         << setw(15) << "rsq"
         << setw(15) << "lambda"
         << setw(15) << "L1norm" << endl;
    for (size_t i = 0; i < this->lambda.size(); i++) {
        cout << setw(15) << "TABLE"
             << setw(15) << this->df[i]
             << setw(15) << this->devRatio[i]
             << setw(15) << this->lambda[i]
             << setw(15) << this->L1norm[i] << endl;
    }

    cout << setw(15) << "BETA"
         << setw(15) << "X_i";
    for (size_t j = 0; j < this->lambda.size(); j++) {
        string colName = string("beta") + std::to_string(j);
        cout << setw(15) << colName;
    }
    cout << endl;

    for (size_t i = 0; i < this->nVars_; i++) {
        double tmpSum = 0;
        for (size_t j = 0; j < this->lambda.size(); j++) {
            tmpSum += this->beta[j][i];
        }
        if (tmpSum > 0) {
            cout << setw(15) << "BETA"
                 << setw(15) << i;
            for (size_t j = 0; j < this->lambda.size(); j++) {
                cout << setw(15) << beta[j][i];
            }
            cout << endl;
        }
    }
    cout << endl;
}

