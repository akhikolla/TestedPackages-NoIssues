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

#include <iostream>  // std::cout
#include "lasso.hpp"
#include "dEploidLasso.hpp"

int main(int argc, char *argv[]) {
    try {
        vector < vector <double> > matrix = lasso::TxtReader(argv[1]).matrix;
        vector < double > wsaf = lasso::TxtReader(argv[2]).vec;
        // vector < vector <double> > matrix =
                // TxtReader ("data/panel_chrom1.txt").matrix;
        // vector < double > wsaf = TxtReader ("data/PG0402-C_chrom1.wsaf").vec;
        // vector < vector <double> > matrix =
                // TxtReader ("data/myX.txt").matrix;
        // vector < double > wsaf = TxtReader ("data/myy.txt").vec;
        Lasso dummy(matrix, wsaf, 250);
        dummy.printResults();
        matrix.clear();
        wsaf.clear();
        return EXIT_SUCCESS;
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
