/* Copyright (C) 2014 - 2018  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of simts R Methods Package
 *
 * The `simts` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `simts` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

#ifndef COVARIANCE_MATRIX
#define COVARIANCE_MATRIX

arma::field<arma::mat> compute_cov_cpp(arma::field<arma::vec> signal_modwt, unsigned int nb_level, std::string compute_v,
                                        bool robust, double eff);

arma::mat fast_cov_cpp(const arma::vec& ci_hi, const arma::vec& ci_lo);

#endif
