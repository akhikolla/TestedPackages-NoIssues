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

#ifndef DWT
#define DWT

arma::field<arma::vec> brick_wall(arma::field<arma::vec> x,  
                                  arma::field<arma::vec> wave_filter, 
                                  std::string method = "modwt") ;

arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name = "haar", 
                               unsigned int nlevels = 4, std::string boundary = "periodic", bool brickwall = true);
                                 
arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name = "haar", 
                                 unsigned int nlevels = 4, std::string boundary = "periodic", bool brickwall = true);
#endif
