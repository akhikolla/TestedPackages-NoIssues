/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of DEploid-r.
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


#ifndef dEploidr_r_random_generator
#define dEploidr_r_random_generator

#include <Rcpp.h>
#include <memory>
#include "random_generator.hpp"

using namespace Rcpp;

class RRandomGenerator : public RandomGenerator {
 public:
  RRandomGenerator(std::shared_ptr<FastFunc> ff):RandomGenerator(ff) {
    this->initializeUnitExponential();
    size_t new_seed = this->sample()*1073741824; //2^30
    this->set_seed(new_seed);
  };
  virtual ~RRandomGenerator() {};

  void initialize() {};

  double sample() {
    RNGScope scope;
    return R::runif(0,1);
  }

  double sampleUnitExponential() {
    RNGScope scope;
    return R::rexp(1);
  }
};

#endif
