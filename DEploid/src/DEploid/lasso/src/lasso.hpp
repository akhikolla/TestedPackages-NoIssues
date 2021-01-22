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
#include <exception>
#include <algorithm>
#include "dbgmacro.hpp"


#ifndef LASSO
#define LASSO

using std::vector;
using std::exception;
using std::string;
using std::endl;
using std::cout;
using std::ifstream;
using std::min;
using std::max;

namespace lasso {
struct InvalidInput : std::exception {
    string throwMsg;

    InvalidInput() {
    }

    explicit InvalidInput(string str) {
        this->throwMsg      = "\033[1;31m" + str + "\033[0m";
    }
    virtual ~InvalidInput() throw() {}
    virtual const char* what() const noexcept {
        return throwMsg.c_str();
    }
};


template <typename T>
vector <T> matrixTimesVec(const vector < vector < T > > &x,
                          const vector <double> &b) {
    vector <double> ret(x.size(), 0.0);
    for (size_t i = 0; i < x.size(); i++) {
        for (size_t k = 0; k < x[i].size(); k++) {
            ret[i] += x[i][k] * b[k];
        }
    }
    return ret;
}


template <typename T>
vector <T> vecDiff(const vector<T> &vecA,
                   const vector<T> &vecB) {
    // assert(vecA.size() == vecB.size());
    vector <T> difference(vecA.size(), (T)0);
    for (size_t i = 0; i < vecA.size(); i++) {
        difference[i] = vecA[i] - vecB[i];
    }
    return difference;
}


template <typename T>
vector <T> vecSum(const vector<T> &vecA, const vector<T> &vecB ) {
    assert(vecA.size() == vecB.size());
    vector <T> tmpSum(vecA.size(), (T)0);
    for (size_t i = 0; i < vecA.size(); i++) {
        tmpSum[i] = vecA[i] + vecB[i];
    }
    return tmpSum;
}


template <typename T>
vector <T> vecProd(const vector<T> &vecA, const vector<T> &vecB) {
    // assert(vecA.size() == vecB.size());
    vector <T> tmpProd(vecA.size(), (T)0);
    for (size_t i = 0; i < vecA.size(); i++) {
        tmpProd[i] = vecA[i] * vecB[i];
    }
    return tmpProd;
}


template <typename T>
T sumOfVec(const vector <T>& array) {
    T tmp = 0;
    for (auto const& value : array) {
        tmp += value;
    }
    return tmp;
}

class TxtReader{
  friend class Lasso;

 public:
    vector < vector < double > > matrix;
    vector < double > vec;
    explicit TxtReader(const char inchar[]);
    ~TxtReader() {}

 private:
};
}  // namespace lasso

struct standardizeVector {
    vector <double> ret;
    double mean;
    double stdv;
    double variance;

    explicit standardizeVector(vector <double> vec);
    ~standardizeVector(){}
};


class Lasso{
#ifdef UNITTEST
  friend class TestDEploidLASSO;
#endif
  friend class DEploidLASSO;
  friend class DEploidIO;
 public:
    // Lasso();
    Lasso(const vector < vector <double> > &x,  // nObs x nVariable
                 const vector < double > &y,
                 size_t nLambda = 100);
    ~Lasso();
    void printResults();

 private:
    // FUNCTIONS
    // COMMON
    void initialization(size_t nLambda);
    void standarization(const vector < vector <double> > &x,
                        const vector < double > &y);
    void checkVariables(const vector < vector <double> > &x);
    void productOfxy();
    void computeL1Norm();
    void computeNullDev(const vector < vector <double> > &x,
                        const vector < double > &y);

    // FOR EACH LAMBDA UPDATE
    void lassoGivenLambda();

    // VARIABLES, GETTERS AND SETTERS
    // COMMON
    // OUTPUT
    vector < vector <double> > beta;  // nLambda x nVars
    vector < double > lambda;  // size of nLambda
    vector < double > devRatio;
    vector < double > intercept;
    vector < int > df;
    vector < double> L1norm;

    // OTHER VARIABLES
    // DATA RELATED
    size_t nObs_;
    size_t nVars_;
    vector < vector <double> > standardized_x_transposed;  // nVariable x nObs
    vector <double> standardized_y;
    vector <double> x_mean;
    vector <double> x_stdv;
    vector <double> x_variance;
    double y_stdv;
    double y_mean;
    double nulldev_;

    // COMPUTATION RELATED
    vector < size_t > indexArray;
    vector < size_t > mm;  // indicator, that kth variable is already in use
    size_t nin;  // number of variables in use
    int maxIteration_;
    double thresh_;
    vector <double> ju;
    vector <double> g;
    size_t dfmax_;
    vector <double> ix;
    int npass_;
    double lowerLimit;
    double upperLimit;

    // FOR EACH LAMBDA UPDATE
    // VARIABLES
    vector < double > betaCurrent;  // size of nVars
    vector < double > coefficentCurrent;  // size of nVars

    double lambdaCurrent_;
    void setLambdaCurrent(const double setTo) {this->lambdaCurrent_ = setTo;}
    double lambdaCurrent() const {return this->lambdaCurrent_;}

    double lambdaPrevious_;
    void setLambdaPrevious(const double setTo) {this->lambdaPrevious_ = setTo;}
    double lambdaPrevious() const {return this->lambdaPrevious_;}

    double lambdaCurrentScaled_;
    void setLambdaCurrentScaled(const double setTo) {
        this->lambdaCurrentScaled_ = setTo;}
    double lambdaCurrentScaled() const {return this->lambdaCurrentScaled_;}

    double rsqCurrent_;
    void setRsqCurrent(const double setTo) {this->rsqCurrent_ = setTo;}
    double rsqCurrent() const { return this->rsqCurrent_; }

    double interceptCurrent_;
    void setInterceptCurrent(const double setTo) {
        this->interceptCurrent_ = setTo;}
    double interceptCurrent() const { return this->interceptCurrent_; }

    int dfCurrent_;
    void setDfCurrent(const int setTo) {this->dfCurrent_ = setTo;}
    int dfCurrent() const {return this->dfCurrent_;}

    // int ninCurrent_;
    int iz, jz;
    // FUNCTIONS
    void computeIntercept();
    void rescaleCoefficents();
    void coefficentToBeta();
    void updateCoefficient(size_t k, double previousCoefficentValue, double gk);
    void updateWithNewVariables();
    void updateWithTheSameVariables();
    void updatingCore();
    void chooseVariables(double tlam);
    double updateYReturnDel(size_t k, double gk, double ak);
    double computeGk(const vector<double> &y, const vector<double> &x);
    double computeGk_abs(const vector<double> &y, const vector<double> &x);
    double rechooseVariables();

    // Debug tools
    bool print_normalized_struff();
    bool print_initial_gk();
    bool print_homogeneous_input();
};
#endif  // LASSO
