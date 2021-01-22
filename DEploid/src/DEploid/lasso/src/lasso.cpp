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
#include <cmath>        // std::abs
#include <limits>       // std::numeric_limits<double>::infinity();
#include <algorithm>    // for min
#include <iomanip>      // std::setw
#include "lasso.hpp"

using std::setw;
// using std::abs;  // THIS IS VERY IMPORTANT! without this line, abs gives int
// NOW using fabs instead

lasso::TxtReader::TxtReader(const char inchar[]) {
    string fileName(inchar);
    ifstream in_file(inchar);
    string tmp_line;
    if (in_file.good()) {
        getline(in_file, tmp_line);
        while (tmp_line.size() > 0) {
            size_t field_start = 0;
            size_t field_end = 0;
            size_t field_index = 0;
            vector <double> contentRow;
            while (field_end < tmp_line.size()) {
                field_end = min(min(tmp_line.find(',', field_start),
                                    tmp_line.find('\t', field_start)),
                                    tmp_line.find('\n', field_start));
                string tmp_str = tmp_line.substr(field_start,
                                                 field_end - field_start);
                contentRow.push_back(strtod(tmp_str.c_str(), NULL));
                if (field_index == 0) {
                    vec.push_back(strtod(tmp_str.c_str(), NULL));
                }
                field_start = field_end+1;
                field_index++;
            }
            this->matrix.push_back(contentRow);
            getline(in_file, tmp_line);
        }
    }

    in_file.close();
}


standardizeVector::standardizeVector(vector <double> vec) {
    assert(ret.size() == 0);
    size_t nObs_ = vec.size();
    this->mean = lasso::sumOfVec(vec) / static_cast<double>(nObs_);
    vector <double> mean_vec = vector <double> (nObs_, mean);
    vector <double> vec_diff = lasso::vecDiff(vec, mean_vec);  // (y-ym)

    vector <double> v_vec = vector <double> (nObs_,
                                         1.0/sqrt(static_cast<double>(nObs_)));
    vector <double> v_times_vec_diff = lasso::vecProd(vec_diff, v_vec);
    // y=v*(y-ym)

    vector <double> tmpProd = lasso::vecProd(v_times_vec_diff,
                                             v_times_vec_diff);
    this->variance = lasso::sumOfVec(tmpProd);
    // this->variance = sumOfVec(tmpProd)/(nObs_);
    this->stdv = sqrt(this->variance);

    for (double tmp : v_times_vec_diff) {
        ret.push_back(tmp/stdv);
    }
    assert(ret.size() == nObs_);

    //// DEBUG
    // dout << "sum of x = " << sumOfVec(ret) << endl;
    // vector <double> tmpdbg = vecProd(ret,ret);
    // double tmpSum = 0.0;
    // for ( double tmp : ret ){
        // tmpSum += tmp*tmp;
    // }
    // dout << "1/N sum of x_sq " << tmpSum/(double)nObs_ <<endl;
}


Lasso::Lasso(const vector < vector <double> > &x,
    const vector < double > &wsaf, size_t nLambda) {
    this->nVars_ = x[0].size();
    this->nObs_ = x.size();

    // TODO(JOE_ZHU): check for x size
    dout<< "Matrix size = "<< this->nObs_ << " " << this->nVars_ << endl;
    dout<< "Vector length = " << wsaf.size() << endl;

    if (this->nObs_ != wsaf.size()) {
        throw lasso::InvalidInput("Lasso matrix and vector size incompatible!");
    }

    // Initialize

    this->initialization(nLambda);
    this->checkVariables(x);
    this->standarization(x, wsaf);
    this->productOfxy();
    this->computeNullDev(x, wsaf);
    for (size_t i = 0; i < this->lambda.size(); i++) {
        this->setLambdaCurrent(1.0 / (3.0+static_cast<double>(i)));
        dout << endl << "****************** current lambda: "
             << this->lambdaCurrent() << endl;
        this->lassoGivenLambda();
        dout << "######### lassoGivenLambda finished at " << npass_;
        this->setLambdaPrevious(this->lambdaCurrentScaled());
        // dout << " this->ninCurrent_ "<< this->ninCurrent_ << endl;
        // FETCH AND UPDATE THE CURRENT INFERENCE RESULTS
        this->beta.push_back(betaCurrent);
        this->lambda[i] = lambdaCurrent();
        this->intercept[i] = interceptCurrent();
        this->devRatio[i] = rsqCurrent();
        this->df[i] = dfCurrent();
    }
    dout << "beta.size = " << beta.size() <<endl;
    this->computeL1Norm();
}


Lasso::~Lasso() {}


void Lasso::productOfxy() {
    this->g = vector <double> (nVars_, 0.0);
    for (size_t j = 0; j < this->nVars_; j++) {
        // skip homogeneous variables
        if (this->ju[j] != 0) {
            this->g[j] = computeGk_abs(standardized_y,
                                       standardized_x_transposed[j]);
        }
    }
    assert(this->g.size() == nVars_);
    assert(this->print_initial_gk());
}


void Lasso::checkVariables(const vector < vector <double> > &x) {
    // check for homogeneous vectors
    this->ju = vector <double> (this->nVars_, 0.0);
    for (size_t i = 0; i < this->nVars_; i++) {
        for (size_t ii = 1; ii < this->nObs_; ii++) {
            if (x[ii][i] != x[0][i]) {
                ju[i] = 1.0;
                break;
            }
        }
    }
    assert(print_homogeneous_input());
}


void Lasso::standarization(const vector < vector <double> > &x,
                                  const vector < double > &y) {
    // standarize x
    assert(standardized_x_transposed.size() == 0);
    assert(x_stdv.size() == 0);
    assert(x_mean.size() == 0);
    for (size_t i = 0; i < nVars_; i++) {
        // Extract the ith variable
        vector <double> var_i;
        for (size_t j = 0; j < this->nObs_; j++) {
            var_i.push_back(x[j][i]);
        }
        standardizeVector vecX(var_i);
        x_stdv.push_back(vecX.stdv);
        x_mean.push_back(vecX.mean);
        x_variance.push_back(1.0);  // set the variance to 1
        standardized_x_transposed.push_back(vecX.ret);
    }
    assert(standardized_x_transposed.size() == nVars_);

    // standarize y
    standardizeVector vecY(y);
    standardized_y = vecY.ret;
    this->y_stdv = vecY.stdv;
    this->y_mean = vecY.mean;

    assert(print_normalized_struff());
}


void Lasso::initialization(size_t nLambda) {
    this->nulldev_ = 0.0;
    this->betaCurrent = vector <double> (this->nVars_, 0.0);

    this->lambda = vector <double> (nLambda, 0.0);
    this->intercept = vector <double> (nLambda, 0.0);
    this->devRatio = vector <double> (nLambda, 0.0);
    this->df = vector <int> (nLambda, 0);
    this->L1norm = vector <double> (nLambda, 0.0);

    this->ix = vector <double> (nVars_, 0.0);
    this->coefficentCurrent = vector <double> (this->nVars_, 0.0);

    // for (size_t i = 0; i < nLambda; i++){
        // this->lambda.push_back(0.0);
        // this->intercept.push_back(0.0);
        // this->devRatio.push_back(0.0);
        // this->df.push_back((int)0);
    // }

    this->setRsqCurrent(0.0);

    // To initialize to nVars, as index are 0-based, index is from 0 to nVars-1
    this->indexArray = vector <size_t> (nVars_, nVars_);

    this->lowerLimit = 0;
    // this->lowerLimit = -std::numeric_limits<double>::infinity();
    this->upperLimit = std::numeric_limits<double>::infinity();

    this->maxIteration_ = 100000;
    this->thresh_ = 1e-7;
    this->dfmax_ = nVars_ + 1;
    this->mm = vector < size_t > (nVars_, (size_t)0);
    this->iz = 0;

    this->setLambdaPrevious(0.0);
    this->npass_ = 0;
    this->setDfCurrent(dfmax_);
    this->nin = 0;
}


void Lasso::computeNullDev(const vector < vector <double> > &x,
                                  const vector < double > &wsaf) {
    double ybar = lasso::sumOfVec(wsaf) / static_cast<double>(wsaf.size());
    vector <double> ybar_vec = vector <double> (wsaf.size(), ybar);
    vector <double> diff = lasso::vecDiff(wsaf, ybar_vec);
    vector <double> tmpSq = lasso::vecProd(diff, diff);
    this->nulldev_ = lasso::sumOfVec(tmpSq);
    dout << "nulldev = " << this->nulldev_ << endl;
}


void Lasso::lassoGivenLambda() {
    /* USE THE FOLLOWING VARIABLES
     *
     * lambdaCurrent()
     * standardized_x_transposed
     * standardized_y
     *
     */

    /*
     * LOCAL INITIALIZATION
     */
    // double rsq=0.0; //TODO, just replace as rsq0

    // ulam is user defined lambdas...
    // no = number of observations
    // ni = number of predictor variables

    // c   parm = penalty member index (0 <= parm <= 1)
    // c        = 0.0 => ridge
    // c        = 1.0 => lasso
    // double beta = 1.0;
    // double bta = beta; // lasso part
    // double omb=1.0-bta; // ridge part, omb = 0
    // double flmin = 1.0; // this is defined in glmnet.R

    this->setDfCurrent(0);
    this->jz = 1;
    this->setLambdaCurrentScaled(lambdaCurrent()/this->y_stdv);
    // tlam=bta*(2.0*alm-alm0)
    // beta, bta = 1
    double tlam = 1.0*(2.0 * lambdaCurrentScaled() - lambdaPrevious());

    this->chooseVariables(tlam);

    this->updatingCore();

    // Map coefficients, coefficient -> beta
    this->coefficentToBeta();

    // Rescale coefficients
    this->rescaleCoefficents();

    // Compute intercept
    this->computeIntercept();
}

void Lasso::chooseVariables(double tlam) {
    dout << "Choose variables, where tlam = " << tlam <<endl;

    for (size_t k = 0; k < this->nVars_; k++) {
        if (ix[k] == 1) {continue;}
        if (ju[k] == 0) {continue;}
        if (g[k] > (tlam*1.0)) {
            dout << "  * will need vairable " << k
                 << ", where g[k] = " << g[k] << " > tlam." << endl;
            ix[k] = 1.0;
        }
    }
    dout << endl;
}


double Lasso::rechooseVariables() {
    dout << "Choose variables, where lambdaCurrentScaled = "
         << lambdaCurrentScaled() <<endl;

    double ixx = 0;
    for (size_t k = 0; k < nVars_; k++) {
        if (ix[k] == 1) {continue;}
        if (ju[k] == 0) {continue;}

        g[k] = this->computeGk_abs(standardized_y,
                                   standardized_x_transposed[k]);
        if (g[k] > lambdaCurrentScaled()) {
            dout << "  * will need vairable " << k
                 << ", where g[k] = " << g[k]
                 << " > lambdaCurrentScaled()." << endl;
            ix[k] = 1;
            ixx = 1;
        }
    }
    dout << endl;
    return ixx;
}


void Lasso::updateWithNewVariables() {
    dout << "Update with new variables" <<endl;
    dout << "###### begin scanning variables ##########" <<endl;
    this->npass_++;
    double dlx = 0.0;

    for (size_t k = 0; k < nVars_; k++) {
        if (ix[k] == 0) {
            dout << "  * Skipping variable " << k << endl;
            continue;
        }

        dout << "  * Variable " << k << endl;
        double ak = this->coefficentCurrent[k];
        double gk = computeGk(this->standardized_y,
                              this->standardized_x_transposed[k]);
        this->updateCoefficient(k, ak, gk);

        // If the coefficient unchanged, move on to the next
        if ( coefficentCurrent[k] == ak ) {
            continue;
        }

        // Updating variables array
        if (mm[k] == 0) {
          this->indexArray[nin] = k;
          dout << "  ** Include indexArray[" << nin
               << "] with variable " << k << endl;

          dout << "  *** Number of choose variables mm[k] was " << mm[k];
          nin += 1;
          mm[k] = nin;
          dout << ", updated to " << mm[k] << endl;
          if (nin > nVars_) {
              dout << "!!!NUMBER OF MAXIMUM VARIABLE REACHED" << endl;
              break;
          }
        }
        dout << "  *** current nin = " << nin << ", at Variable " << k << endl;

        double del = this->updateYReturnDel(k, gk, ak);

        dout << " ** Convergence check, dlx was " << dlx;
        dlx = max(x_variance[k]*del*del, dlx);
        dout << " updated to " << dlx << endl;
    }

    dout << "###### finish scanning variables ##########" <<endl;

    if (dlx >= this->thresh_) {
        updateWithTheSameVariables();
    }
}


void Lasso::updatingCore() {
    double ixx = 0;
    bool keepUpdating = true;
    while (keepUpdating) {
        if (iz*jz != 0) {
            iz = 1;
            updateWithTheSameVariables();
        } else {
            updateWithNewVariables();

            if (nin > nVars_) {
                keepUpdating = false;
                break;
            }

            ixx = this->rechooseVariables();

            if (ixx != 1) {
                keepUpdating = false;
                break;
            }
            if (npass_ > maxIteration_) {
                keepUpdating = false;
                break;
            }
        }
    }
}


void Lasso::updateWithTheSameVariables() {
    bool keepUpdating = true;
    dout << endl;
    dout << "Update With The Same Variables" <<endl;

    while (keepUpdating) {
        this->npass_++;
        double dlx = 0.0;
        dout << "###### begin scanning variables ##########" << endl;
        for (size_t l = 0; l < this->nin; l++) {
            size_t k = indexArray[l];
            dout << "  * Current variable " << k << endl;
            double ak = this->coefficentCurrent[k];
            double gk = computeGk(this->standardized_y,
                                  this->standardized_x_transposed[k]);

            this->updateCoefficient(k, ak, gk);

            if (coefficentCurrent[k] == ak) {
                continue;
            }
            double del = this->updateYReturnDel(k, gk, ak);
            dout << "  ** Convergence check, dlx was " << dlx;
            dlx = max(x_variance[k]*del*del, dlx);
            dout << " updated to " << dlx << endl;
        }
        dout << "###### finish scanning variables ##########" << endl;

        if (dlx < this->thresh_) {
            keepUpdating = false;
        }

        if (npass_ > this->maxIteration_) {
            keepUpdating = false;
            break;
        }
    }
    jz = 0;
}


void Lasso::updateCoefficient(size_t k, double previousCoefficentValue,
                                     double gk) {
    double u = gk + previousCoefficentValue * x_variance[k];
    // u=gk+ak*xv(k)
    // u=gk+ak*xv(k)
    double v = fabs(u) - 1.0 * lambdaCurrentScaled();
    // v=fabs(u)-vp(k)*ab
    // v=fabs(u)-vp(k)*ab

    dout << "  ** Update coefficient" << endl;
    dout << "     Check if " << lambdaCurrentScaled()
         << " is in the radius of u = " << u << ", it gives v = " << v <<endl;
    dout << "     Updating coefficient: " << previousCoefficentValue << " to ";
    coefficentCurrent[k] = 0.0;
    double u_sign = 1.0;
    if (u < 0) {
        u_sign = -1.0;
    }
    if (v > 0.0) {
       coefficentCurrent[k] = max(lowerLimit,
                                   min(upperLimit, v*u_sign/x_variance[k]));
    }
    dout << coefficentCurrent[k] <<endl;
}


double Lasso::updateYReturnDel(size_t k, double gk,
                                      double previousCoefficentValue) {
    double del = this->coefficentCurrent[k] - previousCoefficentValue;
    // dout << "del = " << del ;
    this->rsqCurrent_ += del * (2.0 * gk - del * x_variance[k]);
    dout << "  ** Current rsq is " << this->rsqCurrent()
         << ", with updated ys:";
    for (size_t i = 0; i < this->nObs_; i++) {
        standardized_y[i] -= del * standardized_x_transposed[k][i];
        if (i < 5) {dout << standardized_y[i] << ",";}
    }
    dout << endl;
    return del;
}


void Lasso::computeIntercept() {
    this->setInterceptCurrent(0.0);

    if (this->nin == 0) {
        return;
    }

    double y_remaining = this->y_mean;
    for (size_t i = 0; i < (size_t)this->nin; i++) {
        size_t k = indexArray[i];
        y_remaining -= betaCurrent[k] * x_mean[k];
        // dout << "y_remaining "<<y_remaining<<endl;
    }
    this->setInterceptCurrent(y_remaining);
}


void Lasso::rescaleCoefficents() {
    this->dfCurrent_ = 0;
    for (size_t i = 0; i < (size_t)this->nin; i++) {
        size_t k = indexArray[i];
        if (this->betaCurrent[k] > 0) {
            this->dfCurrent_++;
        }
        this->betaCurrent[k] *= y_stdv;
        this->betaCurrent[k] /= x_stdv[k];
    }
}


void Lasso::coefficentToBeta() {
    this->betaCurrent = this->coefficentCurrent;

    // for ( size_t i = 0; i < this->nVars_; i++){
        // this->betaCurrent[i] = this->coefficentCurrent[i];
    // }
}


double Lasso::computeGk(const vector<double> &y,
                               const vector<double> &x) {
    vector <double> gk_vec = lasso::vecProd(y, x);
    return lasso::sumOfVec(gk_vec);
}


double Lasso::computeGk_abs(const vector<double> &y,
                                   const vector<double> &x) {
    return fabs(computeGk(y, x));
}


void Lasso::computeL1Norm() {
    for (size_t i = 0; i < beta.size(); i++) {
        this->L1norm[i] = lasso::sumOfVec(beta[i]);
    }
}
