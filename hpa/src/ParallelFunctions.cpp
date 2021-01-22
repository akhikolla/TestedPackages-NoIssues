#include "ParallelFunctions.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;

// Initialize some binary functions

  // calculate square
double sqr_parallel(double x)
{
  return(x * x);
}

  // calculate standard normal cdf
double normalCDF(double x)
{
  return std::erfc(-x / std::sqrt(2)) / 2;
}

// calling ::pow and ::exp directly causes error for
// some systems probably because of confusion between
// Rcpp and std implementations so let's provide wrapper
// functions in order to avoid errors

  // wrapper function to calculate pow
double pow_parallel(double x, int y)
{
  return std::pow(x, y);
}

  // wrapper function to calculate exp
double exp_parallel(double x)
{
  return std::exp(x);
}

  // wrapper function to calculate sqrt
double sqrt_parallel(double x)
{
  return std::sqrt(x);
}

// Parallel functions

// Parallel pow of vectors (struct)
struct ParallelVectorPowStruct : public Worker
{
  // source matrix
  const RVector<double> input;
  
  // powers matrix
  const RVector<double> input_powers;
  
  // destination matrix
  RVector<double> output;
  
  // type of power (0 - general, 1 - square, 2 - square root)
  int pow_type;
  
  // initialize with source and destination
  ParallelVectorPowStruct(const NumericVector input, 
                          NumericVector input_powers, 
                          NumericVector output, int pow_type) 
    : input(input), input_powers(input_powers), 
      output(output), pow_type(pow_type) {}
  
  // take the powers of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    
    if (pow_type == 0)
    {
      // perform pow operation
      std::transform(input.begin() + begin, 
                     input.begin() + end, 
                     input_powers.begin(),
                     output.begin() + begin, 
                     ::pow_parallel);
    }
    
    if (pow_type == 1)
    {
      std::transform(input.begin() + begin, 
                     input.begin() + end,
                     output.begin() + begin, 
                     ::sqr_parallel);
    }
    
    if (pow_type == 2)
    {
      // perform pow operation
      std::transform(input.begin() + begin, 
                     input.begin() + end,
                     output.begin() + begin, 
                     ::sqrt_parallel);
    }
  }
};

// Parallel pow of vector elements
NumericVector ParallelVectorPow(NumericVector x, double value = 1) 
{
  int pow_type = 0;
  
  if (value == 0.5)
  {
    pow_type = 2;
  }
  
  if (value == 2)
  {
    return (x * x);
  }
  
  // allocate the output matrix
  NumericVector output(x.size());
  
  if (value == 1)
  {
    return (x);
  }
  
  if(value == 0)
  {
    std::fill(output.begin(), output.end(), 1);
    return (output);
  }
  
  // allocate the powers matrix
  NumericVector input_powers(x.size());
  std::fill(input_powers.begin(), input_powers.end(), value);
  
  // ParallelVectorPowStruct functor
  ParallelVectorPowStruct parallelVectorPowStruct(x, input_powers, 
                                                  output, pow_type);
  
  // call parallelFor to do the work
  parallelFor(0, x.length(), parallelVectorPowStruct);
  
  // return the output matrix
  return (output);
}

// Parallel exp of vectors
struct ParallelVectorExpStruct : public Worker
{
  // source matrix
  const RVector<double> input;
  
  // destination matrix
  RVector<double> output;
  
  // initialize with source and destination
  ParallelVectorExpStruct(const NumericVector input, NumericVector output) 
    : input(input), output(output) {}
  
  // take the exponents of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
      std::transform(input.begin() + begin, 
                     input.begin() + end,
                     output.begin() + begin, 
                     ::exp_parallel);
  }
};

// Parallel exponent of vector elements
NumericVector ParallelVectorExp(NumericVector x) 
{
  // allocate the output matrix
  NumericVector output(x.size());
  
  // ParallelVectorPowStruct functor
  ParallelVectorExpStruct parallelVectorExpStruct(x, output);
  
  // call parallelFor to do the work
  parallelFor(0, x.length(), parallelVectorExpStruct);
  
  // return the output matrix
  return (output);
}

//' Calculate normal pdf in parallel
//' @description Calculate in parallel for each value from vector \code{x} 
//' density function of normal distribution with 
//' mean equal to \code{mean} and standard deviation equal to \code{sd}.
//' @param x numeric vector of quantiles.
//' @param mean double value.
//' @param sd double positive value.
//' @template is_parallel_Template
//' @template dnorm_parallel_examples_Template
//' @export
// [[Rcpp::export]]
NumericVector dnorm_parallel(NumericVector x, double mean = 0, 
                             double sd = 1, bool is_parallel = false) 
{
  if(!is_parallel)
  {
    return(dnorm(x, mean, sd));
  }
  
  NumericVector result = ParallelVectorPow((x - mean) / sd, 2);
  result = -0.5 * result;
  result = ParallelVectorExp(result);
  result = result / (sd * std::sqrt(2 * M_PI));
  
  return(result);
}

// Parallel cdf vectors
struct ParallelVectorNormalCDFStruct : public Worker
{
  // source matrix
  const RVector<double> input;
  
  // destination matrix
  RVector<double> output;
  
  // initialize with source and destination
  ParallelVectorNormalCDFStruct(const NumericVector input, 
                                NumericVector output) 
    : input(input), output(output) {}
  
  // take the normal cdfs of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    
    std::transform(input.begin() + begin, 
                   input.begin() + end,
                   output.begin() + begin, 
                   ::normalCDF);
  }
};

//' Calculate normal cdf in parallel
//' @description Calculate in parallel for each value from vector \code{x} 
//' distribution function of normal distribution with 
//' mean equal to \code{mean} and standard deviation equal to \code{sd}.
//' @param x vector of quantiles: should be numeric vector,
//' not just double value.
//' @param mean double value.
//' @param sd double positive value.
//' @template is_parallel_Template
//' @export
// [[Rcpp::export]]
NumericVector pnorm_parallel(NumericVector x, double mean = 0, 
                             double sd = 1, bool is_parallel = false)
{
  if(!is_parallel)
  {
    return(pnorm(x, mean, sd));
  }
  
  // allocate the output matrix
  NumericVector output(x.size());
  
  // ParallelVectorNormalCDFStruct functor
  ParallelVectorNormalCDFStruct parallelVectorNormalCDFStruct((x - mean) / 
                                                              sd, output);
  
  // call parallelFor to do the work
  parallelFor(0, x.length(), parallelVectorNormalCDFStruct);
  
  // return the output matrix
  return (output);
}
