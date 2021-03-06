// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// h11_DeepState_TestHarness_generation.cpp and h11_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <ctime>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector h11(NumericMatrix x, NumericVector t, NumericVector theta, NumericVector integrControl);

TEST(smam_deepstate_test,h11_test){
  RInside R;
  std::time_t t = std::time(0);
  std::cout << "input starts" << std::endl;
  NumericMatrix x  = RcppDeepState_NumericMatrix();
  std::string x_t = "/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/h11/AFL_h11/afl_inputs/" + std::to_string(t) + "_x.qs";
  qs::c_qsave(x,x_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector t  = RcppDeepState_NumericVector();
  std::string t_t = "/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/h11/AFL_h11/afl_inputs/" + std::to_string(t) + "_t.qs";
  qs::c_qsave(t,t_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "t values: "<< t << std::endl;
  NumericVector theta  = RcppDeepState_NumericVector();
  std::string t_t = "/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/h11/AFL_h11/afl_inputs/" + std::to_string(t) + "_t.qs";
  std::string theta_t = "/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/h11/AFL_h11/afl_inputs/" + std::to_string(t) + "_theta.qs";
  qs::c_qsave(theta,theta_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "theta values: "<< theta << std::endl;
  NumericVector integrControl  = RcppDeepState_NumericVector();
  std::string integrControl_t = "/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/h11/AFL_h11/afl_inputs/" + std::to_string(t) + "_integrControl.qs";
  qs::c_qsave(integrControl,integrControl_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "integrControl values: "<< integrControl << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    h11(x,t,theta,integrControl);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
