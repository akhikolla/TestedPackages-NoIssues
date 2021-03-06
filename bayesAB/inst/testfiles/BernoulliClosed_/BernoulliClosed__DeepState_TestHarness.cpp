// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// BernoulliClosed__DeepState_TestHarness_generation.cpp and BernoulliClosed__DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double BernoulliClosed_(double alpha_1, double beta_1, double alpha_2, double beta_2);

TEST(bayesAB_deepstate_test,BernoulliClosed__test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector alpha_1(1);
  alpha_1[0]  = RcppDeepState_double();
  qs::c_qsave(alpha_1,"/home/akhila/fuzzer_packages/fuzzedpackages/bayesAB/inst/testfiles/BernoulliClosed_/inputs/alpha_1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "alpha_1 values: "<< alpha_1 << std::endl;
  NumericVector beta_1(1);
  beta_1[0]  = RcppDeepState_double();
  qs::c_qsave(beta_1,"/home/akhila/fuzzer_packages/fuzzedpackages/bayesAB/inst/testfiles/BernoulliClosed_/inputs/beta_1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "beta_1 values: "<< beta_1 << std::endl;
  NumericVector alpha_2(1);
  alpha_2[0]  = RcppDeepState_double();
  qs::c_qsave(alpha_2,"/home/akhila/fuzzer_packages/fuzzedpackages/bayesAB/inst/testfiles/BernoulliClosed_/inputs/alpha_2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "alpha_2 values: "<< alpha_2 << std::endl;
  NumericVector beta_2(1);
  beta_2[0]  = RcppDeepState_double();
  qs::c_qsave(beta_2,"/home/akhila/fuzzer_packages/fuzzedpackages/bayesAB/inst/testfiles/BernoulliClosed_/inputs/beta_2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "beta_2 values: "<< beta_2 << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    BernoulliClosed_(alpha_1[0],beta_1[0],alpha_2[0],beta_2[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
