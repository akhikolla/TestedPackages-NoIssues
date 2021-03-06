// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// makeprop_DeepState_TestHarness_generation.cpp and makeprop_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector makeprop(NumericVector mu, NumericVector sample, NumericVector low, NumericVector up);

TEST(yuima_deepstate_test,makeprop_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector mu  = RcppDeepState_NumericVector();
  qs::c_qsave(mu,"/home/akhila/fuzzer_packages/fuzzedpackages/yuima/inst/testfiles/makeprop/inputs/mu.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "mu values: "<< mu << std::endl;
  NumericVector sample  = RcppDeepState_NumericVector();
  qs::c_qsave(sample,"/home/akhila/fuzzer_packages/fuzzedpackages/yuima/inst/testfiles/makeprop/inputs/sample.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "sample values: "<< sample << std::endl;
  NumericVector low  = RcppDeepState_NumericVector();
  qs::c_qsave(low,"/home/akhila/fuzzer_packages/fuzzedpackages/yuima/inst/testfiles/makeprop/inputs/low.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "low values: "<< low << std::endl;
  NumericVector up  = RcppDeepState_NumericVector();
  qs::c_qsave(up,"/home/akhila/fuzzer_packages/fuzzedpackages/yuima/inst/testfiles/makeprop/inputs/up.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "up values: "<< up << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    makeprop(mu,sample,low,up);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
