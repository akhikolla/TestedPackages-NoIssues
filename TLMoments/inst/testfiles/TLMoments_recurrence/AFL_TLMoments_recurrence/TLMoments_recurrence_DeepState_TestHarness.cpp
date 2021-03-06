// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// TLMoments_recurrence_DeepState_TestHarness_generation.cpp and TLMoments_recurrence_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <ctime>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector TLMoments_recurrence(NumericVector x, int maxr, double s, double t);

TEST(TLMoments_deepstate_test,TLMoments_recurrence_test){
  RInside R;
  std::time_t t = std::time(0);
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  std::string x_t = "/home/akhila/fuzzer_packages/fuzzedpackages/TLMoments/inst/testfiles/TLMoments_recurrence/AFL_TLMoments_recurrence/afl_inputs/" + std::to_string(t) + "_x.qs";
  qs::c_qsave(x,x_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  IntegerVector maxr(1);
  maxr[0]  = RcppDeepState_int();
  std::string maxr_t = "/home/akhila/fuzzer_packages/fuzzedpackages/TLMoments/inst/testfiles/TLMoments_recurrence/AFL_TLMoments_recurrence/afl_inputs/" + std::to_string(t) + "_maxr.qs";
  qs::c_qsave(maxr,maxr_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "maxr values: "<< maxr << std::endl;
  NumericVector s(1);
  s[0]  = RcppDeepState_double();
  std::string s_t = "/home/akhila/fuzzer_packages/fuzzedpackages/TLMoments/inst/testfiles/TLMoments_recurrence/AFL_TLMoments_recurrence/afl_inputs/" + std::to_string(t) + "_s.qs";
  qs::c_qsave(s,s_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "s values: "<< s << std::endl;
  NumericVector t(1);
  t[0]  = RcppDeepState_double();
  std::string t_t = "/home/akhila/fuzzer_packages/fuzzedpackages/TLMoments/inst/testfiles/TLMoments_recurrence/AFL_TLMoments_recurrence/afl_inputs/" + std::to_string(t) + "_t.qs";
  qs::c_qsave(t,t_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "t values: "<< t << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    TLMoments_recurrence(x,maxr[0],s[0],t[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
