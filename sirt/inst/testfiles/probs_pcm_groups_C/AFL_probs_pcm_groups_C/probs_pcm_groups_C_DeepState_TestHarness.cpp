// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// probs_pcm_groups_C_DeepState_TestHarness_generation.cpp and probs_pcm_groups_C_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <ctime>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::List probs_pcm_groups_C(Rcpp::NumericMatrix dat, Rcpp::NumericMatrix dat_resp, Rcpp::NumericVector group, Rcpp::NumericMatrix probs, int CC, int TP);

TEST(sirt_deepstate_test,probs_pcm_groups_C_test){
  RInside R;
  std::time_t t = std::time(0);
  std::cout << "input starts" << std::endl;
  NumericMatrix dat  = RcppDeepState_NumericMatrix();
  std::string dat_t = "/home/akhila/fuzzer_packages/fuzzedpackages/sirt/inst/testfiles/probs_pcm_groups_C/AFL_probs_pcm_groups_C/afl_inputs/" + std::to_string(t) + "_dat.qs";
  qs::c_qsave(dat,dat_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "dat values: "<< dat << std::endl;
  NumericMatrix dat_resp  = RcppDeepState_NumericMatrix();
  std::string dat_t = "/home/akhila/fuzzer_packages/fuzzedpackages/sirt/inst/testfiles/probs_pcm_groups_C/AFL_probs_pcm_groups_C/afl_inputs/" + std::to_string(t) + "_dat.qs";
  std::string dat_resp_t = "/home/akhila/fuzzer_packages/fuzzedpackages/sirt/inst/testfiles/probs_pcm_groups_C/AFL_probs_pcm_groups_C/afl_inputs/" + std::to_string(t) + "_dat_resp.qs";
  qs::c_qsave(dat_resp,dat_resp_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "dat_resp values: "<< dat_resp << std::endl;
  NumericVector group  = RcppDeepState_NumericVector();
  std::string group_t = "/home/akhila/fuzzer_packages/fuzzedpackages/sirt/inst/testfiles/probs_pcm_groups_C/AFL_probs_pcm_groups_C/afl_inputs/" + std::to_string(t) + "_group.qs";
  qs::c_qsave(group,group_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "group values: "<< group << std::endl;
  NumericMatrix probs  = RcppDeepState_NumericMatrix();
  std::string probs_t = "/home/akhila/fuzzer_packages/fuzzedpackages/sirt/inst/testfiles/probs_pcm_groups_C/AFL_probs_pcm_groups_C/afl_inputs/" + std::to_string(t) + "_probs.qs";
  qs::c_qsave(probs,probs_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "probs values: "<< probs << std::endl;
  IntegerVector CC(1);
  CC[0]  = RcppDeepState_int();
  std::string CC_t = "/home/akhila/fuzzer_packages/fuzzedpackages/sirt/inst/testfiles/probs_pcm_groups_C/AFL_probs_pcm_groups_C/afl_inputs/" + std::to_string(t) + "_CC.qs";
  qs::c_qsave(CC,CC_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "CC values: "<< CC << std::endl;
  IntegerVector TP(1);
  TP[0]  = RcppDeepState_int();
  std::string TP_t = "/home/akhila/fuzzer_packages/fuzzedpackages/sirt/inst/testfiles/probs_pcm_groups_C/AFL_probs_pcm_groups_C/afl_inputs/" + std::to_string(t) + "_TP.qs";
  qs::c_qsave(TP,TP_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "TP values: "<< TP << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    probs_pcm_groups_C(dat,dat_resp,group,probs,CC[0],TP[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
