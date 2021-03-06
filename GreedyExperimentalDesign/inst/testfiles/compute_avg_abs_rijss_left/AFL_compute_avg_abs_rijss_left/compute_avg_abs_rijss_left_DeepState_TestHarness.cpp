// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// compute_avg_abs_rijss_left_DeepState_TestHarness_generation.cpp and compute_avg_abs_rijss_left_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <ctime>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector compute_avg_abs_rijss_left(NumericMatrix Rij);

TEST(GreedyExperimentalDesign_deepstate_test,compute_avg_abs_rijss_left_test){
  RInside R;
  std::time_t t = std::time(0);
  std::cout << "input starts" << std::endl;
  NumericMatrix Rij  = RcppDeepState_NumericMatrix();
  std::string Rij_t = "/home/akhila/fuzzer_packages/fuzzedpackages/GreedyExperimentalDesign/inst/testfiles/compute_avg_abs_rijss_left/AFL_compute_avg_abs_rijss_left/afl_inputs/" + std::to_string(t) + "_Rij.qs";
  qs::c_qsave(Rij,Rij_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "Rij values: "<< Rij << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    compute_avg_abs_rijss_left(Rij);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
