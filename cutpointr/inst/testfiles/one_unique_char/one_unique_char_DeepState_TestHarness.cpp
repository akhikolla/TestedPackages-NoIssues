// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// one_unique_char_DeepState_TestHarness_generation.cpp and one_unique_char_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

LogicalVector one_unique_char(CharacterVector x);

TEST(cutpointr_deepstate_test,one_unique_char_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector x  = RcppDeepState_CharacterVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/cutpointr/inst/testfiles/one_unique_char/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    one_unique_char(x);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
