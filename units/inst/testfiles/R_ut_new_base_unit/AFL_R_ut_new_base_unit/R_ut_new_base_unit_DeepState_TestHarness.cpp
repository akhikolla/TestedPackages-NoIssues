// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// R_ut_new_base_unit_DeepState_TestHarness_generation.cpp and R_ut_new_base_unit_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <ctime>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

void R_ut_new_base_unit(CharacterVector name);

TEST(units_deepstate_test,R_ut_new_base_unit_test){
  RInside R;
  std::time_t t = std::time(0);
  std::cout << "input starts" << std::endl;
  CharacterVector name  = RcppDeepState_CharacterVector();
  std::string name_t = "/home/akhila/fuzzer_packages/fuzzedpackages/units/inst/testfiles/R_ut_new_base_unit/AFL_R_ut_new_base_unit/afl_inputs/" + std::to_string(t) + "_name.qs";
  qs::c_qsave(name,name_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "name values: "<< name << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    R_ut_new_base_unit(name);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
