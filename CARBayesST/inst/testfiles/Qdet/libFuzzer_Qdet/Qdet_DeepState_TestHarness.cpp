// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// Qdet_DeepState_TestHarness_generation.cpp and Qdet_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector Qdet(const int Nchains, NumericVector rho, NumericVector Wstar_val);

TEST(CARBayesST_deepstate_test,Qdet_test){
  static int rinside_flag = 0;
  if(rinside_flag == 0)
  {
    rinside_flag = 1;
    RInside R;
  } std::time_t current_timestamp = std::time(0);
  std::cout << "input starts" << std::endl;
  IntegerVector Nchains(1);
  Nchains[0]  = RcppDeepState_int();
  std::string Nchains_t = "/home/akhila/fuzzer_packages/fuzzedpackages/CARBayesST/inst/testfiles/Qdet/libFuzzer_Qdet/libfuzzer_inputs/" + std::to_string(current_timestamp) +
          "_Nchains.qs";
  qs::c_qsave(Nchains,Nchains_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "Nchains values: "<< Nchains << std::endl;
  NumericVector rho  = RcppDeepState_NumericVector();
  std::string rho_t = "/home/akhila/fuzzer_packages/fuzzedpackages/CARBayesST/inst/testfiles/Qdet/libFuzzer_Qdet/libfuzzer_inputs/" + std::to_string(current_timestamp) +
          "_rho.qs";
  qs::c_qsave(rho,rho_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "rho values: "<< rho << std::endl;
  NumericVector Wstar_val  = RcppDeepState_NumericVector();
  std::string Wstar_val_t = "/home/akhila/fuzzer_packages/fuzzedpackages/CARBayesST/inst/testfiles/Qdet/libFuzzer_Qdet/libfuzzer_inputs/" + std::to_string(current_timestamp) +
          "_Wstar_val.qs";
  qs::c_qsave(Wstar_val,Wstar_val_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "Wstar_val values: "<< Wstar_val << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    Qdet(Nchains[0],rho,Wstar_val);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
