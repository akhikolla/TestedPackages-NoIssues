// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// apoplasticWaterPotential_DeepState_TestHarness_generation.cpp and apoplasticWaterPotential_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double apoplasticWaterPotential(double RWC, double c, double d);

TEST(medfate_deepstate_test,apoplasticWaterPotential_test){
  static int rinside_flag = 0;
  if(rinside_flag == 0)
  {
    rinside_flag = 1;
    RInside R;
  } std::time_t current_timestamp = std::time(0);
  std::cout << "input starts" << std::endl;
  NumericVector RWC(1);
  RWC[0]  = RcppDeepState_double();
  std::string RWC_t = "/home/akhila/fuzzer_packages/fuzzedpackages/medfate/inst/testfiles/apoplasticWaterPotential/libFuzzer_apoplasticWaterPotential/libfuzzer_inputs/" + std::to_string(current_timestamp) +
          "_RWC.qs";
  qs::c_qsave(RWC,RWC_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "RWC values: "<< RWC << std::endl;
  NumericVector c(1);
  c[0]  = RcppDeepState_double();
  std::string c_t = "/home/akhila/fuzzer_packages/fuzzedpackages/medfate/inst/testfiles/apoplasticWaterPotential/libFuzzer_apoplasticWaterPotential/libfuzzer_inputs/" + std::to_string(current_timestamp) +
          "_c.qs";
  qs::c_qsave(c,c_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "c values: "<< c << std::endl;
  NumericVector d(1);
  d[0]  = RcppDeepState_double();
  std::string d_t = "/home/akhila/fuzzer_packages/fuzzedpackages/medfate/inst/testfiles/apoplasticWaterPotential/libFuzzer_apoplasticWaterPotential/libfuzzer_inputs/" + std::to_string(current_timestamp) +
          "_d.qs";
  qs::c_qsave(d,d_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "d values: "<< d << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    apoplasticWaterPotential(RWC[0],c[0],d[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
