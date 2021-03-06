// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// qgamma_ssd_DeepState_TestHarness_generation.cpp and qgamma_ssd_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <ctime>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double qgamma_ssd(double p, double shape, double scale);

TEST(ssdtools_deepstate_test,qgamma_ssd_test){
  RInside R;
  std::time_t t = std::time(0);
  std::cout << "input starts" << std::endl;
  NumericVector p(1);
  p[0]  = RcppDeepState_double();
  std::string p_t = "/home/akhila/fuzzer_packages/fuzzedpackages/ssdtools/inst/testfiles/qgamma_ssd/AFL_qgamma_ssd/afl_inputs/" + std::to_string(t) + "_p.qs";
  qs::c_qsave(p,p_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "p values: "<< p << std::endl;
  NumericVector shape(1);
  shape[0]  = RcppDeepState_double();
  std::string shape_t = "/home/akhila/fuzzer_packages/fuzzedpackages/ssdtools/inst/testfiles/qgamma_ssd/AFL_qgamma_ssd/afl_inputs/" + std::to_string(t) + "_shape.qs";
  qs::c_qsave(shape,shape_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "shape values: "<< shape << std::endl;
  NumericVector scale(1);
  scale[0]  = RcppDeepState_double();
  std::string scale_t = "/home/akhila/fuzzer_packages/fuzzedpackages/ssdtools/inst/testfiles/qgamma_ssd/AFL_qgamma_ssd/afl_inputs/" + std::to_string(t) + "_scale.qs";
  qs::c_qsave(scale,scale_t,
		"high", "zstd", 1, 15, true, 1);
  std::cout << "scale values: "<< scale << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    qgamma_ssd(p[0],shape[0],scale[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
