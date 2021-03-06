// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// q10_mrme_approx_DeepState_TestHarness_generation.cpp and q10_mrme_approx_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double q10_mrme_approx(NumericVector z, double t, NumericVector theta, NumericVector integrControl, NumericVector err_start, NumericVector err_end, NumericVector err_end_prob);

TEST(smam_deepstate_test,q10_mrme_approx_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector z  = RcppDeepState_NumericVector();
  qs::c_qsave(z,"/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/q10_mrme_approx/inputs/z.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "z values: "<< z << std::endl;
  NumericVector t(1);
  t[0]  = RcppDeepState_double();
  qs::c_qsave(t,"/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/q10_mrme_approx/inputs/t.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "t values: "<< t << std::endl;
  NumericVector theta  = RcppDeepState_NumericVector();
  qs::c_qsave(theta,"/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/q10_mrme_approx/inputs/theta.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "theta values: "<< theta << std::endl;
  NumericVector integrControl  = RcppDeepState_NumericVector();
  qs::c_qsave(integrControl,"/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/q10_mrme_approx/inputs/integrControl.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "integrControl values: "<< integrControl << std::endl;
  NumericVector err_start  = RcppDeepState_NumericVector();
  qs::c_qsave(err_start,"/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/q10_mrme_approx/inputs/err_start.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "err_start values: "<< err_start << std::endl;
  NumericVector err_end  = RcppDeepState_NumericVector();
  qs::c_qsave(err_end,"/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/q10_mrme_approx/inputs/err_end.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "err_end values: "<< err_end << std::endl;
  NumericVector err_end_prob  = RcppDeepState_NumericVector();
  qs::c_qsave(err_end_prob,"/home/akhila/fuzzer_packages/fuzzedpackages/smam/inst/testfiles/q10_mrme_approx/inputs/err_end_prob.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "err_end_prob values: "<< err_end_prob << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    q10_mrme_approx(z,t[0],theta,integrControl,err_start,err_end,err_end_prob);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
