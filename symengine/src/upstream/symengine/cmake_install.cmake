# Install script for directory: /home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/libsymengine.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES
    "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symengine_config.h"
    "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symengine_config_cling.h"
    "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symengine_export.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/add.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/basic.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/basic-inl.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/basic-methods.inc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/complex_double.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/complex.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/complex_mpc.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/constants.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/cwrapper.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/derivative.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/dict.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/diophantine.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/eval_arb.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/eval_double.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/eval.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/eval_mpc.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/eval_mpfr.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/expression.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/fields.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/finitediff.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/flint_wrapper.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/functions.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/infinity.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/integer.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/lambda_double.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/llvm_double.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/logic.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/matrix.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/monomials.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/mp_class.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/mp_wrapper.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/mul.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/nan.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/ntheory.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/number.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/parser.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/parser" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/parser/parser.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/parser" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/parser/parser_stype.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/parser" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/parser/tokenizer.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/basic_conversions.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/cancel.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/uexprpoly.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/uintpoly_flint.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/uintpoly.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/uintpoly_piranha.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/upolybase.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/uratpoly.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/usymenginepoly.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/polys" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/polys/msymenginepoly.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/pow.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/printers" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/printers/codegen.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/printers" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/printers/mathml.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/printers" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/printers/strprinter.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine/printers" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/printers/latex.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/printers.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/rational.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/real_double.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/real_mpfr.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/rings.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/series_flint.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/series_generic.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/series.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/series_piranha.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/series_visitor.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/sets.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/solve.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/subs.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symbol.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symengine_assert.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symengine_casts.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symengine_exception.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/symengine_rcp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/type_codes.inc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/symengine" TYPE FILE FILES "/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/visitor.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/akhila/fuzzer_packages/fuzzedpackages/symengine/src/upstream/symengine/utilities/matchpycpp/cmake_install.cmake")

endif()

