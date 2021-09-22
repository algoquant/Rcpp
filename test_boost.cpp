////////////////////////////
// Functions to test syntax of the Boost library
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp("C:/Develop/R/Rcpp/test_boost.cpp")

#include <Rcpp.h>
// include Boost headers from R package BH  
#include <boost/math/common_factor.hpp>
// [[Rcpp::depends(BH)]]    

using namespace Rcpp;

// [[Rcpp::export]]   
int computeGCD(int a, int b) {
  return boost::math::gcd(a, b);
}  // end computeGCD

