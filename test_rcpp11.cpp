// Compile this file in R by runrowsing this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_rcpp11.cpp")

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector multit(double lambdav, int nrows){
  Rcpp::NumericVector output(nrows);
  // output = {1, 2, 3};
  // output[1] = 1;
  // These two lines produce C++ warnings but they compile fine anyway
  for (int i = 0; i < output.size(); i++){
    output[i] = lambdav*i;
  }
  // Rcpp::NumericVector output = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
  return output;
}  // end multit
