// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_rcpp11.cpp")

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector foo(int nn){
  Rcpp::NumericVector out_put(nn);
  // out_put = {1, 2, 3};
  // out_put[1] = 1;
  // These two lines produce C++ warnings but they compile fine anyway
  for (int i = 0; i < out_put.size(); i++){
    out_put[i] = 2*i;
  }
  // Rcpp::NumericVector out_put = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
  return out_put;
}
