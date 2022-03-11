// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// The function roll_sum() calculates the rolling sum over a vector using Rcpp
//' @export
// [[Rcpp::export]]
NumericVector roll_sum(NumericVector vectorv, int look_back) {
  int nrows = vectorv.size();
  NumericVector rolling_sum(nrows);
  
  rolling_sum[0] = vectorv[0];
  for (int it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vectorv[it];
  }  // end for

  for (int it = look_back; it < nrows; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vectorv[it] - vectorv[it-look_back];
  }  // end for
  
  return rolling_sum;
}  // end roll_sum


// The function roll_sum_arma() calculates the rolling sum over a vector using 
// RcppArmadillo
//' @export
// [[Rcpp::export]]
arma::vec roll_sum_arma(arma::vec& vectorv, arma::uword look_back=11) {
  arma::vec lag_ged = arma::zeros(vectorv.n_elem);
  // arma::uword nrows = vectorv.n_elem;
  // arma::vec rolling_sum(nrows);
  
  // Warmup period
  // rolling_sum[0] = vectorv[0];
  // for (arma::uword it = 1; it < look_back; it++) {
  //   rolling_sum[it] = rolling_sum[it-1] + vectorv[it];
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = look_back; it < nrows; it++) {
  //   rolling_sum[it] = rolling_sum[it-1] + vectorv[it] - vectorv[it-look_back];
  // }  // end for
  
  // return rolling_sum;
  
  vectorv = arma::cumsum(vectorv);
  // lag_ged(1, look_back) = 0;
  lag_ged.subvec(look_back, lag_ged.n_elem-1) = vectorv.subvec(0, vectorv.n_elem-look_back-1);
  
  return (vectorv - lag_ged);
  // return vectorv;
}  // end roll_sum_arma


