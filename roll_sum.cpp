// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// The function roll_sum() calculates the rolling sum over a vector using Rcpp
//' @export
// [[Rcpp::export]]
NumericVector roll_sum(NumericVector vec_tor, int look_back) {
  int len_gth = vec_tor.size();
  NumericVector rolling_sum(len_gth);
  
  rolling_sum[0] = vec_tor[0];
  for (int it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it];
  }  // end for

  for (int it = look_back; it < len_gth; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it] - vec_tor[it-look_back];
  }  // end for
  
  return rolling_sum;
}  // end roll_sum


// The function roll_sum_arma() calculates the rolling sum over a vector using 
// RcppArmadillo
//' @export
// [[Rcpp::export]]
arma::vec roll_sum_arma(arma::vec& vec_tor, arma::uword look_back=11) {
  arma::vec lag_ged = arma::zeros(vec_tor.n_elem);
  // arma::uword len_gth = vec_tor.n_elem;
  // arma::vec rolling_sum(len_gth);
  
  // Warmup period
  // rolling_sum[0] = vec_tor[0];
  // for (arma::uword it = 1; it < look_back; it++) {
  //   rolling_sum[it] = rolling_sum[it-1] + vec_tor[it];
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = look_back; it < len_gth; it++) {
  //   rolling_sum[it] = rolling_sum[it-1] + vec_tor[it] - vec_tor[it-look_back];
  // }  // end for
  
  // return rolling_sum;
  
  vec_tor = arma::cumsum(vec_tor);
  // lag_ged(1, look_back) = 0;
  lag_ged.subvec(look_back, lag_ged.n_elem-1) = vec_tor.subvec(0, vec_tor.n_elem-look_back-1);
  
  return (vec_tor - lag_ged);
  // return vec_tor;
}  // end roll_sum_arma


