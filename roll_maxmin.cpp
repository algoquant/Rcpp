// Rcpp header with information for C++ compiler
#include <RcppArmadillo.h> // include C++ header files
using namespace arma; // use Armadillo C++ namespace
// declare dependency on RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

// export function roll_maxmin() to R
// [[Rcpp::export]]
arma::mat roll_maxmin(const arma::vec& vec_tor, 
                      const arma::uword& look_back) {
  arma::uword num_rows = vec_tor.size();
  arma::mat max_min(num_rows, 2);
  arma::vec sub_vec;
  // startup period
  max_min(0, 0) = vec_tor[0];
  max_min(0, 1) = vec_tor[0];
  for (uword it = 1; it < look_back; it++) {
    sub_vec = vec_tor.subvec(0, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  // remaining periods
  for (uword it = look_back; it < num_rows; it++) {
    sub_vec = vec_tor.subvec(it-look_back+1, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  return max_min;
}  // end roll_maxmin
