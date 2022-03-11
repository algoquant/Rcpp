// Rcpp header with information for C++ compiler
#include <RcppArmadillo.h> // include C++ header files
using namespace arma; // use Armadillo C++ namespace
// declare dependency on RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

// export function roll_maxmin() to R
// [[Rcpp::export]]
arma::mat roll_maxmin(const arma::vec& vectorv, 
                      const arma::uword& look_back) {
  arma::uword nrows = vectorv.size();
  arma::mat max_min(nrows, 2);
  arma::vec sub_vec;
  // startup period
  max_min(0, 0) = vectorv[0];
  max_min(0, 1) = vectorv[0];
  for (uword it = 1; it < look_back; it++) {
    sub_vec = vectorv.subvec(0, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  // remaining periods
  for (uword it = look_back; it < nrows; it++) {
    sub_vec = vectorv.subvec(it-look_back+1, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  return max_min;
}  // end roll_maxmin
