// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// The function lagvec_rcpp() lags a vector using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector lagvec_rcpp(NumericVector& vec_tor, int lagg=1) {
  // return vec_tor[Range(0, vec_tor.size()-1-lagg)];
  int len_gth = vec_tor.size();
  // double firs_t = vec_tor[0];
  Rcpp::NumericVector lagg_ed(len_gth);
  
  // The below copy operation doesn't copy anything (vector remains full of zeros):
  // lagg_ed[Range(lagg, len_gth-1)] = vec_tor[Range(0, len_gth-1-lagg)];
  
  // So instead copy using an explicit loop:
  for (int it = 0; it < (len_gth-lagg); it++) {
    lagg_ed[it+lagg] = vec_tor[it];
  }  // end for
  
  // This copy operation does work:
  lagg_ed[Range(0, lagg-1)] = rep(vec_tor[0], lagg);
  
  return lagg_ed;
}  // end lagvec_rcpp


// The function lagmat_rcpp() lags a vector using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix lagmat_rcpp(NumericMatrix& mat_rix, int lagg=1) {
  // return mat_rix[Range(0, mat_rix.nrow()-1-lagg)];
  int n_rows = mat_rix.nrow();
  int n_cols = mat_rix.ncol();
  // double firs_t = mat_rix[0];
  Rcpp::NumericMatrix lagg_ed(n_rows, n_cols);
  
  // The below copy operation doesn't copy anything (vector remains full of zeros):
  // lagg_ed[Range(lagg, n_rows-1)] = mat_rix[Range(0, n_rows-1-lagg)];
  
  // So instead copy using an explicit loop:
  for (int it = 0; it < (n_rows-lagg); it++) {
    lagg_ed(it+lagg, _) = mat_rix(it, _);
  }  // end for
  
  // Pad:
  for (int it = 0; it < lagg; it++) {
    lagg_ed(it, _) = mat_rix(0, _);
  }  // end for

  return lagg_ed;
}  // end lagmat_rcpp



// The function lagvec_arma() lags a vector using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::vec lagvec_arma(arma::vec& vec_tor, int lagg=1) {
  return arma::join_cols(arma::repelem(vec_tor.subvec(0, 0), lagg, 1), vec_tor.subvec(0, vec_tor.n_elem-1-lagg));
}  // end lagvec_arma


// The function lagmat_arma() lags a vector using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::mat lagmat_arma(arma::mat& mat_rix, int lagg=1) {
  return arma::join_cols(arma::repelem(mat_rix.row(0), lagg, 1), 
                         mat_rix.rows(0, mat_rix.n_rows-1-lagg));
}  // end lagmat_arma



//' @export
// [[Rcpp::export]]
NumericVector lag_ohlc(arma::mat& oh_lc, int lagg=1) {
  arma::colvec clo_se = oh_lc.col(3);
  Rcpp::NumericVector lag_close = Rcpp::NumericVector(clo_se.begin(), clo_se.end());
  // lag_close = lagvec_rcpp(lag_close, lagg=lagg);
  return lagvec_rcpp(lag_close, lagg=lagg);
}  // end lag_ohlc


// The function copy_loop() copies a vector using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector copy_loop(NumericVector& vec_tor, int lagg=1) {
  // return vec_tor[Range(0, vec_tor.size()-1-lagg)];
  // int len_gth = vec_tor.size();
  double firs_t = vec_tor[0];
  Rcpp::NumericVector lagg_ed(lagg);
  
  // Copy operation using an explicit loop:
  for (int it = 0; it < (lagg); it++) {
    lagg_ed[it] = firs_t;
  }  // end for
  
  return lagg_ed;
}  // end copy_loop


// The function copy_range() copies a vector using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector copy_range(NumericVector& vec_tor, int lagg=1) {
  // return vec_tor[Range(0, vec_tor.size()-1-lagg)];
  // int len_gth = vec_tor.size();
  double firs_t = vec_tor[0];
  Rcpp::NumericVector lagg_ed(lagg);
  
  // Copy operation using range:
  lagg_ed[Range(0, lagg-1)] = rep(firs_t, lagg);
  
  return lagg_ed;
}  // end copy_range


