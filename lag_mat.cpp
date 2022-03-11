// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// The function lagvec_rcpp() lags a vector using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector lagvec_rcpp(NumericVector& vectorv, int lagg=1) {
  // return vectorv[Range(0, vectorv.size()-1-lagg)];
  int nrows = vectorv.size();
  // double firs_t = vectorv[0];
  Rcpp::NumericVector lagg_ed(nrows);
  
  // The below copy operation doesn't copy anything (vector remains full of zeros):
  // lagg_ed[Range(lagg, nrows-1)] = vectorv[Range(0, nrows-1-lagg)];
  
  // So instead copy using an explicit loop:
  for (int it = 0; it < (nrows-lagg); it++) {
    lagg_ed[it+lagg] = vectorv[it];
  }  // end for
  
  // This copy operation does work:
  lagg_ed[Range(0, lagg-1)] = rep(vectorv[0], lagg);
  
  return lagg_ed;
}  // end lagvec_rcpp


// The function lagmat_rcpp() lags a vector using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix lagmat_rcpp(NumericMatrix& matrixv, int lagg=1) {
  // return matrixv[Range(0, matrixv.nrow()-1-lagg)];
  int nrows = matrixv.nrow();
  int ncols = matrixv.ncol();
  // double firs_t = matrixv[0];
  Rcpp::NumericMatrix lagg_ed(nrows, ncols);
  
  // The below copy operation doesn't copy anything (vector remains full of zeros):
  // lagg_ed[Range(lagg, nrows-1)] = matrixv[Range(0, nrows-1-lagg)];
  
  // So instead copy using an explicit loop:
  for (int it = 0; it < (nrows-lagg); it++) {
    lagg_ed(it+lagg, _) = matrixv(it, _);
  }  // end for
  
  // Pad:
  for (int it = 0; it < lagg; it++) {
    lagg_ed(it, _) = matrixv(0, _);
  }  // end for

  return lagg_ed;
}  // end lagmat_rcpp



// The function lagvec_arma() lags a vector using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::vec lagvec_arma(arma::vec& vectorv, int lagg=1) {
  return arma::join_cols(arma::repelem(vectorv.subvec(0, 0), lagg, 1), vectorv.subvec(0, vectorv.n_elem-1-lagg));
}  // end lagvec_arma


// The function lagmat_arma() lags a vector using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::mat lagmat_arma(arma::mat& matrixv, int lagg=1) {
  return arma::join_cols(arma::repelem(matrixv.row(0), lagg, 1), 
                         matrixv.rows(0, matrixv.n_rows-1-lagg));
}  // end lagmat_arma



//' @export
// [[Rcpp::export]]
NumericVector lag_ohlc(arma::mat& ohlc, int lagg=1) {
  arma::colvec clo_se = ohlc.col(3);
  Rcpp::NumericVector lag_close = Rcpp::NumericVector(clo_se.begin(), clo_se.end());
  // lag_close = lagvec_rcpp(lag_close, lagg=lagg);
  return lagvec_rcpp(lag_close, lagg=lagg);
}  // end lag_ohlc


// The function copy_loop() copies a vector using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector copy_loop(NumericVector& vectorv, int lagg=1) {
  // return vectorv[Range(0, vectorv.size()-1-lagg)];
  // int nrows = vectorv.size();
  double firs_t = vectorv[0];
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
Rcpp::NumericVector copy_range(NumericVector& vectorv, int lagg=1) {
  // return vectorv[Range(0, vectorv.size()-1-lagg)];
  // int nrows = vectorv.size();
  double firs_t = vectorv[0];
  Rcpp::NumericVector lagg_ed(lagg);
  
  // Copy operation using range:
  lagg_ed[Range(0, lagg-1)] = rep(firs_t, lagg);
  
  return lagg_ed;
}  // end copy_range


