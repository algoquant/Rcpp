////////////////////////////
// Matrix functions using Armadillo and the Standard Template Library (STL)
////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// Use STL
// using namespace std;
// #include <functional>
// #include <string>
// #include <map>
// #include <iostream>
// [[Rcpp::plugins(cpp11)]]

arma::mat mult_mat(arma::vec vector, arma::mat matrix, bool byrow = true);

void mult_mat_ref(arma::vec vector, arma::mat matrix, bool byrow = true);

Rcpp::List calc_eigen(const arma::mat& tseries);

arma::mat calc_inv(const arma::mat& tseries, double eigen_thresh = 0.01, arma::uword dimax = 0);

arma::mat calc_scaled(const arma::mat& tseries, bool use_median = false);



