////////////////////////////
// Functions to test C++ syntax with Armadillo
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_temp.cpp")

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::uword calc_pos(arma::vec datav, arma::uword eigen_max, double eigen_thresh = 0.001) {
  
  return min(eigen_max, arma::sum(datav > eigen_thresh) - 1);
  
}  // end calc_pos


// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword eigen_max = 0) {
  
  // Allocate SVD variables
  arma::vec svd_val;  // Singular values
  arma::mat svd_u, svd_v;  // Singular matrices
  // Calculate the SVD
  arma::svd(svd_u, svd_val, svd_v, tseries);
  // Calculate the number of non-small singular values
  arma::uword svd_num = arma::sum(svd_val > eigen_thresh*arma::sum(svd_val));
  
  if (eigen_max == 0) {
    // Set eigen_max
    eigen_max = svd_num - 1;
  } else {
    // Adjust eigen_max
    eigen_max = min(eigen_max - 1, svd_num - 1);
  }  // end if
  
  // Remove all small singular values
  svd_val = svd_val.subvec(0, eigen_max);
  svd_u = svd_u.cols(0, eigen_max);
  svd_v = svd_v.cols(0, eigen_max);
  
  // Calculate the regularized inverse from the SVD decomposition
  return svd_v*arma::diagmat(1/svd_val)*svd_u.t();
  
}  // end calc_inv


