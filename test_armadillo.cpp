////////////////////////////
// Functions to test C++ syntax with Armadillo
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_armadillo.cpp")

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


////////////////////////////////////////////////
// Tests misc

//' @export
// [[Rcpp::export]]
arma::mat floor_it(arma::mat& data, double minval) {
  
  arma::mat mydata = data;
  mydata.transform([&minval](double x) {return max(x, minval);});
  return mydata;
  
}  // end floor_it





////////////////////////////////////////////////
// Included to facilitate Tests - remove later - don't migrate


////////////////////////////////////////////////
// Test versions to be migrated to package HighFreq




////////////////////////////////////////////////
// Old stuff - can be deleted later

// Calculate the rolling maximum or minimum of streaming data using a lambda decay factor
//' @export
// [[Rcpp::export]]
arma::mat run_maxmin(arma::mat tseries, double lambda, bool calc_max = true) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat maxmin = tseries;
  arma::mat means = tseries;
  double lambda1 = 1-lambda;
  
  if (calc_max) {
    // Perform loop over rows
    for (arma::uword it = 1; it < num_rows; it++) {
      // Calculate the mean as a weighted sum
      means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
      // Calculate the max from a weighted sum
      maxmin.row(it) = arma::max(maxmin.row(it), means.row(it-1) + lambda*(maxmin.row(it-1) - means.row(it-1)));
    }  // end for
  } else {
    // Perform loop over rows
    for (arma::uword it = 1; it < num_rows; it++) {
      // Calculate the mean as a weighted sum
      means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
      // Calculate the max from a weighted sum
      maxmin.row(it) = arma::min(maxmin.row(it), means.row(it-1) + lambda*(maxmin.row(it-1) - means.row(it-1)));
    }  // end for
  }  // end if
  
  return maxmin;
  
}  // end run_maxmin


// Calculate the rolling maximum of streaming data using a lambda decay factor
//' @export
// [[Rcpp::export]]
arma::colvec armax(arma::colvec tseries, arma::colvec tseries2) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::colvec maxs = tseries;
  
  // Perform loop over rows
  for (arma::uword it = 0; it < num_rows; it++) {
    // Calculate the mean as a weighted sum
    // means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
    // Calculate the max from a weighted sum
    // maxs.row(it) = arma::max(maxs.row(it), means.row(it-1) + lambda*(maxs.row(it-1) - means.row(it-1)));
    maxs.row(it) = arma::max(tseries.row(it), tseries2.row(it));
  }  // end for
  
  return maxs;
  
}  // end armax

