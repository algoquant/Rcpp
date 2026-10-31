// The z-scores are calculated using an exponential moving maximum 
// EMM volatility and the exponential moving average EMA prices.
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/run_zscores.cpp")

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

// using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;


// Calculated the z-scores using an exponential moving maximum 
// EMM volatility and the exponential moving average EMA prices.
// [[Rcpp::export]]
arma::mat run_zscorex(const arma::mat& timeser, double lambdaf) {
  
  arma::uword nrows = timeser.n_rows;
  // arma::uword ncols = timeser.n_cols;
  arma::mat meanm = arma::zeros(nrows, 1);
  // arma::mat vars = arma::ones(nrows, 1);
  arma::mat vols = arma::ones(nrows, 1);
  arma::mat zscores = arma::zeros(nrows, 1);
  double lambda1 = 1-lambdaf;
  double difx = 0;
  
  // Perform loop over the rows
  meanm(0) = timeser(0);
  vols(0) = abs(timeser(1) - timeser(0));
  // vars(0) = pow(vols(0), 2);
  // vars(0) = 0;
  for (arma::uword it = 1; it < nrows; it++) {
    // Variance is the weighted sum of the past variance and the square of the data minus its mean
    // vars(it) = lambdaf*vars(it-1) + lambda1*pow(timeser(it) - meanm(it-1), 2);
    // vols(it) = sqrt(vars(it));
    // Calculate the means using the decay factor
    meanm(it) = lambdaf*meanm(it-1) + lambda1*timeser(it);
    // Calculate the max using the decay factor
    difx = abs(timeser(it) - timeser(it-1));
    vols(it) = max(difx, lambdaf*vols(it-1) + lambda1*difx);
    // Calculate the z-scores
    zscores(it) = (timeser(it) - meanm(it-1))/vols(it);
  }  // end for
  
  // arma::mat zscores = (timeser - meanm)/vars;
  
  return arma::join_rows(meanm, vols, zscores);
  
}  // end run_zscorex


