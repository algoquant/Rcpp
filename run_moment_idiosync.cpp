// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

////////////////////////////
// The function run_moment_idiosync() simulates an EMA momentum strategy 
// using the idiosyncratic returns. 
// The momentum weights are equal to the EMA idiosyncratic returns divided 
// by their trailing EMA variance.
// 
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/run_moment_idiosync.cpp")
// 
// Copyright: 2026 Jerzy Pawlowski
////////////////////////////

#include "RcppArmadillo.h"
#include <vector>
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
// Use STL
using namespace std;
// For eigen solver SymEigsSolver
using namespace arma::newarp;


// Prompts:
// Modify the function run_reg() above to create a new function.
// Create the function called run_moment_idiosync() to calculate the EMA idiosyncratic returns divided by their trailing EMA variance of a time series.
// run_moment_idiosync() should accept multiple columns of returns data.
// The first column is the market returns, and the remaining columns are the individual asset returns.
// In run_moment_idiosync() change the order of the loops: first loop over rows, second loop over stocks.
// Don't use intercept - the market is only single predictor


////////////////////////////////////////////////////////////
//' Simulate an EMA momentum strategy using the idiosyncratic returns.
//'
//' @param \code{returns} A \emph{matrix} of returns data where the first column
//'   is the market returns and the remaining columns are individual asset
//'   returns.
//'   
//' @param \code{lambdaf} A decay factor which multiplies past estimates.
//'   
//' @return A \emph{matrix} with the same number of rows as the input
//'   \code{returns}, and number of columns equal to the number of individual
//'   stocks (excluding the market column). Each column contains the
//'   standardized idiosyncratic returns for each asset.
//'
//' @details
//'   The function \code{run_moment_idiosync()} calculates idiosyncratic returns
//'   by regressing each individual asset return against the market return:
//'   \deqn{
//'     r_{i,t} = \beta_i r_{m,t} + \epsilon_{i,t}
//'   }
//'   Where \eqn{r_{i,t}} are the individual asset returns, \eqn{r_{m,t}} are
//'   the market returns, and \eqn{\epsilon_{i,t}} are the idiosyncratic returns.
//'   Note that no intercept term is included in the regression.
//'   
//'   It then calculates the trailing EMA variance of the idiosyncratic returns
//'   and standardizes them by dividing by their volatility:
//'   \deqn{
//'     z_{i,t} = \frac{\epsilon_{i,t}}{\sigma_{i,t}}
//'   }
//'
//' @export
// [[Rcpp::export]]
arma::mat run_moment_idiosync(const arma::mat& returns, // Returns matrix: first column = market, rest = individual stocks
                       double lambdaf) { // Decay factor

  arma::uword nrows = returns.n_rows;
  arma::uword ncols = returns.n_cols;
  // Number of individual stocks
  arma::uword nstocks = ncols - 1;
  
  if (ncols < 2) {
    Rcpp::stop("Returns matrix must have at least 2 columns (market + individual stocks)");
  }
  
  // Extract market returns (first column)
  arma::colvec retm = returns.col(0);
  
  // Initialize output matrix for the Kelly ratios of the idiosyncratic returns
  arma::mat pnlm = arma::zeros(nrows, nstocks);
  
  double lambda1 = 1 - lambdaf;
  // Decay factors squared for the variance calculations
  double lambda2 = pow(lambdaf, 2);
  double lambda21 = 1 - lambda2;
  
  // Initialize variables for all the stocks
  // Covariance of market with itself (scalar)
  double varm = retm(0) * retm(0); 
  // Covariances between the stocks and the market
  arma::rowvec covm = arma::zeros<arma::rowvec>(nstocks);
  // Beta coefficients for each stock
  arma::rowvec betav = arma::zeros<arma::rowvec>(nstocks); 
  // Idiosyncratic returns for each stock
  arma::rowvec retid = arma::zeros<arma::rowvec>(nstocks); 
  // Mean of the idiosyncratic returns for each stock
  arma::rowvec retidm = arma::zeros<arma::rowvec>(nstocks); 
  // Variance of idiosyncratic returns for each stock
  arma::rowvec varid = arma::zeros<arma::rowvec>(nstocks); 
  // The vector of weights for each stock
  arma::rowvec weightv = arma::zeros<arma::rowvec>(nstocks); 
  
  // Initialize for first observation for all the stocks
  // varm = retm(0) * retm(0);
  // Current stock return
  double rets = 0;
  // Loop over the stocks
  for (arma::uword coln = 0; coln < nstocks; coln++) {
    rets = returns(0, coln + 1);
    covm(coln) = rets * retm(0);
    betav(coln) = covm(coln) / varm;
    retid(coln) = rets - betav(coln) * retm(0);
    retidm(coln) = retid(coln);
    varid(coln) = retid(coln) * retid(coln);
  } // end for stocks
  
  // Loop over rows first
  for (arma::uword it = 1; it < (nrows-1); it++) {

    // Update market variance (same for all stocks)
    varm = lambda2 * varm + lambda21 * retm(it) * retm(it);
    
    // Loop over stocks second
    for (arma::uword coln = 0; coln < nstocks; coln++) {

      // Current stock return
      rets = returns(it, coln + 1); 
      // Calculate the idiosyncratic return (residual)
      retid(coln) = rets - betav(coln) * retm(it);
      // Calculate the momentum PnLs by multiplying the Kelly weights times the idiosyncratic returns
      pnlm(it, coln) = weightv(coln) * retid(coln);
      
      // Update trailing mean and variance of idiosyncratic returns
      retidm(coln) = lambdaf * retidm(coln) + lambda1 * retid(coln);
      double retidx = retid(coln) - retidm(coln);
      varid(coln) = lambda2 * varid(coln) + lambda21 * retidx * retidx;
      // Update the covariance
      covm(coln) = lambda2 * covm(coln) + lambda21 * rets * retm(it);
      // Update the beta
      betav(coln) = covm(coln) / varm;
      // Calculate the Kelly weights equal to the mean idiosyncratic returns divided by their variance
      weightv(coln) = retidm(coln) / varid(coln);
      
    } // end for stocks

  } // end for rows
  
  return pnlm;
  
}  // end run_moment_idiosync

