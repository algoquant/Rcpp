////////////////////////////
// Test the functions in HighFreq.cpp
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_HighFreq.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;

////////////////////////////////////////////////////////////
// Rcpp and RcppArmadillo functions for package HighFreq
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// Functions for matrix algebra
////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat run_rego(const arma::mat& response, 
                  const arma::mat& predictor,
                  double lambda, 
                  std::string method = "none") {
  
  arma::uword nrows = predictor.n_rows;
  arma::uword ncols = predictor.n_cols;
  arma::mat means_resp = arma::zeros<mat>(nrows, 1);
  arma::mat means_pred = arma::zeros<mat>(nrows, ncols);
  arma::mat vars = arma::square(predictor);
  arma::mat covars = arma::zeros<mat>(nrows, ncols);
  arma::mat betas = arma::zeros<mat>(nrows, ncols);
  arma::mat alphas = arma::zeros<mat>(nrows, 1);
  arma::mat resids = arma::zeros<mat>(nrows, 1);
  arma::mat varz = arma::ones<mat>(nrows, 1);
  arma::mat meanz = arma::zeros<mat>(nrows, 1);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  means_resp.row(0) = response.row(0);
  means_pred.row(0) = predictor.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    means_resp.row(it) = lambda1*response.row(it) + lambda*means_resp.row(it-1);
    means_pred.row(it) = lambda1*predictor.row(it) + lambda*means_pred.row(it-1);
    // cout << "Calculating vars: " << it << endl;
    vars.row(it) = lambda1*(vars.row(it) - arma::square(means_pred.row(it))) + lambda*vars.row(it-1);
    // cout << "Calculating covars: " << it << endl;
    covars.row(it) = lambda1*((response.row(it)-means_resp.row(it))*(predictor.row(it)-means_pred.row(it))) + lambda*covars.row(it-1);
    // cout << "Calculating betas: " << it << endl;
    // Calculate the alphas and betas.
    betas.row(it) = lambda1*covars.row(it)/vars.row(it) + lambda*betas.row(it-1);
    alphas.row(it) = lambda1*(means_resp.row(it) - arma::dot(betas.row(it), means_pred.row(it))) + lambda*alphas.row(it-1);
    // cout << "Calculating resids: " << it << endl;
    // Calculate the residuals.
    resids.row(it) = lambda1*(response.row(it) - arma::dot(betas.row(it), predictor.row(it))) + lambda*resids.row(it-1);
    // Calculate the mean and variance of the residuals.
    meanz.row(it) = lambda1*resids.row(it) + lambda*meanz.row(it-1);
    varz.row(it) = lambda1*arma::square(resids.row(it) - resids.row(it-1)) + lambda*varz.row(it-1);
  }  // end for
  
  if (method == "scale") {
    // Divide the residuals by their volatility
    resids = resids/sqrt(varz);
  } else if (method == "standardize") {
    // De-mean the residuals and divide them by their volatility
    resids = (resids - meanz)/sqrt(varz);
  }  // end if
  
  return join_rows(resids, alphas, betas);
  
}  // end run_rego

