// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
////////////////////////////
// The function sim_momentum() simulates a momentum strategy.
// 
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/sim_momentum.cpp")
// 
// Copyright: 2026 Jerzy Pawlowski
////////////////////////////

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
// Use STL
using namespace std;


// [[Rcpp::export]]
arma::mat sim_momentum(const arma::mat& rets, // Asset returns
                         const double lambdaf, // Returns decay factor
                         const double volt = 0.01) { // Volatility target for scaling the weights
  
  arma::uword nrows = rets.n_rows;
  arma::uword ncols = rets.n_cols;
  arma::rowvec retsc; // Row of scaled asset returns
  arma::rowvec retm(ncols, fill::zeros); // EMA of asset returns
  arma::rowvec retv(ncols); // Row of asset returns
  arma::rowvec varv(ncols, fill::zeros); // Variance of asset returns
  arma::rowvec weightv(ncols, fill::zeros); // Row of portfolio weights
  double vol2 = pow(volt, 2);
  // arma::colvec varportf(nrows); // Portfolio variance
  arma::mat invmat(ncols, ncols, fill::ones); // Inverse covariance matrix
  arma::mat weightm(nrows, ncols, fill::zeros); // Portfolio weights
  arma::colvec pnls(nrows, fill::zeros); // Strategy PnLs
  double pnlm = 0; // EMA of strategy PnLs
  double pnlvar = 1e-8; // Variance of strategy PnLs
  double pnlvol = std::sqrt(pnlvar); // Volatility of strategy PnLs
  double weightd; // Sum of squared weights
  double lambda1 = 1-lambdaf;
  double lambda2 = pow(lambdaf, 2);
  double lambda21 = 1-lambda2;
  
  // Initialize EMA of returns with first row of returns
  retm = rets.row(0);
  // Calculate the centered returns
  arma::rowvec retc = retm;
  // Initialize the variance of the returns
  varv = arma::square(rets.row(1) - rets.row(0));
  
  if (!(rets.has_nan() || rets.has_inf())) {
    // No NA or Inf values
    for (arma::uword rn = 1; rn < nrows; rn++) {
      // Calculate the Kelly weights
      weightv = retm*vol2/(varv + 1e-8); // Add small number to prevent division by zero
      // Scale the weights so their sum of squares is equal to one
      // weightd = std::sqrt(arma::sum(arma::square(weightv)));
      // weightv = weightv/weightd;
      // Scale the weights by the volatility target divided by the PnLs volatility
      weightv = weightv*volt/std::sqrt(pnlvar);
      // Calculate the strategy PnLs - the products of the weights times the returns
      retv = rets.row(rn);
      pnls(rn) = arma::dot(retv, weightv);
      // Update the PnL variance and EMA
      pnlvar = lambda2*pnlvar + lambda21*pow(pnls(rn) - pnlm, 2);
      pnlm = lambdaf*pnlm + lambda1*pnls(rn);
      // Calculate the centered returns
      retc = (retv - retm);
      // Update the variance of the returns
      varv = lambda2*varv + lambda21*arma::trans(retc)*retc;
      // Update the EMA of the returns
      retm = lambdaf*retm + lambda1*retv;
    }  // end for
  } else {
    // cout << "NA or Inf values\n";
    // Loop over the columns because of the NA or Inf values.
    // Calculate covariance with NA or Inf values
    double retp = rets(0, 0);
    retm = rets.row(0);
    for (arma::uword rn = 1; rn < nrows; rn++) {
      // cout << "Row " << rn << "\n";
      retv = rets.row(rn);
      pnlvol = std::sqrt(pnlvar);
      pnls(rn) = 0; // Initialize the strategy PnL for this row
      for (arma::uword cn = 0; cn < ncols; cn++) {
        // cout << "Column " << cn << "\n";
        retp = retv(cn); // Current return
        if (std::isnan(retm(cn)) || std::isinf(retm(cn)) || std::isnan(retp) || std::isinf(retp)) {
          // cout << "NA or Inf values\n";
          // Set the mean and variance equal to the current return
          retm(cn) = retp;
          varv(cn) = pow(retp, 2);
          // Set the Kelly weight to zero
          weightv(cn) = 0;
        } else {
          // Calculate the Kelly weight
          // cout << "return: " << retm(cn) << "\n";
          weightv(cn) = retm(cn)*vol2/(varv(cn) + 1e-8); // Add small number to prevent division by zero
          // cout << "Weight: " << weightv(cn) << "\n";
          // Scale the weights by the volatility target divided by the PnLs volatility
          // cout << "Weights: " << sum(weightv) << "\n";
          weightv(cn) = weightv(cn)*volt/pnlvol;
          // cout << "Weights: " << sum(weightv) << "\n";
          // Calculate the strategy PnLs - the products of the weights times the returns
          pnls(rn) += weightv(cn)*retv(cn);
          // cout << "PnL: " << pnls(rn) << "\n";
          // Update the variance of returns
          varv(cn) = lambda2*varv(cn) + lambda21*pow(retp - retm(cn), 2);
          // cout << "Variance: " << varv(cn) << "\n";
          // Update the mean returns using the decay factor
          retm(cn) = lambdaf*retm(cn) + lambda1*retp;
          // cout << "Mean: " << retm(cn) << "\n";
        }  // end if
      }  // end column for
      // Scale the weights so their sum of squares is equal to one
      // weightd = std::sqrt(arma::sum(arma::square(weightv)));
      // weightv = weightv/weightd;
      // Update the PnL variance and EMA
      // cout << "C1\n";
      pnlvar = lambda2*pnlvar + lambda21*pow(pnls(rn) - pnlm, 2);
      // cout << "PnL variance: " << pnlvar << "\n";
      // cout << "C2\n";
      pnlm = lambdaf*pnlm + lambda1*pnls(rn);
      // cout << "PnL EMA: " << pnlm << "\n";

    }  // end row for
    // retm = rets;
  }  // end if
  
  // Return the strategy PnLs and the portfolio weights
  // return arma::join_rows(pnls, weightm, varportf);
  // return arma::join_rows(pnls, weightm);
  return pnls;
  
}  // end sim_momentum

