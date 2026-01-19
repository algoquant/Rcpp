// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Copyright: 2026 Jerzy Pawlowski

// This is another version of run_var_ohlc() created using Copilot.

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/run_var_ohlc.cpp")

#include "RcppArmadillo.h"
#include <vector>
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
// Use STL
using namespace std;
// For eigen solver SymEigsSolver
using namespace arma::newarp;

// Copilot prompt:
// Using the functions calc_var_ohlc() and run_var() as starting points,
// create the function run_var_ohlc(), for calculating the EMA variance from OHLC prices.
////////////////////////////////////////////////////////////
//' Calculate the exponential moving average variance of streaming \emph{OHLC}
//' prices using different price range estimators.
//'
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} of \emph{OHLC}
//'   prices.
//'   
//' @param \code{lambdaf} A decay factor which multiplies past variance estimates.
//'   
//' @param \code{method} A \emph{character string} representing the price range
//'   estimator for calculating the variance.  The estimators include:
//'   \itemize{
//'     \item "close" close-to-close estimator,
//'     \item "rogers_satchell" Rogers-Satchell estimator,
//'     \item "garman_klass" Garman-Klass estimator,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "yang_zhang" Yang-Zhang estimator,
//'    }
//'    (The default is the \code{method = "yang_zhang"}.)
//'    
//' @param \code{scalit} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scalit = TRUE}).
//'
//' @param \code{index} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{index = 0}).
//'   
//' @return A \emph{vector} with the same number of rows as the input argument
//'   \code{ohlc}, containing the running variance estimates.
//'
//' @details
//'   The function \code{run_var_ohlc()} calculates the exponential moving
//'   average variance from \emph{OHLC} prices using different price range
//'   estimators. It combines the methods for calculating the variance from OHLC
//'   prices from the function \code{calc_var_ohlc()} with the exponential
//'   weighting approach from the function \code{run_var()}.
//'
//' @examples
//' \dontrun{
//' # Extract the log OHLC prices of SPY
//' ohlc <- log(HighFreq::SPY)
//' # Calculate SPY variance without scaling to account for overnight price jumps
//' lambdaf <- 0.9 # Decay factor
//' vars <- HighFreq::run_var_ohlc(ohlc, lambdaf=lambdaf, method="yang_zhang", scalit=FALSE)
//' 
//' @export
// [[Rcpp::export]]
arma::colvec run_var_ohlc(const arma::mat& ohlc, 
                          double lambdaf,
                          std::string method = "yang_zhang", 
                          bool scalit = true,
                          arma::colvec index = 0) {
  
  arma::uword nrows = ohlc.n_rows;
  
  // Return zeros if not enough data
  if (nrows < 3) {
    return arma::zeros(nrows);
  }  // end if
  
  // Set time index to ones if not scaling
  if (!scalit || (index.n_rows == 1)) {
    index = arma::ones(nrows);
  }  // end if
  
  // Initialize variance estimates
  arma::colvec vars = arma::zeros(nrows);
  // double lambda2 = pow(lambdaf, 2);
  double lambda2 = lambdaf;
  double lambda21 = 1 - lambda2;
  double coeff = 0.34/(1.34 + (nrows+1)/(nrows-1));
  
  // Extract OHLC prices
  arma::colvec openp = ohlc.col(0);
  arma::colvec highp = ohlc.col(1);
  arma::colvec lowp = ohlc.col(2);
  arma::colvec closep = ohlc.col(3);
  
  // Initialize first variance estimate
  if (nrows >= 2) {
    double var0;
    if (method == "close") {
      var0 = pow(closep(1) - closep(0), 2);
    } else {
      // For other methods, use a simple estimate for the first value
      double opcl = (openp(1) - closep(0)) / index(1);
      double clop = (closep(1) - openp(1)) / index(1);
      double hiop = (highp(1) - openp(1)) / index(1);
      double lowop = (lowp(1) - openp(1)) / index(1);
      double hilow = (highp(1) - lowp(1)) / index(1);
      
      if (method == "rogers_satchell") {
        var0 = -(hiop - clop) * hiop - (lowop - clop) * lowop;
      } else if (method == "garman_klass") {
        var0 = 0.5 * pow(hilow, 2) - (2*log(2)-1) * pow(clop, 2);
      } else {
        var0 = pow(opcl, 2) + pow(clop, 2);
      }
    }
    vars(0) = var0;
  }
  
  // Perform loop over the rows to calculate running variance
  for (arma::uword it = 1; it < nrows; it++) {
    double current_var;
    
    if (method == "close") {
      current_var = pow((closep(it) - closep(it-1)), 2);
    } else {
      // Calculate current period returns
      double opcl = (openp(it) - closep(it-1)) / index(it);
      double clop = (closep(it) - openp(it)) / index(it);
      double hiop = (highp(it) - openp(it)) / index(it);
      double lowop = (lowp(it) - openp(it)) / index(it);
      double hilow = (highp(it) - lowp(it)) / index(it);
      double clhi = (closep(it) - highp(it)) / index(it);
      double clow = (closep(it) - lowp(it)) / index(it);
      
      if (method == "rogers_satchell") {
        current_var = -(clhi * hiop + clow * lowop);
      } else if (method == "garman_klass") {
        current_var = 0.5 * pow(hilow, 2) - (2*log(2)-1) * pow(clop, 2);
      } else if (method == "garman_klass_yz") {
        current_var = 0.5 * pow(hilow, 2) - (2*log(2)-1) * pow(clop, 2) + pow(opcl, 2);
      } else if (method == "yang_zhang") {
        current_var = pow(opcl, 2) + coeff * pow(clop, 2) + 
          (coeff-1) * (clhi * hiop + clow * lowop);
      } else {
        current_var = pow(opcl, 2) + pow(clop, 2);
      }
    }
    
    // Update running variance using exponential weighting
    vars(it) = lambda2 * vars(it-1) + lambda21 * current_var;
  }  // end for
  
  return vars;
  
}  // end run_var_ohlc 


