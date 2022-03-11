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


// [[Rcpp::export]]
arma::mat lagit(const arma::mat& tseries, 
                 arma::sword lagg = 1, 
                 bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_rows-1);
  arma::uword ncols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros<mat>(lagg, ncols), 
                             tseries.rows(0, nrows-lagg));
    } else {
      // Pad front with first element of tseries
      return arma::join_cols(arma::repmat(tseries.rows(0, 0), lagg, 1), 
                             tseries.rows(0, nrows-lagg));
    }  // end if
  } else {
    // Negative lag
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(tseries.rows(-lagg, nrows), 
                             arma::zeros<mat>(-lagg, ncols));
    } else {
      // Pad back with last element of tseries
      return arma::join_cols(tseries.rows(-lagg, nrows), 
                             arma::repmat(tseries.rows(nrows, nrows), -lagg, 1));
    }  // end if
  }  // end if
  
  // Old code below
  // if (lagg > 0)
  //   // Positive lag
  //   return arma::join_cols(arma::repelem(tseries.row(0), lagg, 1), 
  //                          tseries.rows(0, nrows-lagg));
  // else
  //   // Negative lag
  //   return arma::join_cols(tseries.rows(-lagg, nrows), 
  //                          arma::repelem(tseries.row(nrows), -lagg, 1));
  
}  // end lagit


// [[Rcpp::export]]
arma::mat run_mean(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat means = arma::zeros<mat>(nrows, 1);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  means(0) = tseries(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    means(it) = lambda1*tseries(it) + lambda*means(it-1);
  }  // end for
  
  return means;
  
}  // end run_mean



// [[Rcpp::export]]
arma::mat run_var(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat vars = arma::square(tseries);
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  means.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the variance as the weighted sum of squared returns minus the squared means
    vars.row(it) = lambda1*(vars.row(it) - arma::square(means.row(it))) + lambda*vars.row(it-1);
  }  // end for
  
  return vars;
  
}  // end run_var



// [[Rcpp::export]]
arma::mat run_covar(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat vars = arma::square(tseries);
  arma::mat covar = arma::zeros<mat>(nrows, 1);
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  means.row(0) = tseries.row(0);
  covar(0) = tseries(0, 0)*tseries(0, 1);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the covariance as the weighted sum of products of returns
    vars.row(it) = lambda1*(vars.row(it) - arma::square(means.row(it))) + lambda*vars.row(it-1);
    covar.row(it) = lambda1*((tseries(it, 0)-means(it, 0))*(tseries(it, 1)-means(it, 1))) + lambda*covar.row(it-1);
  }  // end for
  
  return arma::join_rows(covar, vars);
  
}  // end run_covar



////////////////////////////////////////////////////////////
//' Calculate the running variance of streaming \emph{OHLC} price data.
//' 
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} with \emph{OHLC}
//'   price data.
//'   
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
//'
//' @return A single-column \emph{matrix} of variance estimates, with the same
//'   number of rows as the input \code{ohlc} price data.
//'
//' @details
//'   The function \code{run_var_ohlc()} calculates a single-column
//'   \emph{matrix} of variance estimates of streaming \emph{OHLC} price data.
//'   
//'   The function \code{run_var_ohlc()} calculates the variance from the
//'   differences between the \emph{Open}, \emph{High}, \emph{Low}, and
//'   \emph{Close} prices, using the \emph{Yang-Zhang} range volatility
//'   estimator:
//'   \deqn{
//'     \sigma^2_t = (1-\lambda) ((O_t - C_{t-1})^2 + 0.134 (C_t - O_t)^2 + 
//'     0.866 ((H_i - O_i) (H_i - C_i) + (L_i - O_i) (L_i - C_i))) + 
//'     \lambda \sigma^2_{t-1}
//'   }
//'   It recursively weighs the current variance estimate with the past
//'   estimates \eqn{\sigma^2_{t-1}}, using the decay factor \eqn{\lambda}.
//'
//'   The function \code{run_var_ohlc()} does not calculate the logarithm of
//'   the prices.
//'   So if the argument \code{ohlc} contains dollar prices then
//'   \code{run_var_ohlc()} calculates the dollar variance.
//'   If the argument \code{ohlc} contains the log prices then
//'   \code{run_var_ohlc()} calculates the percentage variance.
//'   
//'   The function \code{run_var_ohlc()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Extract the log OHLC prices of VTI
//' ohlc <- log(rutils::etfenv$VTI)
//' # Calculate the running variance
//' var_running <- HighFreq::run_var_ohlc(ohlc, lambda=0.8)
//' # Calculate the rolling variance
//' var_rolling <- HighFreq::roll_var_ohlc(ohlc, look_back=5, method="yang_zhang", scale=FALSE)
//' datav <- cbind(var_running, var_rolling)
//' colnames(datav) <- c("running", "rolling")
//' colnamev <- colnames(datav)
//' datav <- xts::xts(datav, index(ohlc))
//' # dygraph plot of VTI running versus rolling volatility
//' dygraphs::dygraph(sqrt(datav[-(1:111), ]), main="Running and Rolling Volatility of VTI") %>%
//'   dyOptions(colors=c("red", "blue"), strokeWidth=1) %>%
//'   dyLegend(show="always", width=500)
//' # Compare the speed of running versus rolling volatility
//' library(microbenchmark)
//' summary(microbenchmark(
//'   running=HighFreq::run_var_ohlc(ohlc, lambda=0.8),
//'   rolling=HighFreq::roll_var_ohlc(ohlc, look_back=5, method="yang_zhang", scale=FALSE),
//'   times=10))[, c(1, 4, 5)]
//' }
//' @export
// [[Rcpp::export]]
arma::mat run_var_ohlc(const arma::mat& ohlc, 
                       double lambda) {
  
  // Allocate variance matrix
  arma::uword nrows = ohlc.n_rows;
  arma::mat vars = arma::zeros<mat>(nrows, 1);
  double lambda1 = 1-lambda;
  double coeff = 0.134;
  
  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::mat clo_se = ohlc.col(3);
  arma::mat open_close(clo_se.n_rows, 1);
  open_close = (ohlc.col(0) - lagit(clo_se, 1, false));
  arma::mat close_open = (clo_se - ohlc.col(0));
  arma::mat close_high = (clo_se - ohlc.col(1));
  arma::mat close_low = (clo_se - ohlc.col(2));
  arma::mat high_low = (ohlc.col(1) - ohlc.col(2));
  arma::mat high_open = (ohlc.col(1) - ohlc.col(0));
  arma::mat low_open = (ohlc.col(2) - ohlc.col(0));
  
  // Perform loop over the rows
  vars.row(0) = arma::square(open_close.row(0)) + coeff*arma::square(close_open.row(0)) +
    (coeff-1)*(close_high.row(0)*high_open.row(0) + close_low.row(0)*low_open.row(0));
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the variance as the weighted sum of squared returns minus the squared means
    vars.row(it) = lambda1*(arma::square(open_close.row(it)) + coeff*arma::square(close_open.row(it)) +
      (coeff-1)*(close_high.row(it)*high_open.row(it) + close_low.row(it)*low_open.row(it))) + lambda*vars.row(it-1);
  }  // end for
  
  return vars;
  
}  // end run_var_ohlc



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
  
  arma::uword nrows = tseries.n_rows;
  arma::mat maxmin = tseries;
  arma::mat means = tseries;
  double lambda1 = 1-lambda;
  
  if (calc_max) {
    // Perform loop over rows
    for (arma::uword it = 1; it < nrows; it++) {
      // Calculate the mean as a weighted sum
      means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
      // Calculate the max from a weighted sum
      maxmin.row(it) = arma::max(maxmin.row(it), means.row(it-1) + lambda*(maxmin.row(it-1) - means.row(it-1)));
    }  // end for
  } else {
    // Perform loop over rows
    for (arma::uword it = 1; it < nrows; it++) {
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
  
  arma::uword nrows = tseries.n_rows;
  arma::colvec maxs = tseries;
  
  // Perform loop over rows
  for (arma::uword it = 0; it < nrows; it++) {
    // Calculate the mean as a weighted sum
    // means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
    // Calculate the max from a weighted sum
    // maxs.row(it) = arma::max(maxs.row(it), means.row(it-1) + lambda*(maxs.row(it-1) - means.row(it-1)));
    maxs.row(it) = arma::max(tseries.row(it), tseries2.row(it));
  }  // end for
  
  return maxs;
  
}  // end armax

