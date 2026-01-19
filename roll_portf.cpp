// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
////////////////////////////
// RcppArmadillo functions for backtesting portfolio momentum strategies.
// 
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/roll_portf.cpp")
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


////////////////////////////////////////////////////////////
// Define C++ enum type for different methods of regularization,
// methodsfor calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum methodenum {moment, least_squares, quantile, nonparametric, regular, sharpem, 
                 maxsharpe, maxsharpemed, minvarlin, minvarquad, kellym, robustm, 
                 sumsq, sumone, voltarget, voleqw};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
methodenum calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return methodenum::moment;
  else if (method == "least_squares" || method == "l")
    return methodenum::least_squares;
  else if (method == "quantile" || method == "q")
    return methodenum::quantile;
  else if (method == "nonparametric" || method == "n")
    return methodenum::nonparametric;
  else if (method == "regular")
    return methodenum::regular;
  else if (method == "sharpem")
    return methodenum::sharpem;
  else if (method == "maxsharpe")
    return methodenum::maxsharpe;
  else if (method == "maxsharpemed")
    return methodenum::maxsharpemed;
  else if (method == "minvarlin")
    return methodenum::minvarlin;
  else if (method == "minvarquad")
    return methodenum::minvarquad;
  else if (method == "kellym")
    return methodenum::kellym;
  else if (method == "robustm")
    return methodenum::robustm;
  else if (method == "sumsq")
    return methodenum::sumsq;
  else if (method == "voltarget")
    return methodenum::voltarget;
  else 
    return methodenum::moment;
}  // end calc_method



// [[Rcpp::export]]
Rcpp::List param_portf(std::string method = "sharpem",  // Type of portfolio optimization model
                       double singmin = 1e-5, // Threshold level for discarding small singular values
                       arma::uword dimax = 0, // Number of eigen vectors for dimension reduction
                       arma::uword lagg = 1, // Lag of the calibration interval relative to the test interval
                       double confl = 0.1, // Confidence level for calculating the quantiles of returns
                       double alphac = 0.0, // Shrinkage intensity of returns
                       bool rankw = false, // Should the weights be ranked?
                       bool centerw = false, // Should the weights be centered?
                       std::string scalew = "voltarget", // Method for scaling the weights
                       double voltarget = 0.01) { // Volatility target for scaling the weights
  
  Rcpp::List controll = Rcpp::List::create(Rcpp::Named("method") = method,
                                           Rcpp::Named("singmin") = singmin,
                                           Rcpp::Named("dimax") = dimax,
                                           Rcpp::Named("lagg") = lagg,
                                           Rcpp::Named("confl") = confl,
                                           Rcpp::Named("alphac") = alphac,
                                           Rcpp::Named("rankw") = rankw,
                                           Rcpp::Named("centerw") = centerw,
                                           Rcpp::Named("scalew") = scalew,
                                           Rcpp::Named("voltarget") = voltarget);
  
  return controll;
  
}  // end param_portf


// [[Rcpp::export]]
arma::uvec calc_ranks_stl(arma::vec timeser) {
  
  std::size_t ndata = timeser.size();
  // size_t ndata = sizeof(timeser);
  
  // Define index of integers along timeser
  arma::uvec indeks(ndata);
  // Define the ranks of the vector elements
  arma::uvec ranks(ndata);
  // Fill the vectors with a sequence of consecutive integers
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(indeks.begin(), indeks.end(), 0);
  std::iota(ranks.begin(), ranks.end(), 0);
  
  // Calculate the permutation index by sorting the sequence of consecutive integers 
  // according to the order of timeser.
  std::sort(indeks.begin(), indeks.end(), 
            // Lambda comparison function defines sorting order.
            // The brackets [] are used to pass in variables from the outer scope of the lambdaf function.
            // The "&" passes the outer scope variables by reference.
            [&timeser](int i1, int i2) {return timeser[i1] < timeser[i2];});
  
  // Calculate the ranks (inverse permutation index) by sorting the sequence of consecutive integers 
  // according to the order of the permutation index indeks.
  std::sort(ranks.begin(), ranks.end(), 
            // Lambda comparison function defines sorting order.
            // The brackets [] are used to pass in variables from the outer scope of the lambda function.
            // The "&" passes the outer scope variables by reference.
            [&indeks](int i1, int i2) {return indeks[i1] < indeks[i2];});
  
  return ranks;
  
}  // end calc_ranks_stl




// [[Rcpp::export]]
arma::mat calc_covar(const arma::mat& timeser,
                     std::string method = "moment", 
                     double confl = 0.75) {
  
  arma::uword ncols = timeser.n_cols;
  // Return zeros if not enough data
  // if (timeser.n_rows < 3) {
  //   return arma::zeros<rowvec>(ncols);
  // }  // end if
  
  // Apply different calculation methods for covariance
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::cov(timeser);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, confl};
    arma::mat quantiles = arma::quantile(timeser, levels);
    return (quantiles.row(1) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // MAD
    arma::mat medians = arma::median(timeser);
    arma::mat mads(1, ncols);
    // Loop over columns of timeser
    for (arma::uword it = 0; it < ncols; it++) {
      mads.col(it) = arma::median(arma::abs(timeser.col(it) - arma::as_scalar(medians.col(it))));
    }  // end for
    // timeser.each_row() -= arma::median(timeser, 0);
    // return 1.4826*arma::median(arma::abs(timeser), 0);
    return 1.4826*mads;
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(ncols);
  }  // end default
  }  // end switch
  
}  // end calc_covar


// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& matrixv, 
                   arma::uword dimax = 0, // Number of eigen vectors for dimension reduction
                   double singmin = 0.0) { // Threshold for discarding small eigen values
  
  // Allocate Eigen variables
  arma::uword ncols = matrixv.n_cols;
  arma::vec eigenval(ncols, fill::zeros); // Eigen values
  arma::mat eigenvec(ncols, ncols, fill::zeros); // Eigen vectors
  // Calculate the eigen decomposition - the eigenvalues are in ascending order
  arma::eig_sym(eigenval, eigenvec, matrixv);
  // Calculate the number of non-small singular values
  arma::uword neigen = arma::sum(eigenval > singmin*arma::sum(eigenval));
  
  // If no regularization then set dimax to neigen
  if (dimax == 0) {
    // Set dimax
    dimax = neigen;
  } else {
    // Adjust dimax
    dimax = std::min(dimax, neigen);
  }  // end if
  
  // Remove all small singular values
  eigenval = eigenval.subvec(ncols-dimax, ncols-1);
  eigenvec = eigenvec.cols(ncols-dimax, ncols-1);
  
  // Calculate the reduced inverse from the eigen decomposition
  return eigenvec*arma::diagmat(1/eigenval)*eigenvec.t();
  
}  // end calc_inv



// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& returns, // Asset returns
                       Rcpp::List controll) { // List of portfolio optimization parameters
  
  // Unpack the control list of portfolio optimization parameters
  // Type of portfolio optimization model
  std::string method = Rcpp::as<std::string>(controll["method"]);
  // Threshold level for discarding small singular values
  double singmin = Rcpp::as<double>(controll["singmin"]);
  // Dimension reduction
  arma::uword dimax = Rcpp::as<int>(controll["dimax"]);
  // Confidence level for calculating the quantiles of returns
  double confl = Rcpp::as<double>(controll["confl"]);
  // Shrinkage intensity of returns
  double alphac = Rcpp::as<double>(controll["alphac"]);
  // Should the weights be ranked?
  bool rankw = Rcpp::as<int>(controll["rankw"]);
  // Should the weights be centered?
  bool centerw = Rcpp::as<int>(controll["centerw"]);
  // Method for scaling the weights
  std::string scalew = Rcpp::as<std::string>(controll["scalew"]);
  // Volatility target for scaling the weights
  double voltarget = Rcpp::as<double>(controll["voltarget"]);
  
  // Initialize the variables
  arma::uword ncols = returns.n_cols;
  arma::vec weightv(ncols, fill::zeros);
  // If no regularization then set dimax to ncols
  if (dimax == 0)  dimax = ncols;
  // Calculate the covariance matrix
  arma::mat covmat = calc_covar(returns);
  
  // Apply different calculation methods for the weights
  switch(calc_method(method)) {
  case methodenum::maxsharpe: {
    // Mean returns of columns
    arma::vec colmeans = arma::trans(arma::mean(returns, 0));
    // Shrink colmeans to the mean of returns
    colmeans = ((1-alphac)*colmeans + alphac*arma::mean(colmeans));
    // Calculate the weights using reduced inverse
    weightv = calc_inv(covmat, dimax, singmin)*colmeans;
    break;
  }  // end maxsharpe
  case methodenum::maxsharpemed: {
    // Median returns of columns
    arma::vec colmeans = arma::trans(arma::median(returns, 0));
    // Shrink colmeans to the median of returns
    colmeans = ((1-alphac)*colmeans + alphac*arma::median(colmeans));
    // Calculate the weights using reduced inverse
    weightv = calc_inv(covmat, dimax, singmin)*colmeans;
    break;
  }  // end maxsharpemed
  case methodenum::minvarlin: {
    // Minimum variance weights under linear constraint
    // Multiply reduced inverse times unit vector
    weightv = calc_inv(covmat, dimax, singmin)*arma::ones(ncols);
    break;
  }  // end minvarlin
  case methodenum::minvarquad: {
    // Minimum variance weights under quadratic constraint
    // Calculate highest order principal component
    arma::vec eigenval;
    arma::mat eigenvec;
    arma::eig_sym(eigenval, eigenvec, covmat);
    weightv = eigenvec.col(ncols-1);
    break;
  }  // end minvarquad
  case methodenum::sharpem: {
    // Momentum weights equal to Sharpe ratios
    // Mean returns of columns
    arma::vec colmeans = arma::trans(arma::mean(returns, 0));
    // Standard deviation of columns
    arma::vec colsd = arma::sqrt(covmat.diag());
    colsd.replace(0, 1);
    // Momentum weights equal to Sharpe ratios
    weightv = colmeans/colsd;
    break;
  }  // end sharpem
  case methodenum::kellym: {
    // Momentum weights equal to Kelly ratios
    // Mean returns of columns
    arma::vec colmeans = arma::trans(arma::mean(returns, 0));
    // Variance of columns
    arma::vec colvar = covmat.diag();
    colvar.replace(0, 1);
    // Momentum weights equal to Kelly ratios
    weightv = colmeans/colvar;
    break;
  }  // end kellym
  case methodenum::robustm: {
    // Momentum weights equal to robust Sharpe ratios
    // Median returns of columns
    arma::vec colmeans = arma::trans(arma::median(returns, 0));
    // Standard deviation of columns
    arma::vec colsd = arma::sqrt(covmat.diag());
    colsd.replace(0, 1);
    // Momentum weights equal to robust Sharpe ratios
    colmeans = colmeans/colsd;
    break;
  }  // end robustm
  case methodenum::quantile: {
    // Momentum weights equal to sum of quantiles for columns
    arma::vec levels = {confl, 1-confl};
    weightv = arma::conv_to<vec>::from(arma::sum(arma::quantile(returns, levels, 0), 0));
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(ncols);
  }  // end default
  }  // end switch
  
  if (rankw == TRUE) {
    // Convert the weights to their ranks
    weightv = arma::conv_to<vec>::from(calc_ranks_stl(weightv));
  }  // end if
  
  if (centerw == TRUE) {
    // Center the weights so their sum is equal to zero
    weightv = (weightv - arma::mean(weightv));
  }  // end if
  
  // Apply different scaling methods for the weights
  switch(calc_method(scalew)) {
  case methodenum::voltarget: {
    // Scale the weights so the portfolio has the volatility equal to voltarget
    weightv = weightv*voltarget/arma::stddev(returns*weightv);
    break;
  }  // end voltarget
  case methodenum::voleqw: {
    // Scale the weights to the volatility of the equal weight portfolio
    weightv = weightv*arma::stddev(arma::mean(returns, 1))/arma::stddev(returns*weightv);
    break;
  }  // end voleqw
  case methodenum::sumone: {
    // Scale the weights so their sum is equal to one
    weightv = weightv/arma::sum(weightv*arma::ones(ncols));
    break;
  }  // end sumone
  case methodenum::sumsq: {
    // Scale the weights so their sum of squares is equal to one
    weightv = weightv/std::sqrt(arma::sum(arma::square(weightv)));
    break;
  }  // end sumsq
  default : {
    // No scaling
    break;
  }  // end default
  }  // end switch
  
  return weightv;
  
}  // end calc_weights




// [[Rcpp::export]]
arma::mat roll_portf(const arma::mat& retx, // Asset excess returns
                     const arma::mat& retp, // Asset returns
                     arma::uvec startp, // Start points
                     arma::uvec endd, // End points
                     Rcpp::List controll, // List of portfolio optimization model parameters
                     double lambdaf = 0.0, // Decay factor for averaging the portfolio weights
                     double coeff = 1.0, // Multiplier of strategy returns
                     double bidask = 0.0) { // The bid-ask spread
  
  double lambda1 = 1-lambdaf;
  arma::uword nweights = retp.n_cols;
  arma::vec weightv(nweights, fill::zeros);
  arma::vec weightp = arma::ones(nweights)/std::sqrt(nweights); // Past weights
  arma::mat pnls = arma::zeros(retp.n_rows, 1);
  
  // Unpack the control list of portfolio optimization parameters
  // Lag of the calibration interval
  arma::uword lagg = Rcpp::as<int>(controll["lagg"]);
  
  // Perform loop over the end points
  for (arma::uword it = lagg; it < endd.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate the portfolio weights
    weightv = coeff*calc_weights(retx.rows(startp(it-lagg), endd(it-lagg)), controll);
    // cout << "calc_weights done" << endl;
    // Calculate the weights as the weighted sum with past weights
    weightv = lambda1*weightv + lambdaf*weightp;
    // Calculate out-of-sample returns
    pnls.rows(endd(it-1)+1, endd(it)) = retp.rows(endd(it-1)+1, endd(it))*weightv;
    // cout << "Out-of-sample returns done" << endl;
    // Add transaction costs
    pnls.row(endd(it-1)+1) -= bidask*sum(abs(weightv - weightp))/2;
    // Copy to the past weights
    weightp = weightv;
  }  // end for
  
  // Return the strategy pnls
  return pnls;
  
}  // end roll_portf


