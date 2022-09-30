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


arma::uvec calc_ranks_stl(arma::vec tseries) {
  
  size_t ndata = tseries.size();
  // size_t ndata = sizeof(tseries);
  
  // Define index of integers along tseries
  arma::uvec indeks(ndata);
  // Define the ranks of the vector elements
  arma::uvec ranks(ndata);
  // Fill the vectors with a sequence of consecutive integers
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(indeks.begin(), indeks.end(), 0);
  std::iota(ranks.begin(), ranks.end(), 0);
  
  // Calculate the permutation index by sorting the sequence of consecutive integers 
  // according to the order of tseries.
  std::sort(indeks.begin(), indeks.end(), 
            // Lambda comparison function defines sorting order.
            // The brackets [] are used to pass in variables from the outer scope of the lambda function.
            // The "&" passes the outer scope variables by reference.
            [&tseries](int i1, int i2) {return tseries[i1] < tseries[i2];});
  
  // Calculate the ranks (inverse permutation index) by sorting the sequence of consecutive integers 
  // according to the order of the permutation index indeks.
  std::sort(ranks.begin(), ranks.end(), 
            // Lambda comparison function defines sorting order.
            // The brackets [] are used to pass in variables from the outer scope of the lambda function.
            // The "&" passes the outer scope variables by reference.
            [&indeks](int i1, int i2) {return indeks[i1] < indeks[i2];});
  
  return ranks;
  
}  // end calc_ranks_stl


arma::mat calc_covar(const arma::mat& tseries,
                     std::string method = "moment", 
                     double confl = 0.75) {
  
  double ncols = tseries.n_cols;
  // Return zeros if not enough data
  // if (tseries.n_rows < 3) {
  //   return arma::zeros<rowvec>(ncols);
  // }  // end if
  
  // Apply different calculation methods for covariance
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::cov(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(1) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // MAD
    arma::mat medians = arma::median(tseries);
    arma::mat mads(1, ncols);
    // Loop over columns of tseries
    for (arma::uword it = 0; it < ncols; it++) {
      mads.col(it) = arma::median(arma::abs(tseries.col(it) - arma::as_scalar(medians.col(it))));
    }  // end for
    // tseries.each_row() -= arma::median(tseries, 0);
    // return 1.4826*arma::median(arma::abs(tseries), 0);
    return 1.4826*mads;
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(ncols);
  }  // end default
  }  // end switch
  
}  // end calc_covar



arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword dimax = 0) {
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD
  arma::svd(svdu, svdval, svdv, tseries);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > eigen_thresh*arma::sum(svdval));
  
  if (dimax == 0) {
    // Set dimax
    dimax = svdnum - 1;
  } else {
    // Adjust dimax
    dimax = std::min(dimax - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, dimax);
  svdu = svdu.cols(0, dimax);
  svdv = svdv.cols(0, dimax);
  
  // Calculate the regularized inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_inv


arma::vec calc_weights(const arma::mat& returns, // Asset returns
                       std::string method = "maxsharpe", // Method for calculating the weights
                       double eigen_thresh = 1e-5, // Threshold level for discarding small singular values
                       arma::uword dimax = 0, // Dimension reduction
                       double confl = 0.1, // Confidence level for calculating the quantiles of returns
                       double alpha = 0.0, // Return shrinkage intensity
                       bool rankw = false, // Rank the weights
                       bool centerw = false, // Center the weights
                       std::string scalew = "voltarget", // Method for scaling the weights
                       double vol_target = 0.01) { // Target volatility for scaling the weights
  
  // Initialize
  arma::uword ncols = returns.n_cols;
  arma::vec weights(ncols, fill::zeros);
  // No regularization so set dimax to ncols
  if (dimax == 0)  dimax = ncols;
  // Calculate the covariance matrix
  arma::mat covmat = calc_covar(returns);
  
  // Apply different calculation methods for weights
  switch(calc_method(method)) {
  case methodenum::maxsharpe: {
    // Mean returns of columns
    arma::vec colmeans = arma::trans(arma::mean(returns, 0));
    // Shrink colmeans to the mean of returns
    colmeans = ((1-alpha)*colmeans + alpha*arma::mean(colmeans));
    // Calculate weights using regularized inverse
    weights = calc_inv(covmat, eigen_thresh, dimax)*colmeans;
    break;
  }  // end maxsharpe
  case methodenum::maxsharpemed: {
    // Median returns of columns
    arma::vec colmeans = arma::trans(arma::median(returns, 0));
    // Shrink colmeans to the median of returns
    colmeans = ((1-alpha)*colmeans + alpha*arma::median(colmeans));
    // Calculate weights using regularized inverse
    weights = calc_inv(covmat, eigen_thresh, dimax)*colmeans;
    break;
  }  // end maxsharpemed
  case methodenum::minvarlin: {
    // Minimum variance weights under linear constraint
    // Multiply regularized inverse times unit vector
    weights = calc_inv(covmat, eigen_thresh, dimax)*arma::ones(ncols);
    break;
  }  // end minvarlin
  case methodenum::minvarquad: {
    // Minimum variance weights under quadratic constraint
    // Calculate highest order principal component
    arma::vec eigenval;
    arma::mat eigenvec;
    arma::eig_sym(eigenval, eigenvec, covmat);
    weights = eigenvec.col(ncols-1);
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
    weights = colmeans/colsd;
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
    weights = colmeans/colvar;
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
    weights = arma::conv_to<vec>::from(arma::sum(arma::quantile(returns, levels, 0), 0));
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(ncols);
  }  // end default
  }  // end switch
  
  if (rankw == TRUE) {
    // Convert the weights to their ranks
    weights = arma::conv_to<vec>::from(calc_ranks_stl(weights));
  }  // end if
  
  if (centerw == TRUE) {
    // Center the weights so their sum is equal to zero
    weights = (weights - arma::mean(weights));
  }  // end if
  
  // Apply different scaling methods for weights
  switch(calc_method(scalew)) {
  case methodenum::voltarget: {
    // Scale the weights so the portfolio has the volatility equal to vol_target
    weights = weights*vol_target/arma::stddev(returns*weights);
    break;
  }  // end voltarget
  case methodenum::voleqw: {
    // Scale the weights to the volatility of the equal weight portfolio
    weights = weights*arma::stddev(arma::mean(returns, 1))/arma::stddev(returns*weights);
    break;
  }  // end voleqw
  case methodenum::sumone: {
    // Scale the weights so their sum of squares is equal to one
    weights = weights/arma::sum(weights*arma::ones(ncols));
    break;
  }  // end sumone
  case methodenum::sumsq: {
    // Scale the weights so their sum of squares is equal to one
    weights = weights/std::sqrt(arma::sum(square(weights)));
    break;
  }  // end sumsq
  default : {
    // No scaling
    break;
  }  // end default
  }  // end switch
  
  return weights;
  
}  // end calc_weights



// [[Rcpp::export]]
arma::mat back_test(const arma::mat& excess, // Asset excess returns
                    const arma::mat& returns, // Asset returns
                    arma::uvec startp, // Start points
                    arma::uvec endp, // End points
                    double lambda = 0.0, // Decay factor for averaging the portfolio weights
                    std::string method = "sharpem", // Method for calculating the weights
                    double eigen_thresh = 1e-5, // Threshold level for discarding small singular values
                    arma::uword dimax = 0, // Dimension reduction
                    double confl = 0.1, // Confidence level for calculating the quantiles of returns
                    double alpha = 0.0, // Return shrinkage intensity
                    bool rankw = false, // Rank the weights
                    bool centerw = false, // Center the weights
                    std::string scalew = "voltarget", // Method for scaling the weights
                    double vol_target = 0.01, // Target volatility for scaling the weights
                    double coeff = 1.0, // Multiplier strategy returns
                    double bid_offer = 0.0) {
  
  double lambda1 = 1-lambda;
  arma::uword nweights = returns.n_cols;
  arma::vec weights(nweights, fill::zeros);
  arma::vec weights_past = ones(nweights)/sqrt(nweights);
  arma::mat pnls = zeros(returns.n_rows, 1);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < endp.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate portfolio weights
    weights = coeff*calc_weights(excess.rows(startp(it-1), endp(it-1)), method, eigen_thresh, dimax, 
                                 confl, alpha, rankw, centerw, scalew, vol_target);
    // Calculate the weights as the weighted sum with past weights
    weights = lambda1*weights + lambda*weights_past;
    // Calculate out-of-sample returns
    pnls.rows(endp(it-1)+1, endp(it)) = returns.rows(endp(it-1)+1, endp(it))*weights;
    // Add transaction costs
    pnls.row(endp(it-1)+1) -= bid_offer*sum(abs(weights - weights_past))/2;
    // Copy the weights
    weights_past = weights;
  }  // end for
  
  // Return the strategy pnls
  return pnls;
  
}  // end back_test

