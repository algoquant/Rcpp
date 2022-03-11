// Comment: switch command needs break; or return;

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_switch.cpp")


////////////////////////////////////////////////////////////
// Define C++ enum type for the different methods for regularization,
// calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum method {moment, least_squares, quantile, nonparametric, ranksharpe, 
              max_sharpe, max_sharpe_median, min_var, min_varpca, rank, rankrob};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
method calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return method::moment;
  else if (method == "least_squares" || method == "l")
    return method::least_squares;
  else if (method == "quantile" || method == "q")
    return method::quantile;
  else if (method == "nonparametric" || method == "n")
    return method::nonparametric;
  else if (method == "ranksharpe")
    return method::ranksharpe;
  else if (method == "max_sharpe")
    return method::max_sharpe;
  else if (method == "max_sharpe_median")
    return method::max_sharpe_median;
  else if (method == "min_var")
    return method::min_var;
  else if (method == "min_varpca")
    return method::min_varpca;
  else if (method == "rank")
    return method::rank;
  else if (method == "rankrob")
    return method::rankrob;
  else 
    return method::moment;
}  // end calc_method


//' @export
// [[Rcpp::export]]
arma::mat calc_inv(arma::mat returns,
                   double eigen_thresh = 0.001, 
                   int eigen_max = 0) {
  
  arma::mat covmat = arma::cov(returns);
  
  if (eigen_max == 0) {
    // Calculate the inverse using arma::pinv()
    return arma::pinv(covmat, eigen_thresh);
  } else {
    // Calculate the inverse using eigen decomposition
    arma::mat eigen_vec;
    arma::vec eigen_val;
    arma::eig_sym(eigen_val, eigen_vec, covmat);
    eigen_vec = eigen_vec.cols(eigen_vec.n_cols-eigen_max, eigen_vec.n_cols-1);
    eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-eigen_max, eigen_val.n_elem-1);
    return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  }  // end if
  
}  // end calc_inv


//' @export
// [[Rcpp::export]]
arma::vec calc_switch(arma::mat returns, // Portfolio returns
                       std::string method = "ranksharpe",
                       double eigen_thresh = 0.001,
                       int eigen_max = 0,
                       double confi_level = 0.1,
                       double alpha = 0.0,
                       bool scale = true,
                       double vol_target = 0.01) {
  // Initialize
  arma::vec weights(returns.n_cols);
  if (eigen_max == 0)  eigen_max = returns.n_cols;
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case method::ranksharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
  }  // end ranksharpe
  case method::max_sharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, eigen_max);
    // weights = calc_inv(returns, eigen_max=eigen_max)*meancols;
    weights = calc_inv(returns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*meancols;
  }  // end max_sharpe
  case method::max_sharpe_median: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::median(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, eigen_max);
    weights = calc_inv(returns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*meancols;
  }  // end max_sharpe_median
  case method::min_var: {
    cout << "Method parameter: " << method << endl;
    // Apply regularized inverse to unit vector
    weights = calc_inv(returns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*arma::ones(returns.n_cols);
  }  // end min_var
  case method::min_varpca: {
    cout << "Method parameter: " << method << endl;
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, arma::cov(returns));
    weights = eigen_vec.col(0);
  }  // end min_varpca
  case method::rank: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
  }  // end rank
  case method::rankrob: {
    cout << "Method parameter: " << method << endl;
    // Median returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, eigen_max);
    // weights = calc_inv(returns, eigen_max)*meancols;
    // weights = calc_inv(returns, eigen_max)*meancols;
    // // Standard deviation by columns
    // arma::vec sd_cols = meancols;
    // for (arma::uword it=0; it < returns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((returns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    // level;
    weights = (weights - arma::mean(weights));
  }  // end rankrob
  case method::quantile: {
    cout << "Method parameter: " << method << endl;
    // Sum of quantiles for columns
    arma::vec levels = {confi_level, 1-confi_level};
    weights = conv_to< vec >::from(arma::sum(arma::quantile(returns, levels, 0), 0));
    // Weights equal to ranks
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(weights)));
    weights = (weights - arma::mean(weights));
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(returns.n_cols);
  }  // end default
  }  // end switch
  
  return weights;
  
}  // end calc_switch



//' @export
// [[Rcpp::export]]
arma::vec calc_weights(arma::mat returns, // Portfolio returns
                       std::string method = "ranksharpe",
                       double eigen_thresh = 0.001,
                       int eigen_max = 0,
                       double confi_level = 0.1,
                       double alpha = 0.0,
                       bool scale = true,
                       double vol_target = 0.01) {
  // Initialize
  arma::vec weights(returns.n_cols, fill::zeros);
  if (eigen_max == 0)  eigen_max = returns.n_cols;
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case method::ranksharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end ranksharpe
  case method::max_sharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, eigen_max);
    // weights = calc_inv(returns, eigen_max=eigen_max)*meancols;
    weights = calc_inv(returns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*meancols;
    break;
  }  // end max_sharpe
  case method::max_sharpe_median: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::median(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, eigen_max);
    weights = calc_inv(returns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*meancols;
    break;
  }  // end max_sharpe_median
  case method::min_var: {
    cout << "Method parameter: " << method << endl;
    // Apply regularized inverse to unit vector
    weights = calc_inv(returns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*arma::ones(returns.n_cols);
    break;
  }  // end min_var
  case method::min_varpca: {
    cout << "Method parameter: " << method << endl;
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, arma::cov(returns));
    weights = eigen_vec.col(0);
    break;
  }  // end min_varpca
  case method::rank: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end rank
  case method::rankrob: {
    cout << "Method parameter: " << method << endl;
    // Median returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, eigen_max);
    // weights = calc_inv(returns, eigen_max)*meancols;
    // weights = calc_inv(returns, eigen_max)*meancols;
    // // Standard deviation by columns
    // arma::vec sd_cols = meancols;
    // for (arma::uword it=0; it < returns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((returns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    // level;
    weights = (weights - arma::mean(weights));
    break;
  }  // end rankrob
  case method::quantile: {
    cout << "Method parameter: " << method << endl;
    // Sum of quantiles for columns
    arma::vec levels = {confi_level, 1-confi_level};
    weights = conv_to< vec >::from(arma::sum(arma::quantile(returns, levels, 0), 0));
    // Weights equal to ranks
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(weights)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(returns.n_cols);
  }  // end default
  }  // end switch
  
  if (scale == TRUE) {
    // return weights/sqrt(sum(square(weights)));
    // return weights/sum(weights);
    // Returns of equally weighted portfolio
    // arma::vec meanrows = arma::mean(returns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = returns*weights;
    // Scale weights to equally weighted portfolio and return them
    // return weights*arma::stddev(arma::mean(returns, 1))/arma::stddev(returns*weights);
    // Scale weights so the resulting portfolio has a volatility equal to vol_target
    return weights*vol_target/arma::stddev(returns*weights);
  }  // end if
  
  return weights;
  
}  // end calc_weights


