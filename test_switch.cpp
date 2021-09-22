// Comment: switch command needs break; or return;

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_switch.cpp")


////////////////////////////////////////////////////////////
// Define C++ enum type for the different methods for regularization,
// calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum meth_od {moment, least_squares, quantile, nonparametric, rank_sharpe, 
              max_sharpe, max_sharpe_median, min_var, min_varpca, rank, rankrob};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
meth_od calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return meth_od::moment;
  else if (method == "least_squares" || method == "l")
    return meth_od::least_squares;
  else if (method == "quantile" || method == "q")
    return meth_od::quantile;
  else if (method == "nonparametric" || method == "n")
    return meth_od::nonparametric;
  else if (method == "rank_sharpe")
    return meth_od::rank_sharpe;
  else if (method == "max_sharpe")
    return meth_od::max_sharpe;
  else if (method == "max_sharpe_median")
    return meth_od::max_sharpe_median;
  else if (method == "min_var")
    return meth_od::min_var;
  else if (method == "min_varpca")
    return meth_od::min_varpca;
  else if (method == "rank")
    return meth_od::rank;
  else if (method == "rankrob")
    return meth_od::rankrob;
  else 
    return meth_od::moment;
}  // end calc_method


//' @export
// [[Rcpp::export]]
arma::mat calc_inv(arma::mat re_turns,
                   double eigen_thresh = 0.001, 
                   int eigen_max = 0) {
  
  arma::mat cov_mat = arma::cov(re_turns);
  
  if (eigen_max == 0) {
    // Calculate the inverse using arma::pinv()
    return arma::pinv(cov_mat, eigen_thresh);
  } else {
    // Calculate the inverse using eigen decomposition
    arma::mat eigen_vec;
    arma::vec eigen_val;
    arma::eig_sym(eigen_val, eigen_vec, cov_mat);
    eigen_vec = eigen_vec.cols(eigen_vec.n_cols-eigen_max, eigen_vec.n_cols-1);
    eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-eigen_max, eigen_val.n_elem-1);
    return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  }  // end if
  
}  // end calc_inv


//' @export
// [[Rcpp::export]]
arma::vec calc_switch(arma::mat re_turns, // Portfolio returns
                       std::string method = "rank_sharpe",
                       double eigen_thresh = 0.001,
                       int eigen_max = 0,
                       double confi_level = 0.1,
                       double alpha = 0.0,
                       bool scale = true,
                       double vol_target = 0.01) {
  // Initialize
  arma::vec weight_s(re_turns.n_cols);
  if (eigen_max == 0)  eigen_max = re_turns.n_cols;
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case meth_od::rank_sharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    weight_s = (weight_s - arma::mean(weight_s));
  }  // end rank_sharpe
  case meth_od::max_sharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-alpha)*mean_cols + alpha*arma::mean(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, eigen_max);
    // weight_s = calc_inv(re_turns, eigen_max=eigen_max)*mean_cols;
    weight_s = calc_inv(re_turns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*mean_cols;
  }  // end max_sharpe
  case meth_od::max_sharpe_median: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-alpha)*mean_cols + alpha*arma::median(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, eigen_max);
    weight_s = calc_inv(re_turns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*mean_cols;
  }  // end max_sharpe_median
  case meth_od::min_var: {
    cout << "Method parameter: " << method << endl;
    // Apply regularized inverse to unit vector
    weight_s = calc_inv(re_turns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*arma::ones(re_turns.n_cols);
  }  // end min_var
  case meth_od::min_varpca: {
    cout << "Method parameter: " << method << endl;
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, arma::cov(re_turns));
    weight_s = eigen_vec.col(0);
  }  // end min_varpca
  case meth_od::rank: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    weight_s = (weight_s - arma::mean(weight_s));
  }  // end rank
  case meth_od::rankrob: {
    cout << "Method parameter: " << method << endl;
    // Median returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // mean_cols = ((1-alpha)*mean_cols + alpha*arma::mean(mean_cols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, eigen_max);
    // weight_s = calc_inv(re_turns, eigen_max)*mean_cols;
    // weight_s = calc_inv(re_turns, eigen_max)*mean_cols;
    // // Standard deviation by columns
    // arma::vec sd_cols = mean_cols;
    // for (arma::uword it=0; it < re_turns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((re_turns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    // level;
    weight_s = (weight_s - arma::mean(weight_s));
  }  // end rankrob
  case meth_od::quantile: {
    cout << "Method parameter: " << method << endl;
    // Sum of quantiles for columns
    arma::vec level_s = {confi_level, 1-confi_level};
    weight_s = conv_to< vec >::from(arma::sum(arma::quantile(re_turns, level_s, 0), 0));
    // Weights equal to ranks
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(weight_s)));
    weight_s = (weight_s - arma::mean(weight_s));
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(re_turns.n_cols);
  }  // end default
  }  // end switch
  
  return weight_s;
  
}  // end calc_switch



//' @export
// [[Rcpp::export]]
arma::vec calc_weights(arma::mat re_turns, // Portfolio returns
                       std::string method = "rank_sharpe",
                       double eigen_thresh = 0.001,
                       int eigen_max = 0,
                       double confi_level = 0.1,
                       double alpha = 0.0,
                       bool scale = true,
                       double vol_target = 0.01) {
  // Initialize
  arma::vec weight_s(re_turns.n_cols, fill::zeros);
  if (eigen_max == 0)  eigen_max = re_turns.n_cols;
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case meth_od::rank_sharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    weight_s = (weight_s - arma::mean(weight_s));
    break;
  }  // end rank_sharpe
  case meth_od::max_sharpe: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-alpha)*mean_cols + alpha*arma::mean(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, eigen_max);
    // weight_s = calc_inv(re_turns, eigen_max=eigen_max)*mean_cols;
    weight_s = calc_inv(re_turns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*mean_cols;
    break;
  }  // end max_sharpe
  case meth_od::max_sharpe_median: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-alpha)*mean_cols + alpha*arma::median(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, eigen_max);
    weight_s = calc_inv(re_turns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*mean_cols;
    break;
  }  // end max_sharpe_median
  case meth_od::min_var: {
    cout << "Method parameter: " << method << endl;
    // Apply regularized inverse to unit vector
    weight_s = calc_inv(re_turns, eigen_thresh=eigen_thresh, eigen_max=eigen_max)*arma::ones(re_turns.n_cols);
    break;
  }  // end min_var
  case meth_od::min_varpca: {
    cout << "Method parameter: " << method << endl;
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, arma::cov(re_turns));
    weight_s = eigen_vec.col(0);
    break;
  }  // end min_varpca
  case meth_od::rank: {
    cout << "Method parameter: " << method << endl;
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    weight_s = (weight_s - arma::mean(weight_s));
    break;
  }  // end rank
  case meth_od::rankrob: {
    cout << "Method parameter: " << method << endl;
    // Median returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // mean_cols = ((1-alpha)*mean_cols + alpha*arma::mean(mean_cols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, eigen_max);
    // weight_s = calc_inv(re_turns, eigen_max)*mean_cols;
    // weight_s = calc_inv(re_turns, eigen_max)*mean_cols;
    // // Standard deviation by columns
    // arma::vec sd_cols = mean_cols;
    // for (arma::uword it=0; it < re_turns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((re_turns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    // level;
    weight_s = (weight_s - arma::mean(weight_s));
    break;
  }  // end rankrob
  case meth_od::quantile: {
    cout << "Method parameter: " << method << endl;
    // Sum of quantiles for columns
    arma::vec level_s = {confi_level, 1-confi_level};
    weight_s = conv_to< vec >::from(arma::sum(arma::quantile(re_turns, level_s, 0), 0));
    // Weights equal to ranks
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(weight_s)));
    weight_s = (weight_s - arma::mean(weight_s));
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(re_turns.n_cols);
  }  // end default
  }  // end switch
  
  if (scale == TRUE) {
    // return weight_s/sqrt(sum(square(weight_s)));
    // return weight_s/sum(weight_s);
    // Returns of equally weighted portfolio
    // arma::vec mean_rows = arma::mean(re_turns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = re_turns*weight_s;
    // Scale weight_s to equally weighted portfolio and return them
    // return weight_s*arma::stddev(arma::mean(re_turns, 1))/arma::stddev(re_turns*weight_s);
    // Scale weight_s so the resulting portfolio has a volatility equal to vol_target
    return weight_s*vol_target/arma::stddev(re_turns*weight_s);
  }  // end if
  
  return weight_s;
  
}  // end calc_weights


