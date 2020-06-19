////////////////////////////
// Function back_test() for backtesting
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/back_test.cpp")

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

// Below are functions copied from package HighFreq

//' @export
// [[Rcpp::export]]
arma::mat calc_inv_t(const arma::mat& re_turns, 
                     const arma::uword& max_eigen = 0,
                     const double& to_l = 0.001) {
  
  arma::mat cov_mat = cov(re_turns);
  
  if (max_eigen == 0) {
    // Calculate the inverse using arma::pinv()
    return arma::pinv(cov_mat, to_l);
  } else {
    // Calculate the inverse using eigen decomposition
    arma::mat eigen_vec;
    arma::vec eigen_val;
    arma::eig_sym(eigen_val, eigen_vec, cov_mat);
    eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
    eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
    return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  }  // end if

}  // end calc_inv_t




//' @export
// [[Rcpp::export]]
arma::vec calc_weights_t(const arma::mat& re_turns, 
                       const std::string& typ_e = "max_sharpe",
                       int max_eigen = 1,
                       const double& pro_b = 0.1,
                       const double& al_pha = 0.0,
                       const bool scal_e = true) {
  // Initialize
  arma::vec weight_s(re_turns.n_cols);
  if (max_eigen == 1)  max_eigen = re_turns.n_cols;
  
  // Calculate weights depending on typ_e
  if (typ_e == "max_sharpe") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv_t(re_turns, max_eigen);
    weight_s = calc_inv_t(re_turns, max_eigen)*mean_cols;
  } else if (typ_e == "max_sharpe_median") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::median(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv_t(re_turns, max_eigen);
    weight_s = calc_inv_t(re_turns, max_eigen)*mean_cols;
  } else if (typ_e == "min_var") {
    // Apply regularized inverse to unit vector
    weight_s = calc_inv_t(re_turns, max_eigen)*arma::ones(re_turns.n_cols);
  } else if (typ_e == "min_varpca") {
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, cov(re_turns));
    weight_s = eigen_vec.col(0);
  } else if (typ_e == "rank") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    weight_s = (weight_s - arma::mean(weight_s));
  } else if (typ_e == "rankrob") {
    // Median returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv_t(re_turns, max_eigen);
    // weight_s = calc_inv_t(re_turns, max_eigen)*mean_cols;
    // weight_s = calc_inv_t(re_turns, max_eigen)*mean_cols;
    // // Standard deviation by columns
    // arma::vec sd_cols = mean_cols;
    // for (arma::uword it=0; it < re_turns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((re_turns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    // pro_b;
    weight_s = (weight_s - arma::mean(weight_s));
  } else if (typ_e == "quan_tile") {
    // Sum of quantiles for columns
    arma::vec prob_s = {pro_b, 1-pro_b};
    weight_s = conv_to< vec >::from(arma::sum(arma::quantile(re_turns, prob_s, 0), 0));
    // Weights equal to ranks
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(weight_s)));
    weight_s = (weight_s - arma::mean(weight_s));
  } else {
    cout << "Warning: Incorrect typ_e argument: " << typ_e << endl;
    return arma::ones(re_turns.n_cols);
  }  // end if
  
  if (scal_e == TRUE) {
    // Returns of equally weighted portfolio
    // arma::vec mean_rows = arma::mean(re_turns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = re_turns*weight_s;
    // Scale weight_s to equally weighted portfolio and return them
    // Return weight_s/sqrt(sum(square(weight_s)));
    // Return weight_s/sum(weight_s);
    return weight_s*arma::stddev(arma::mean(re_turns, 1))/arma::stddev(re_turns*weight_s);
  }  // end if
  
  return weight_s;
}  // end calc_weights_t



//' @export
// [[Rcpp::export]]
arma::mat back_test_t(const arma::mat& ex_cess, // Portfolio excess returns
                    const arma::mat& re_turns, // Portfolio returns
                    const arma::uvec& start_points, 
                    const arma::uvec& end_points, 
                    const std::string& typ_e = "max_sharpe",
                    const arma::uword& max_eigen = 1,
                    const double& pro_b = 0.1,
                    const double& al_pha = 0,
                    const bool& scal_e = true,
                    const double& co_eff = 1.0,
                    const double& bid_offer = 0.0) {
  arma::vec pnl_s = zeros(re_turns.n_rows);
  arma::vec weights_past = zeros(re_turns.n_cols);
  arma::vec weight_s(re_turns.n_cols);
  
  // Perform loop over the end_points
  for (arma::uword it=1; it < end_points.size(); it++) {
    cout << "it: " << it << endl;
    // Calculate portfolio weights
    weight_s = co_eff*calc_weights_t(ex_cess.rows(start_points(it-1), end_points(it-1)), typ_e, max_eigen, pro_b, al_pha, scal_e);
    cout << "weights: " << weight_s << endl;
    // Calculate out-of-sample returns
    pnl_s.subvec(end_points(it-1)+1, end_points(it)) = re_turns.rows(end_points(it-1)+1, end_points(it))*weight_s;
    // Add transaction costs
    pnl_s.row(end_points(it-1)+1) -= bid_offer*sum(abs(weight_s - weights_past))/2;
    weights_past = weight_s;
  }  // end for
  // Return the strategy returns
  return pnl_s;
}  // end back_test_t

