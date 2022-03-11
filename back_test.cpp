////////////////////////////
// Copies of functions HighFreq::calc_weights() and 
// HighFreq::back_test() for testing.
// 
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/back_test.cpp")

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& returns,
                   double precision = 0.001, 
                   int max_eigen = 0) {
  
  arma::mat covmat = cov(returns);
  
  if (max_eigen == 0) {
    // Calculate the inverse using arma::pinv()
    return arma::pinv(covmat, precision);
  } else {
    // Calculate the inverse using eigen decomposition
    arma::mat eigen_vec;
    arma::vec eigen_val;
    arma::eig_sym(eigen_val, eigen_vec, covmat);
    eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
    eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
    return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  }  // end if
  
}  // end calc_inv




//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& returns, 
                       const std::string& model_type = "ranksharpe",
                       double precision = 0.001,
                       int max_eigen = 0,
                       const double& probv = 0.1,
                       const double& alpha = 0.0,
                       const bool scalit = true,
                       double vo_l = 0.01) {
  // Initialize
  arma::vec weights(returns.n_cols);
  if (max_eigen == 0)  max_eigen = returns.n_cols;
  
  // Calculate weights depending on model_type
  if (model_type == "ranksharpe") {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    // weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    // weights = (weights - arma::mean(weights));
    weights = meancols;
  } else if (model_type == "max_sharpe") {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, max_eigen);
    // weights = calc_inv(returns, max_eigen=max_eigen)*meancols;
    weights = calc_inv(returns, precision=precision, max_eigen=max_eigen)*meancols;
  } else if (model_type == "max_sharpe_median") {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::median(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, max_eigen);
    weights = calc_inv(returns, precision=precision, max_eigen=max_eigen)*meancols;
  } else if (model_type == "min_var") {
    // Apply regularized inverse to unit vector
    weights = calc_inv(returns, precision=precision, max_eigen=max_eigen)*arma::ones(returns.n_cols);
  } else if (model_type == "min_varpca") {
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, cov(returns));
    weights = eigen_vec.col(0);
  } else if (model_type == "rank") {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
  } else if (model_type == "rankrob") {
    // Median returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, max_eigen);
    // weights = calc_inv(returns, max_eigen)*meancols;
    // weights = calc_inv(returns, max_eigen)*meancols;
    // // Standard deviation by columns
    // arma::vec sd_cols = meancols;
    // for (arma::uword it=0; it < returns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((returns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    // probv;
    weights = (weights - arma::mean(weights));
  } else if (model_type == "quan_tile") {
    // Sum of quantiles for columns
    arma::vec probs = {probv, 1-probv};
    weights = conv_to< vec >::from(arma::sum(arma::quantile(returns, probs, 0), 0));
    // Weights equal to ranks
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(weights)));
    weights = (weights - arma::mean(weights));
  } else {
    cout << "Warning: Incorrect model_type argument: " << model_type << endl;
    return arma::ones(returns.n_cols);
  }  // end if
  
  if (scalit == TRUE) {
    // return weights/sqrt(sum(square(weights)));
    // return weights/sum(weights);
    // Returns of equally weighted portfolio
    // arma::vec meanrows = arma::mean(returns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = returns*weights;
    // Scale weights to equally weighted portfolio and return them
    // return weights*arma::stddev(arma::mean(returns, 1))/arma::stddev(returns*weights);
    // Scale weights so the resulting portfolio has a volatility equal to vo_l
    return weights*vo_l/arma::stddev(returns*weights);
  }  // end if
  
  return weights;
}  // end calc_weights




// First sort the datav and then unsort it back to original
//' @export
// [[Rcpp::export]]
arma::vec sort_back(const arma::vec& datav) {

  // Reverse sort index
  arma::uvec indeks = arma::sort_index(arma::sort_index(datav));
  // Sort the datav
  arma::vec sort_ed = arma::sort(datav);
  // Reverse sort the datav
  sort_ed = sort_ed.elem(indeks);
  
  return sort_ed;
}  // end sort_back


//' @export
// [[Rcpp::export]]
arma::uvec calc_ranks(const arma::vec& datav) {
  return (arma::sort_index(arma::sort_index(datav)) + 1);
}  // end calc_ranks



//' @export
// [[Rcpp::export]]
arma::vec calc_ranks_m(const arma::vec& datav) {
  // Ranks
  arma::vec ranks = conv_to< vec >::from(arma::sort_index(arma::sort_index(datav)));
  return (ranks - arma::mean(ranks));
}  // end calc_ranks_m



//' @export
// [[Rcpp::export]]
arma::vec weight_returns(const arma::mat& returns, const arma::vec& weights) {
  return returns*weights;
}  // end weight_returns




//' @export
// [[Rcpp::export]]
arma::mat back_test(const arma::mat& excess, // Portfolio excess returns
                    const arma::mat& returns, // Portfolio returns
                    const arma::uvec& startpoints, 
                    const arma::uvec& endpoints, 
                    const std::string& model_type = "ranksharpe",
                    double precision = 0.001,
                    int max_eigen = 0,
                    const double& probv = 0.1,
                    const double& alpha = 0,
                    const bool& scalit = true,
                    double vo_l = 0.01,
                    const double& coeff = 1.0,
                    const double& bid_offer = 0.0) {
  
  arma::vec weights(returns.n_cols);
  arma::vec weights_past = zeros(returns.n_cols);
  arma::mat pnls(returns.n_rows, returns.n_cols);
  
  // Perform loop over the endpoints
  for (arma::uword it = 1; it < endpoints.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate portfolio weights
    weights = coeff*calc_weights(excess.rows(startpoints(it-1), endpoints(it-1)), model_type, precision, max_eigen, probv, alpha, scalit, vo_l);
    // Calculate out-of-sample returns
    pnls.rows(endpoints(it-1)+1, endpoints(it)) = weights.t();
    // Add transaction costs
    // pnls.row(endpoints(it-1)+1) -= bid_offer*sum(abs(weights - weights_past))/2;
    // weights_past = weights;
  }  // end for
  
  // Return the strategy pnls
  return pnls;
  
}  // end back_test


