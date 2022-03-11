// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////////
// RcppArmadillo functions for Kalman filter
////////////////////////////


//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& matrixv, const arma::uword& max_eigen) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  
  arma::eig_sym(eigen_val, eigen_vec, cov(matrixv));
  eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
  eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
  // arma::mat eigen_valmat = diagmat(eigen_val);
  
  return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  
}  // end calc_inv


//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& returns, 
                       const arma::uword& max_eigen,
                       const double& alpha
) {
  arma::vec alpha_s(returns.n_cols);
  alpha_s.fill(alpha);
  arma::vec alphas_b(returns.n_cols);
  alphas_b.fill(1-alpha);
  arma::mat inverse = calc_inv(returns, max_eigen);
  arma::vec weights = arma::trans(arma::mean(returns, 0));
  arma::vec mean_s(weights.n_elem);
  mean_s.fill(arma::mean(weights));
  
  // shrink weights to the mean of weights
  weights = (alphas_b % weights + alpha_s % mean_s);
  // apply regularized inverse
  weights = inverse*weights;
  // scale weights and return them
  // return weights/sqrt(sum(square(weights)));
  return weights/sum(weights);
}  // end calc_weights



//' @export
// [[Rcpp::export]]
arma::mat roll_portf(const arma::mat& excess, // portfolio returns
                     const arma::mat& returns, // portfolio returns
                     const arma::uvec& startpoints, 
                     const arma::uvec& endpoints, 
                     const arma::uword& max_eigen,
                     const double& alpha
) {
  arma::vec sreturns = zeros(returns.n_rows);
  arma::vec weights(returns.n_cols);
  
  // perform a loop over the endpoints
  for (arma::uword it = 1; it < endpoints.size(); it++) {
    // cout << "it: " << it << endl;
    // subset the returns
    // arma::mat sub_returns = excess.rows(startpoints[it-1], endpoints[it-1]);
    // calculate portfolio weights
    weights = calc_weights(excess.rows(startpoints[it-1], endpoints[it-1]), max_eigen, alpha);
    // sub_returns = returns.rows(endpoints[it-1]+1, endpoints[it]);
    sreturns.subvec(endpoints[it-1]+1, endpoints[it]) = returns.rows(endpoints[it-1]+1, endpoints[it])*weights;
    // arma::mat foo = returns.rows(endpoints[it-1]+1, endpoints[it])*weights;
  }  // end for
  // return the strategy returns
  return sreturns;
}  // end roll_portf


