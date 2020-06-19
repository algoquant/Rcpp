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
arma::mat calc_inv(const arma::mat& mat_rix, const arma::uword& max_eigen) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  
  arma::eig_sym(eigen_val, eigen_vec, cov(mat_rix));
  eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
  eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
  // arma::mat eigen_valmat = diagmat(eigen_val);
  
  return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  
}  // end calc_inv


//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& re_turns, 
                       const arma::uword& max_eigen,
                       const double& al_pha
) {
  arma::vec alpha_s(re_turns.n_cols);
  alpha_s.fill(al_pha);
  arma::vec alphas_b(re_turns.n_cols);
  alphas_b.fill(1-al_pha);
  arma::mat in_verse = calc_inv(re_turns, max_eigen);
  arma::vec weight_s = arma::trans(arma::mean(re_turns, 0));
  arma::vec mean_s(weight_s.n_elem);
  mean_s.fill(arma::mean(weight_s));
  
  // shrink weight_s to the mean of weight_s
  weight_s = (alphas_b % weight_s + alpha_s % mean_s);
  // apply regularized inverse
  weight_s = in_verse*weight_s;
  // scale weight_s and return them
  // return weight_s/sqrt(sum(square(weight_s)));
  return weight_s/sum(weight_s);
}  // end calc_weights



//' @export
// [[Rcpp::export]]
arma::mat roll_portf(const arma::mat& ex_cess, // portfolio returns
                     const arma::mat& re_turns, // portfolio returns
                     const arma::uvec& start_points, 
                     const arma::uvec& end_points, 
                     const arma::uword& max_eigen,
                     const double& al_pha
) {
  arma::vec sre_turns = zeros(re_turns.n_rows);
  arma::vec weight_s(re_turns.n_cols);
  
  // perform a loop over the end_points
  for (arma::uword it = 1; it < end_points.size(); it++) {
    // cout << "it: " << it << endl;
    // subset the returns
    // arma::mat sub_returns = ex_cess.rows(start_points[it-1], end_points[it-1]);
    // calculate portfolio weights
    weight_s = calc_weights(ex_cess.rows(start_points[it-1], end_points[it-1]), max_eigen, al_pha);
    // sub_returns = re_turns.rows(end_points[it-1]+1, end_points[it]);
    sre_turns.subvec(end_points[it-1]+1, end_points[it]) = re_turns.rows(end_points[it-1]+1, end_points[it])*weight_s;
    // arma::mat foo = re_turns.rows(end_points[it-1]+1, end_points[it])*weight_s;
  }  // end for
  // return the strategy returns
  return sre_turns;
}  // end roll_portf


