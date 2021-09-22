// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// The function calc_eigen() calculates the eigen decomposition 
// of the matrix re_turns.
//' @export
// [[Rcpp::export]]
List calc_eigen(const arma::mat& mat_rix) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  arma::eig_sym(eigen_val, eigen_vec, cov(mat_rix));
  // reverse the order of elements from largest eigenvalue to smallest, similar to R
  return List::create(Named("values") = arma::flipud(eigen_val),
                      Named("vectors") = arma::fliplr(eigen_vec));
}  // end calc_eigen

