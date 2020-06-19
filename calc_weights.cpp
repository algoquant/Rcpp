#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


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
                       const std::string& typ_e = "max_sharpe",
                       int max_eigen = 1,
                       const double& quan_tile = 0.1,
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
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
  } else if (typ_e == "min_var") {
    // Apply regularized inverse to unit vector
    weight_s = calc_inv(re_turns, max_eigen)*arma::ones(re_turns.n_cols);
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
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    // weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
    // weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
    // // Standard deviation by columns
    // arma::vec sd_cols = mean_cols;
    // for (arma::uword it=0; it < re_turns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((re_turns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    // weight_s = conv_to< vec >::from(arma::sort_index(mean_cols));
    // quan_tile;
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
}  // end calc_weights



//' @export
// [[Rcpp::export]]
arma::vec sort_back(const arma::vec& da_ta) {

  // Reverse sort index
  arma::uvec in_dex = arma::sort_index(arma::sort_index(da_ta));
  // Sort the da_ta
  arma::vec sort_ed = arma::sort(da_ta);
  // Reverse sort the da_ta
  sort_ed = sort_ed.elem(in_dex);
  
  return sort_ed;
}  // end sort_back


//' @export
// [[Rcpp::export]]
arma::uvec calc_ranks(const arma::vec& da_ta) {
  return (arma::sort_index(arma::sort_index(da_ta)) + 1);
}  // end calc_ranks



//' @export
// [[Rcpp::export]]
arma::vec calc_ranks_m(const arma::vec& da_ta) {
  // Ranks
  arma::vec rank_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(da_ta)));
  return (rank_s - arma::mean(rank_s));
}  // end calc_ranks_m



