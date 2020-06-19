////////////////////////////
// Functions to test C++ syntax with Armadillo
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_arma.cpp")

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


////////////////////////////////////////////////////////////
//' The function order_index() calculates the order index of a vector
//' using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::uvec order_index(const arma::vec& vec_tor) {
  return arma::sort_index(vec_tor);
}  // end order_index



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
Rcpp::List calc_lm(const arma::vec& res_ponse, const arma::mat& de_sign) {
  // Add column for intercept to explanatory matrix
  arma::mat design_p = join_rows(ones(de_sign.n_rows), de_sign);
  int num_rows = de_sign.n_rows, num_cols = design_p.n_cols;
  int deg_free = (num_rows - num_cols);
  
  // fit the model res_ponse ~ de_sign, and calculate alpha and beta coefficients
  arma::colvec co_eff = arma::solve(design_p, res_ponse);
  // Calculate residuals
  arma::colvec resid_uals = res_ponse - design_p*co_eff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (num_rows-1)*arma::var(res_ponse);
  double res_sumsq = arma::dot(resid_uals, resid_uals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate R-squared and F-statistic
  double r_squared = exp_sumsq/tot_sumsq;
  double f_stat = (exp_sumsq*deg_free)/(res_sumsq*(num_cols-1));
  // arma::rowvec stat_s=join_horiz(r_squared, f_stat);
  Rcpp::NumericVector stat_s(2);
  stat_s(0) = r_squared;
  stat_s(1) = f_stat;
  stat_s.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // Calculate standard errors of beta coefficients
  arma::colvec std_err = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(design_p)*design_p)));
  // Calculate t-values and p-values of beta coefficients
  arma::colvec beta_tvals = co_eff/std_err;
  arma::colvec beta_pvals = 2*Rcpp::pt(-abs(wrap(beta_tvals)), deg_free);
  NumericMatrix coeff_mat = wrap(join_rows(join_rows(join_rows(co_eff, std_err), beta_tvals), beta_pvals));
  Rcpp::colnames(coeff_mat) = Rcpp::CharacterVector::create("coeff", "std_err", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeff_mat,
                            // Named("residuals") = resid_uals,
                            Named("z_score") = resid_uals(num_rows-1)/arma::stddev(resid_uals),
                            Named("stats") = stat_s);
}  // end calc_lm



//' @export
// [[Rcpp::export]]
double calc_alpha(arma::mat& re_turns, 
                  const arma::mat& in_dex, 
                  const std::string& typ_e = "jensen",
                  const double& be_ta = 1.0) {
  // Initialize
  double al_pha = 0;

  // Calculate al_pha depending on typ_e
  if (typ_e == "jensen") {
    // Mean returns by columns
    Rcpp::NumericMatrix co_eff = calc_lm(re_turns, in_dex)["coefficients"];
    return co_eff(0, 0);
  } else if (typ_e == "wilcoxon") {
    re_turns = (re_turns - in_dex);
    // arma::uword num_rows = (re_turns.n_rows);
    arma::uvec rank_s = (arma::sort_index(arma::sort_index(re_turns)) + 1);
    return  dot(sign(re_turns), rank_s);
  } else if (typ_e == "kruskal_wallis") {
    arma::uword num_rows = (re_turns.n_rows);
    arma::mat combine_d = join_cols(re_turns, in_dex);
    arma::uvec rank_s = (arma::sort_index(arma::sort_index(combine_d)) + 1);
    // Apply regularized inverse to unit vector
    // weight_s = calc_inv(re_turns, max_eigen)*arma::ones(re_turns.n_cols);
    return  (0.0 + sum(rank_s.subvec(0, num_rows-1)) - sum(rank_s.subvec(num_rows, 2*num_rows-1)));
  } else if (typ_e == "rank") {
    // Mean returns by columns
    // arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Standard deviation by columns
    // arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    // sd_cols.replace(0, 1);
    // mean_cols = mean_cols/sd_cols;
    // al_pha equal to ranks of Sharpe
    // weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    // weight_s = (weight_s - arma::mean(weight_s));
    return al_pha;
  } else if (typ_e == "rankrob") {
    // Mean returns by columns
    // arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
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
    // al_pha equal to ranks of Sharpe
    // weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    // weight_s = conv_to< vec >::from(arma::sort_index(mean_cols));
    // pro_b;
    // weight_s = (weight_s - arma::mean(weight_s));
    return al_pha;
  } else {
    cout << "Warning: Incorrect typ_e argument: " << typ_e << endl;
    return al_pha;
  }  // end if
  
  return al_pha;
  
}  // end calc_alpha


// Calculate the top and bottom quantiles of the columns of a matrix
//' @export
// [[Rcpp::export]]
arma::vec calc_quant(const arma::mat& re_turns, 
                     const double& pro_b = 0.1) {
  
  arma::vec prob_s = {pro_b, 1-pro_b};
  // return arma::quantile(re_turns, prob_s, 0);
  return conv_to< colvec >::from(arma::sum(arma::quantile(re_turns, prob_s, 0), 0));
  
}  // end calc_quant



// Calculate the top and bottom quantiles of a vector,
// and return a vector of zeros, ones, and minus ones, 
// with the top quantile elements equal to 1, and the 
// bottom equal to -1.
//' @export
// [[Rcpp::export]]
arma::vec calc_top_bottom(const arma::vec& re_turns, 
                          const double& pro_b = 0.1) {
  
  // arma::uword num_elem = re_turns.n_elem;
  arma::vec weight_s = zeros(re_turns.n_elem);
  arma::vec prob_s = {pro_b, 1-pro_b};
  arma::vec quantile_s = arma::quantile(re_turns, prob_s);
  
  // Bottom quantile
  arma::uvec in_dex = find(re_turns <= quantile_s(0));
  weight_s(in_dex).fill(-1);
  
  // Top quantile
  // arma::rowvec quan_tile = arma::quantile(re_turns, pro_b);
  in_dex = find(re_turns >= quantile_s(1));
  // weight_s.elem(in_dex).ones();
  weight_s(in_dex).fill(1);
  
  return weight_s;
  
}  // end calc_top_bottom



// Calculate the top and bottom quantiles of the columns of a matrix.
// Add the quantiles for each column, and rank the vector of sums.
// Return a vector of zeros, ones, and minus ones, with the top 
// ranks equal to 1, and the bottom equal to -1.
//' @export
// [[Rcpp::export]]
arma::vec calc_top_bottom_columns(const arma::vec& re_turns, 
                      const double& pro_b = 0.1) {
  
  arma::uword num_cols = re_turns.n_cols;
  arma::uword thresh_old = round(pro_b*num_cols);
  arma::vec prob_s = {pro_b, 1-pro_b};
  arma::vec weight_s = zeros(num_cols);

  // return arma::quantile(re_turns, prob_s, 0);
  arma::rowvec sum_quant = arma::sum(arma::quantile(re_turns, prob_s, 0), 0);
  arma::uvec rank_s = (arma::sort_index(arma::sort_index(sum_quant)));
  // arma::uvec in_dex = find((rank_s <= thresh_old) || (rank_s >= (num_cols - thresh_old)));
  arma::uvec in_dex = find(rank_s >= (num_cols - thresh_old));
  // weight_s.elem(in_dex).ones();
  weight_s(in_dex).fill(1);
  in_dex = find(rank_s <= thresh_old);
  weight_s(in_dex).fill(-1);
  
  // in_dex = find((rank_s > thresh_old) || (rank_s < (num_cols - thresh_old)));
  // weight_s.elem(in_dex) = 0;
  // weight_s.head(thresh_old) = 1;
  // weight_s.tail(thresh_old) = 1;
  // return conv_to< colvec >::from(weight_s);
  return weight_s;
  
}  // end calc_top_bottom_columns



//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& re_turns, 
                       const std::string& typ_e = "quan_tile",
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
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
  } else if (typ_e == "quan_tile") {
    // Sum of quantiles by columns
    arma::vec prob_s = {pro_b, 1-pro_b};
    weight_s = conv_to< vec >::from(arma::sum(arma::quantile(re_turns, prob_s, 0), 0));
    // Weights equal to ranks
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(weight_s)));
    weight_s = (weight_s - arma::mean(weight_s));
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
    // Median returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
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
    // pro_b;
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


