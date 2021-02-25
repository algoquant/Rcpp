////////////////////////////
// Function roll_reg() for performing rolling PCA and robust regressions
// Add option to perform rolling predictions using the regressions.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/roll_reg.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(const arma::vec& res_ponse, const arma::mat& de_sign) {
  // add column for intercept to explanatory matrix
  arma::mat design_p = join_rows(ones(de_sign.n_rows), de_sign);
  int num_rows = de_sign.n_rows, num_cols = design_p.n_cols;
  int deg_free = (num_rows - num_cols);
  
  // fit the model res_ponse ~ de_sign, and calculate alpha and beta coefficients
  arma::colvec co_eff = arma::solve(design_p, res_ponse);
  // calculate residuals
  arma::colvec resid_uals = res_ponse - design_p*co_eff;
  
  // calculate TSS, RSS, and ESS
  double tot_sumsq = (num_rows-1)*arma::var(res_ponse);
  double res_sumsq = arma::dot(resid_uals, resid_uals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // calculate R-squared and F-statistic
  double r_squared = exp_sumsq/tot_sumsq;
  double f_stat = (exp_sumsq*deg_free)/(res_sumsq*(num_cols-1));
  // arma::rowvec stat_s=join_horiz(r_squared, f_stat);
  Rcpp::NumericVector stat_s(2);
  stat_s(0) = r_squared;
  stat_s(1) = f_stat;
  stat_s.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // calculate standard errors of beta coefficients
  arma::colvec std_err = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(design_p)*design_p)));
  // calculate t-values and p-values of beta coefficients
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
arma::mat roll_reg(const arma::vec& res_ponse, 
                   const arma::mat& de_sign, 
                   const arma::uword& look_back) {
  arma::uword num_rows = de_sign.n_rows;
  arma::uword num_cols = de_sign.n_cols;
  arma::mat reg_stats(num_rows, (num_cols + 1));
  arma::vec sub_response;
  arma::mat sub_design;
  // Rcpp::List lm_list;
  
  // warmup period
  reg_stats.rows(0, num_cols+1) = zeros(num_cols+2, (num_cols + 1));
  for (arma::uword it = (num_cols+2); it < look_back; it++) {
    sub_response = res_ponse.subvec(0, it);
    sub_design = de_sign.rows(0, it);
    arma::mat co_eff = calc_lm(sub_response, sub_design)["coefficients"];
    reg_stats.row(it) = conv_to< rowvec >::from(co_eff.col(1));
  }  // end for
  
  // remaining periods
  for (arma::uword it = look_back; it < num_rows; it++) {
    sub_response = res_ponse.subvec(it-look_back+1, it);
    sub_design = de_sign.rows(it-look_back+1, it);
    arma::mat co_eff = calc_lm(sub_response, sub_design)["coefficients"];
    reg_stats.row(it) = conv_to< rowvec >::from(co_eff.col(1));
  }  // end for
  
  return reg_stats;
}  // end roll_reg



////////////////////////////////////////////////////////////
//' The function roll_wsum_mat() calculates the rolling weighted sum 
//' over a vector of data using \emph{RcppArmadillo}.
//' @export
// [[Rcpp::export]]
arma::mat roll_wsum_mat(arma::mat vec_tor, const arma::mat& weight_s) {
  uword n_rows = vec_tor.n_rows;
  uword look_back = weight_s.n_rows;
  // arma::mat vec_tor(n_rows);
  // arma::vec rev_weights = arma::reverse(weight_s);
  // arma::vec rev_weights = weight_s;
  
  // warmup period
  // vec_tor.n_rows(0, look_back-2) = vec_tor.subvec(0, look_back-2);
  
  // remaining periods
  for (uword it = look_back; it < n_rows; it++) {
    vec_tor.row(it-1) = trans(weight_s)*vec_tor.rows(it-look_back, it-1);
  }  // end for
  
  return vec_tor;
}  // end roll_wsum_mat


