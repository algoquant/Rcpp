////////////////////////////////////////////////
// Rcpp functions for regression and rolling statistics
////////////////////////////////////////////////
// You can compile this file as follows:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/lm_arma.cpp")

// Rcpp header with information for C++ compiler
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// The function calc_lm() performs multivariate linear regression, and 
// calculates the alpha and beta coefficients and their t-values and p-values, 
// and the R-squared and F-statistic.
// It uses RcppArmadillo.
// Adapted from:
// http://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
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



//' Perform a rolling regression over a time series of prices, 
//' \emph{RcppArmadillo}.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//' @param wei_ghts A numeric \emph{vector} of weights.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_scale()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vec_tor, filter=wei_ghts, 
//'   method="convolution", sides=1)}, but it's about six times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vec_tor <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' wei_ghts <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_zscores(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # compare with original
//' all.equal(vec_tor, as.numeric(weight_ed))
//' # Second example
//' # create exponentially decaying weights
//' wei_ghts <- exp(-0.2*1:11)
//' wei_ghts <- wei_ghts/sum(wei_ghts)
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_zscores(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=vec_tor, filter=wei_ghts, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::vec& res_ponse, 
                       const arma::mat& de_sign, 
                       const arma::uword& look_back) {
  arma::uword num_rows = de_sign.n_rows;
  arma::vec z_scores(num_rows);
  arma::vec sub_response;
  arma::mat sub_design;
  Rcpp::List lm_list;
  
  // startup period
  for (uword it = 1; it < look_back; it++) {
    sub_response = res_ponse.subvec(0, it);
    sub_design = de_sign.rows(0, it);
    lm_list = calc_lm(sub_response, sub_design);
    z_scores[it] = lm_list["z_score"];
  }  // end for
  
  // remaining periods
  for (uword it = look_back; it < num_rows; it++) {
    sub_response = res_ponse.subvec(it-look_back+1, it);
    sub_design = de_sign.rows(it-look_back+1, it);
    lm_list = calc_lm(sub_response, sub_design);
    z_scores[it] = lm_list["z_score"];
  }  // end for
  
  return z_scores;
}  // end roll_zscores




//' Calculate the rolling maximum and minimum over a \emph{vector} of data, 
//' using \emph{RcppArmadillo}.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_scale()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vec_tor, filter=wei_ghts, 
//'   method="convolution", sides=1)}, but it's about six times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vec_tor <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' wei_ghts <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_maxmin(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # compare with original
//' all.equal(vec_tor, as.numeric(weight_ed))
//' # Second example
//' # create exponentially decaying weights
//' wei_ghts <- exp(-0.2*1:11)
//' wei_ghts <- wei_ghts/sum(wei_ghts)
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_maxmin(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=vec_tor, filter=wei_ghts, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_maxmin(const arma::vec& vec_tor, 
                      const arma::uword& look_back) {
  arma::uword num_rows = vec_tor.size();
  arma::mat max_min(num_rows, 2);
  arma::vec sub_vec;
  
  // startup period
  max_min(0, 0) = vec_tor[0];
  max_min(0, 1) = vec_tor[0];
  for (uword it = 1; it < look_back; it++) {
    sub_vec = vec_tor.subvec(0, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  
  // remaining periods
  for (uword it = look_back; it < num_rows; it++) {
    sub_vec = vec_tor.subvec(it-look_back+1, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  
  return max_min;
}  // end roll_maxmin



//' Cumulate the values of a numeric vector, and reset
//' the count to zero after every FALSE element.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_cum()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vec_tor, filter=wei_ghts, 
//'   method="convolution", sides=1)}, but it's about \emph{6} times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vec_tor <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' wei_ghts <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_cum(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # compare with original
//' all.equal(vec_tor, as.numeric(weight_ed))
//' # Second example
//' # create exponentially decaying weights
//' wei_ghts <- exp(-0.2*1:11)
//' wei_ghts <- wei_ghts/sum(wei_ghts)
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_cum(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=vec_tor, filter=wei_ghts, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::ivec roll_cum(const arma::ivec& vec_tor, const arma::sword& m_ax) {
  uword len_gth = vec_tor.n_elem;
  arma::ivec roll_sum(len_gth);
  
  // startup period
  roll_sum[0] = vec_tor[0];
  // remaining periods
  for (uword it = 1; it < len_gth; it++) {
    roll_sum[it] = roll_sum[it-1] + vec_tor[it];
    if (roll_sum[it] > m_ax)
      roll_sum[it] = m_ax;
    if (roll_sum[it] < (-m_ax))
      roll_sum[it] = -m_ax;
  }  // end for
  
  return roll_sum;
}  // end roll_cum


