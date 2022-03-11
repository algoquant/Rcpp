////////////////////////////////////////////////
// Rcpp functions for regression and rolling statistics
////////////////////////////////////////////////
// You can compile this file as follows:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/lm_arma.cpp")

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
Rcpp::List calc_lm(const arma::vec& response, const arma::mat& design) {
  // add column for intercept to explanatory matrix
  arma::mat designp = join_rows(ones(design.n_rows), design);
  int nrows = design.n_rows, ncols = designp.n_cols;
  int deg_free = (nrows - ncols);
  
  // fit the model response ~ design, and calculate alpha and beta coefficients
  arma::colvec coeff = arma::solve(designp, response);
  // calculate residuals
  arma::colvec residuals = response - designp*coeff;
  
  // calculate TSS, RSS, and ESS
  double tot_sumsq = (nrows-1)*arma::var(response);
  double res_sumsq = arma::dot(residuals, residuals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // calculate R-squared and F-statistic
  double rsquared = exp_sumsq/tot_sumsq;
  double fstat = (exp_sumsq*deg_free)/(res_sumsq*(ncols-1));
  // arma::rowvec stats=join_horiz(rsquared, fstat);
  Rcpp::NumericVector stats(2);
  stats(0) = rsquared;
  stats(1) = fstat;
  stats.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // calculate standard errors of beta coefficients
  arma::colvec stderr = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(designp)*designp)));
  // calculate t-values and p-values of beta coefficients
  arma::colvec beta_tvals = coeff/stderr;
  arma::colvec beta_pvals = 2*Rcpp::pt(-abs(wrap(beta_tvals)), deg_free);
  NumericMatrix coeffmat = wrap(join_rows(join_rows(join_rows(coeff, stderr), beta_tvals), beta_pvals));
  Rcpp::colnames(coeffmat) = Rcpp::CharacterVector::create("coeff", "stderr", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeffmat,
                            // Named("residuals") = residuals,
                            Named("z_score") = residuals(nrows-1)/arma::stddev(residuals),
                            Named("stats") = stats);
}  // end calc_lm



//' Perform a rolling regression over a time series of prices, 
//' \emph{RcppArmadillo}.
//' 
//' @param vectorv A numeric \emph{vector} of data.
//' @param wei_ghts A numeric \emph{vector} of weights.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vectorv}.
//'
//' @details The function \code{roll_scale()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vectorv, filter=wei_ghts, 
//'   method="convolution", sides=1)}, but it's about six times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vectorv <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' wei_ghts <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_zscores(vectorv=vectorv, wei_ghts=rev(wei_ghts))
//' # compare with original
//' all.equal(vectorv, as.numeric(weight_ed))
//' # Second example
//' # create exponentially decaying weights
//' wei_ghts <- exp(-0.2*1:11)
//' wei_ghts <- wei_ghts/sum(wei_ghts)
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_zscores(vectorv=vectorv, wei_ghts=rev(wei_ghts))
//' # calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=vectorv, filter=wei_ghts, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::vec& response, 
                       const arma::mat& design, 
                       const arma::uword& look_back) {
  arma::uword nrows = design.n_rows;
  arma::vec z_scores(nrows);
  arma::vec sub_response;
  arma::mat sub_design;
  Rcpp::List lm_list;
  
  // startup period
  for (uword it = 1; it < look_back; it++) {
    sub_response = response.subvec(0, it);
    sub_design = design.rows(0, it);
    lm_list = calc_lm(sub_response, sub_design);
    z_scores[it] = lm_list["z_score"];
  }  // end for
  
  // remaining periods
  for (uword it = look_back; it < nrows; it++) {
    sub_response = response.subvec(it-look_back+1, it);
    sub_design = design.rows(it-look_back+1, it);
    lm_list = calc_lm(sub_response, sub_design);
    z_scores[it] = lm_list["z_score"];
  }  // end for
  
  return z_scores;
}  // end roll_zscores




//' Calculate the rolling maximum and minimum over a \emph{vector} of data, 
//' using \emph{RcppArmadillo}.
//' 
//' @param vectorv A numeric \emph{vector} of data.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vectorv}.
//'
//' @details The function \code{roll_scale()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vectorv, filter=wei_ghts, 
//'   method="convolution", sides=1)}, but it's about six times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vectorv <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' wei_ghts <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_maxmin(vectorv=vectorv, wei_ghts=rev(wei_ghts))
//' # compare with original
//' all.equal(vectorv, as.numeric(weight_ed))
//' # Second example
//' # create exponentially decaying weights
//' wei_ghts <- exp(-0.2*1:11)
//' wei_ghts <- wei_ghts/sum(wei_ghts)
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_maxmin(vectorv=vectorv, wei_ghts=rev(wei_ghts))
//' # calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=vectorv, filter=wei_ghts, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_maxmin(const arma::vec& vectorv, 
                      const arma::uword& look_back) {
  arma::uword nrows = vectorv.size();
  arma::mat max_min(nrows, 2);
  arma::vec sub_vec;
  
  // startup period
  max_min(0, 0) = vectorv[0];
  max_min(0, 1) = vectorv[0];
  for (uword it = 1; it < look_back; it++) {
    sub_vec = vectorv.subvec(0, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  
  // remaining periods
  for (uword it = look_back; it < nrows; it++) {
    sub_vec = vectorv.subvec(it-look_back+1, it);
    max_min(it, 0) = sub_vec.max();
    max_min(it, 1) = sub_vec.min();
  }  // end for
  
  return max_min;
}  // end roll_maxmin



//' Cumulate the values of a numeric vector, and reset
//' the count to zero after every FALSE element.
//' 
//' @param vectorv A numeric \emph{vector} of data.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vectorv}.
//'
//' @details The function \code{roll_cum()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vectorv, filter=wei_ghts, 
//'   method="convolution", sides=1)}, but it's about \emph{6} times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vectorv <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' wei_ghts <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_cum(vectorv=vectorv, wei_ghts=rev(wei_ghts))
//' # compare with original
//' all.equal(vectorv, as.numeric(weight_ed))
//' # Second example
//' # create exponentially decaying weights
//' wei_ghts <- exp(-0.2*1:11)
//' wei_ghts <- wei_ghts/sum(wei_ghts)
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_cum(vectorv=vectorv, wei_ghts=rev(wei_ghts))
//' # calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=vectorv, filter=wei_ghts, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::ivec roll_cum(const arma::ivec& vectorv, const arma::sword& m_ax) {
  uword nrows = vectorv.n_elem;
  arma::ivec roll_sum(nrows);
  
  // startup period
  roll_sum[0] = vectorv[0];
  // remaining periods
  for (uword it = 1; it < nrows; it++) {
    roll_sum[it] = roll_sum[it-1] + vectorv[it];
    if (roll_sum[it] > m_ax)
      roll_sum[it] = m_ax;
    if (roll_sum[it] < (-m_ax))
      roll_sum[it] = -m_ax;
  }  // end for
  
  return roll_sum;
}  // end roll_cum


