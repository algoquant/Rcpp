////////////////////////////////////////////////
// Rcpp functions for regression and rolling statistics
////////////////////////////////////////////////
// You can compile this file as follows:
//  Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/lm_arma.cpp")

// Rcpp header with information for C++ compiler
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]



////////////////////////////////////////////////////////////
// Define C++ enum type for the different methods for regularization,
// calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum methodenum {moment, least_squares, quantile, nonparametric, regular, ranksharpe, 
                 max_sharpe, max_sharpe_median, min_var, min_varpca, rank, rankrob};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
methodenum calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return methodenum::moment;
  else if (method == "least_squares" || method == "l")
    return methodenum::least_squares;
  else if (method == "quantile" || method == "q")
    return methodenum::quantile;
  else if (method == "nonparametric" || method == "n")
    return methodenum::nonparametric;
  else if (method == "regular")
    return methodenum::regular;
  else if (method == "ranksharpe")
    return methodenum::ranksharpe;
  else if (method == "max_sharpe")
    return methodenum::max_sharpe;
  else if (method == "max_sharpe_median")
    return methodenum::max_sharpe_median;
  else if (method == "min_var")
    return methodenum::min_var;
  else if (method == "min_varpca")
    return methodenum::min_varpca;
  else if (method == "rank")
    return methodenum::rank;
  else if (method == "rankrob")
    return methodenum::rankrob;
  else 
    return methodenum::moment;
}  // end calc_method


//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword eigen_max = 0) {
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD
  arma::svd(svdu, svdval, svdv, tseries);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > eigen_thresh*arma::sum(svdval));
  
  if (eigen_max == 0) {
    // Set eigen_max
    eigen_max = svdnum - 1;
  } else {
    // Adjust eigen_max
    eigen_max = min(eigen_max - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, eigen_max);
  svdu = svdu.cols(0, eigen_max);
  svdv = svdv.cols(0, eigen_max);
  
  // Calculate the shrinkage inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_inv




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


//' @export
// [[Rcpp::export]]
arma::mat calc_reg(const arma::mat& response, 
                   const arma::mat& predictor,
                   bool intercept = true,
                   std::string method = "least_squares",
                   double eigen_thresh = 1e-5,
                   arma::uword eigen_max = 0,
                   double conf_lev = 0.1,
                   double alpha = 0.0) {
  
  // Add column for intercept to predictor matrix
  arma::uword nrows = predictor.n_rows;
  arma::mat predictori = predictor;
  if (intercept)
    predictori = join_rows(ones(nrows), predictor);
  
  arma::uword ncols = predictori.n_cols;
  arma::uword deg_free = (nrows - ncols);
  arma::vec coeff;
  arma::vec tvals;
  arma::mat reg_data = arma::zeros<mat>(2*ncols+1, 1);
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case methodenum::least_squares: {
    // Calculate regression coefficients for the model response ~ predictor
    coeff = arma::solve(predictori, response);
    break;
  }  // end least_squares
  case methodenum::regular: {
    // Calculate shrinkage regression coefficients
    coeff = calc_inv(predictori, eigen_thresh, eigen_max)*response;
    break;
  }  // end regular
  case methodenum::quantile: {
    // Not implemented yet
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return reg_data;
  }  // end default
  }  // end switch
  
  // Calculate residuals
  arma::mat residuals = response - predictori*coeff;
  
  // Calculate TSS, RSS, and ESS
  // double tot_sumsq = (nrows-1)*arma::var(response);
  double res_sumsq = arma::dot(residuals, residuals);
  // double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate standard errors of the beta coefficients
  arma::mat stderr = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(predictori)*predictori)));
  // Calculate t-values of the beta coefficients
  tvals = coeff/stderr;
  
  // Calculate z-score
  mat zscore = residuals(nrows-1, 0)/arma::stddev(residuals);
  
  // Combine regression data
  reg_data.rows(0, ncols-1) = coeff;
  reg_data.rows(ncols, 2*ncols-1) = tvals;
  reg_data.row(2*ncols) = zscore;
  
  return reg_data.t();
  
}  // end calc_reg


//' Perform a rolling regression over a time series of prices, 
//' \emph{RcppArmadillo}.
//' 
//' @param vectorv A numeric \emph{vector} of data.
//' @param weightv A numeric \emph{vector} of weights.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vectorv}.
//'
//' @details The function \code{roll_scale()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vectorv, filter=weightv, 
//'   method="convolution", sides=1)}, but it's about six times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vectorv <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' weightv <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weighted <- HighFreq::roll_zscores(vectorv=vectorv, weights=rev(weightv))
//' # compare with original
//' all.equal(vectorv, as.numeric(weighted))
//' # Second example
//' # create exponentially decaying weights
//' weightv <- exp(-0.2*1:11)
//' weightv <- weightv/sum(weightv)
//' # calculate rolling weighted sum
//' weighted <- HighFreq::roll_zscores(vectorv=vectorv, weights=rev(weightv))
//' # calculate rolling weighted sum using filter()
//' filtered <- filter(x=vectorv, filter=weightv, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filtered[-(1:11)]), as.numeric(weighted[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::vec& response, 
                       const arma::mat& design, 
                       const arma::uword& look_back) {
  arma::uword nrows = design.n_rows;
  arma::vec zscores(nrows);
  arma::vec sub_response;
  arma::mat sub_design;
  Rcpp::List lm_list;
  
  // startup period
  for (uword it = 1; it < look_back; it++) {
    sub_response = response.subvec(0, it);
    sub_design = design.rows(0, it);
    lm_list = calc_lm(sub_response, sub_design);
    zscores[it] = lm_list["z_score"];
  }  // end for
  
  // remaining periods
  for (uword it = look_back; it < nrows; it++) {
    sub_response = response.subvec(it-look_back+1, it);
    sub_design = design.rows(it-look_back+1, it);
    lm_list = calc_lm(sub_response, sub_design);
    zscores[it] = lm_list["z_score"];
  }  // end for
  
  return zscores;
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
//'   the standard \emph{R} function \code{filter(x=vectorv, filter=weightv, 
//'   method="convolution", sides=1)}, but it's about six times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vectorv <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' weightv <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weighted <- HighFreq::roll_maxmin(vectorv=vectorv, weights=rev(weightv))
//' # compare with original
//' all.equal(vectorv, as.numeric(weighted))
//' # Second example
//' # create exponentially decaying weights
//' weightv <- exp(-0.2*1:11)
//' weightv <- weightv/sum(weightv)
//' # calculate rolling weighted sum
//' weighted <- HighFreq::roll_maxmin(vectorv=vectorv, weights=rev(weightv))
//' # calculate rolling weighted sum using filter()
//' filtered <- filter(x=vectorv, filter=weightv, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filtered[-(1:11)]), as.numeric(weighted[-(1:11)]))
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
//'   the standard \emph{R} function \code{filter(x=vectorv, filter=weightv, 
//'   method="convolution", sides=1)}, but it's about \emph{6} times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vectorv <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' weightv <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weighted <- HighFreq::roll_cum(vectorv=vectorv, weightv=rev(weightv))
//' # compare with original
//' all.equal(vectorv, as.numeric(weighted))
//' # Second example
//' # create exponentially decaying weights
//' weightv <- exp(-0.2*1:11)
//' weightv <- weightv/sum(weightv)
//' # calculate rolling weighted sum
//' weighted <- HighFreq::roll_cum(vectorv=vectorv, weightv=rev(weightv))
//' # calculate rolling weighted sum using filter()
//' filtered <- filter(x=vectorv, filter=weightv, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filtered[-(1:11)]), as.numeric(weighted[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::ivec roll_cum(const arma::ivec& vectorv, const arma::sword& maxv) {
  uword nrows = vectorv.n_elem;
  arma::ivec roll_sum(nrows);
  
  // startup period
  roll_sum[0] = vectorv[0];
  // remaining periods
  for (uword it = 1; it < nrows; it++) {
    roll_sum[it] = roll_sum[it-1] + vectorv[it];
    if (roll_sum[it] > maxv)
      roll_sum[it] = maxv;
    if (roll_sum[it] < (-maxv))
      roll_sum[it] = -maxv;
  }  // end for
  
  return roll_sum;
}  // end roll_cum


