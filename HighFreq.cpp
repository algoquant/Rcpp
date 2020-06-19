// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////////
// Rcpp and RcppArmadillo functions for package HighFreq
////////////////////////////


////////////////////////////
// Functions for statistics
////////////////////////////


// The function vari_ance() calculates the variance of a vector using Rcpp.
//' @export
// [[Rcpp::export]]
double vari_ance(NumericVector vec_tor) {
  return sum(pow(vec_tor - sum(vec_tor)/vec_tor.size(), 2))/(vec_tor.size()-1);
}  // end vari_ance

// double vari_ance(NumericVector vec_tor);



////////////////////////////
// Functions for rolling statistics
////////////////////////////

// The function roll_sum() calculates the rolling sum over a vector using Rcpp.
//' @export
// [[Rcpp::export]]
NumericVector roll_sum(NumericVector vec_tor, int look_back) {
  int len_gth = vec_tor.size();
  NumericVector rolling_sum(len_gth);
  
  rolling_sum[0] = vec_tor[0];
  for (int it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it];
  }  // end for
  
  for (int it = look_back; it < len_gth; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it] - vec_tor[it-look_back];
  }  // end for
  
  return rolling_sum;
}  // end roll_sum



// The function roll_wsum() calculates the weighted rolling sum over a vector
// using Rcpp.
//' @export
// [[Rcpp::export]]
NumericVector roll_wsum(NumericVector vec_tor, NumericVector wei_ghts) {
  int len_gth = vec_tor.size();
  int look_back = wei_ghts.size();
  NumericVector rolling_sum(len_gth);
  
  rolling_sum[look_back-1] = sum(wei_ghts * vec_tor[Range(0, (look_back-1))]);
  for (int i = look_back; i < len_gth; ++i) {
    //    vec_tor[i] = the_ta*(eq_price - rolling_sum[i-1]) + vol_at*r_norm[i-1];
    rolling_sum[i] = sum(wei_ghts * vec_tor[Range(i-look_back+1, i)]);
  }
  
  return rolling_sum;
}


//' Calculate a time series of variance estimates over a rolling look-back interval
//' for an \emph{OHLC} time series of prices, using different range estimators
//' for variance.
//' 
//' Currently only works for vectors
//'
//' @param oh_lc An \emph{OHLC} time series of prices in \emph{xts} format.
//' @param calc_method \emph{character} string representing method for estimating
//'   variance.  The methods include:
//'   \itemize{
//'     \item "close" close to close,
//'     \item "garman_klass" Garman-Klass,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "rogers_satchell" Rogers-Satchell,
//'     \item "yang_zhang" Yang-Zhang,
//'    }
//'    (default is \code{"yang_zhang"})
//' @param look_back The size of the look-back interval, equal to the number of rows
//'   of data used for calculating the variance.
//' @param sca_le \emph{Boolean} argument: should the returns be divided by the
//'   number of seconds in each period? (default is \code{TRUE})
//'
//' @return An \emph{xts} time series with a single column and the same number of
//'   rows as the argument \code{oh_lc}.
//'
//' @details The function \code{roll_var()} calculates a time series of rolling 
//'   variance variance estimates of percentage returns, from over a
//'   \emph{vector} of returns, using several different variance estimation
//'   methods based on the range of \emph{OHLC} prices.
//'
//'   If \code{sca_le} is \code{TRUE} (the default), then the variance is divided
//'   by the squared differences of the time index (which scales the variance to
//'   units of variance per second squared.) This is useful for example, when
//'   calculating intra-day variance from minutely bar data, because dividing
//'   returns by the number of seconds decreases the effect of overnight price
//'   jumps.
//'
//'   If \code{sca_le} is \code{TRUE} (the default), then the variance is
//'   expressed in the scale of the time index of the \emph{OHLC} time series.
//'   For example, if the time index is in seconds, then the variance is given in
//'   units of variance per second squared.  If the time index is in days, then
//'   the variance is equal to the variance per day squared.
//'
//'   The time index of the \code{oh_lc} time series is assumed to be in
//'   \emph{POSIXct} format, so that its internal value is equal to the number of
//'   seconds that have elapsed since the \emph{epoch}.
//'
//'   The methods \code{"close"}, \code{"garman_klass_yz"}, and
//'   \code{"yang_zhang"} do account for close-to-open price jumps, while the
//'   methods \code{"garman_klass"} and \code{"rogers_satchell"} do not account
//'   for close-to-open price jumps.
//'
//'   The default method is \code{"yang_zhang"}, which theoretically has the
//'   lowest standard error among unbiased estimators.
//'
//'   The function \code{roll_var()} performs the same calculations as the
//'   function \code{volatility()} from package
//'   \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR}, but
//'   it's a little faster because it uses function RcppRoll::roll_sd(), and it
//'   performs less data validation.
//'
//' @examples
//' \dontrun{
//' # create minutely OHLC time series of random prices
//' oh_lc <- HighFreq::random_ohlc()
//' # calculate variance estimates for oh_lc over a 21 period interval
//' var_rolling <- HighFreq::roll_var(oh_lc, look_back=21)
//' # calculate variance estimates for SPY
//' var_rolling <- HighFreq::roll_var(HighFreq::SPY, calc_method="yang_zhang")
//' # calculate SPY variance without accounting for overnight jumps
//' var_rolling <- HighFreq::roll_var(HighFreq::SPY, calc_method="rogers_satchell")
//' }

//' @export
// [[Rcpp::export]]

NumericVector roll_var(NumericVector vec_tor, int look_back) {
  int len_gth = vec_tor.size();
  NumericVector var_vec(len_gth);
  NumericVector roll_mean(len_gth);
  
  var_vec[0] = 0;
  for (int it = 1; it < vec_tor.size(); it++) {
    var_vec[it] = vari_ance(vec_tor[Range(std::max(it-look_back+1, 0), it)]);
  }  // end for
  
  return var_vec;
}  // end roll_var



////////////////////////////
// Functions for matrix algebra
////////////////////////////


// The function inv_reg() calculates the regularized inverse 
// of the covariance matrix, by truncating the number of 
// eigen-vectors to max_eigen.
//' @export
// [[Rcpp::export]]
arma::mat inv_reg(const arma::mat& re_turns, const arma::uword& max_eigen) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  
  arma::eig_sym(eigen_val, eigen_vec, cov(re_turns));
  eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
  eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
  // arma::mat eigen_valmat = diagmat(eigen_val);
  
  return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  
}  // end inv_reg


////////////////////////////
// Functions for simulation
////////////////////////////



// The function garch_proc() simulates a GARCH model using Rcpp.
//' @export
// [[Rcpp::export]]
NumericMatrix garch_proc(int len_gth, 
                         double om_ega, 
                         double al_pha, 
                         double be_ta, 
                         NumericVector r_norm) {
  NumericVector vari_ance(len_gth);
  NumericVector re_turns(len_gth);
  vari_ance[0] = om_ega/(1-al_pha-be_ta);
  re_turns[0] = sqrt(vari_ance[0])*r_norm[0];
  
  for (int i = 1; i < len_gth; i++) {
    re_turns[i] = sqrt(vari_ance[i-1])*r_norm[i];
    vari_ance[i] = om_ega + al_pha*pow(re_turns[i], 2) + be_ta*vari_ance[i-1];
  }
  return cbind(re_turns, vari_ance);
}  // end garch_proc



// The function rcpp_ou_proc() simulates an Ornstein-Uhlenbeck process using Rcpp.
//' @export
// [[Rcpp::export]]
NumericVector rcpp_ou_proc(int len_gth, double eq_price, double vol_at, double the_ta, NumericVector r_norm) {
  NumericVector price_s(len_gth);
  NumericVector re_turns(len_gth);
  price_s[0] = eq_price;
  for (int i = 1; i < len_gth; ++i) {
    re_turns[i] = the_ta*(eq_price - price_s[i-1]) + vol_at*r_norm[i-1];
    price_s[i] = price_s[i-1] * exp(re_turns[i]);
  }
  return price_s;
}


////////////////////////////
// Functions for backtests
////////////////////////////


// The function sharpe_weights_reg() calculates the maximum 
// Sharpe ratio portfolio weights for the matrix re_turns.
// It uses the regularized inverse of the covariance matrix.
//' @export
// [[Rcpp::export]]
arma::vec sharpe_weights_reg(const arma::mat& re_turns, 
                             const arma::vec alpha_s, 
                             const arma::vec alphas_b, 
                             const arma::uword& max_eigen) {
  arma::mat in_verse = inv_reg(re_turns, max_eigen);
  arma::vec weight_s = arma::trans(arma::mean(re_turns, 0));
  arma::vec mean_s(weight_s.n_elem);
  mean_s.fill(arma::mean(weight_s));
  
  // shrink weight_s to the mean of weight_s
  weight_s = (alphas_b % weight_s + alpha_s % mean_s);
  // apply regularized inverse
  weight_s = in_verse*weight_s;
  return weight_s/sqrt(sum(square(weight_s)));
}  // end sharpe_weights_reg



// The function roll_portf() performs a loop over the 
// end_points, subsets the re_turns matrix, and calculates 
// the PCA variances using eigen decomposition.
//' @export
// [[Rcpp::export]]
arma::mat roll_portf(const arma::mat& ex_cess, // portfolio returns
                     const arma::mat& re_turns, // portfolio returns
                     const arma::uvec& start_points, 
                     const arma::uvec& end_points, 
                     const double& al_pha, 
                     const arma::uword& max_eigen) {
  arma::vec sre_turns = zeros(re_turns.n_rows);
  arma::vec weight_s(re_turns.n_cols);
  arma::vec alpha_s(re_turns.n_cols);
  alpha_s.fill(al_pha);
  arma::vec alphas_b(re_turns.n_cols);
  alphas_b.fill(1-al_pha);
  
  // perform a loop over the end_points
  for (arma::uword i = 1; i < end_points.size(); i++) {
    // subset the returns
    arma::mat sub_returns = ex_cess.rows(start_points[i-1], end_points[i-1]);
    // calculate portfolio weights
    weight_s = sharpe_weights_reg(sub_returns, alpha_s, alphas_b, max_eigen);
    // sub_returns = re_turns.rows(end_points[i-1]+1, end_points[i]);
    sre_turns.subvec(end_points[i-1]+1, end_points[i]) = re_turns.rows(end_points[i-1]+1, end_points[i])*weight_s;
    // arma::mat foo = re_turns.rows(end_points[i-1]+1, end_points[i])*weight_s;
  }  // end for
  // return the strategy returns
  return sre_turns;
}  // end roll_portf

