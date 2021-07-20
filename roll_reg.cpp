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
arma::uvec calc_endpoints(arma::uword length, arma::uword step = 1, arma::uword stub = 0) {
  
  arma::uword extra = length % step;
  arma::uvec end_p;
  
  if ((stub == 0) & (extra == 0)) {
    // No stub interval
    end_p = arma::regspace<uvec>(step, step, length);
  } else if ((stub == 0) & (extra > 0)) {
    // Add stub interval at end
    end_p = arma::regspace<uvec>(step, step, length + step);
    end_p.back() = length;
  } else if ((stub > 0) & (extra == 0)) {
    // Add initial stub interval equal to stub
    end_p = arma::regspace<uvec>(stub, step, length + step);
    end_p.back() = length;
  } else if ((stub > 0) & (extra > 0) & (stub == extra)) {
    // Add initial stub interval equal to stub without stub at end
    end_p = arma::regspace<uvec>(stub, step, length);
  } else {
    // Add initial stub interval equal to stub and with extra stub at end
    end_p = arma::regspace<uvec>(stub, step, length + step);
    end_p.back() = length;
  }  // end if
  
  // Subtract 1 from end_p because C++ indexing starts at 0
  end_p = end_p - 1;
  return end_p;
  
}  // end calc_endpoints


//' @export
// [[Rcpp::export]]
arma::uvec calc_startpoints(arma::uvec end_p, arma::uword look_back) {
  
  arma::uword num_points = end_p.n_elem;
  arma::uvec start_p = arma::join_cols(arma::zeros<uvec>(look_back), 
                                       end_p.subvec(0, num_points - look_back - 1) + 1);
  
  return start_p;
  
}  // end calc_startpoints




// Performs the same regression calculations as the
// function lm() from package stats and returns a list.
//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(arma::vec response, arma::mat design) {
  
  // Add column for intercept to explanatory matrix
  int num_rows = design.n_rows;
  arma::mat design_p = join_rows(ones(num_rows), design);
  int num_cols = design_p.n_cols;
  int deg_free = (num_rows - num_cols);
  
  // Calculate alpha and beta coefficients for the model response ~ design
  arma::colvec co_eff = arma::solve(design_p, response);
  // Calculate residuals
  arma::colvec resid_uals = response - design_p*co_eff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (num_rows-1)*arma::var(response);
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
  arma::colvec t_vals = co_eff/std_err;
  arma::colvec p_vals = 2*Rcpp::pt(-abs(wrap(t_vals)), deg_free);
  NumericMatrix coeff_mat = wrap(join_rows(join_rows(join_rows(co_eff, std_err), t_vals), p_vals));
  Rcpp::colnames(coeff_mat) = Rcpp::CharacterVector::create("coeff", "std_err", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeff_mat,
                            // Named("residuals") = resid_uals,
                            Named("z_score") = resid_uals(num_rows-1)/arma::stddev(resid_uals),
                            Named("stats") = stat_s);
  
}  // end calc_lm



// Extract regression elements into a vector.
//' @export
// [[Rcpp::export]]
arma::vec calc_lm_vec(arma::vec response, arma::mat design) {
  
  Rcpp::List reg_list = calc_lm(response, design);
  arma::mat co_eff = reg_list["coefficients"];
  arma::vec z_score = reg_list["z_score"];
  // arma::vec reg_data = co_eff.col(0);
  // arma::vec reg_data = join_rows(co_eff.col(0), co_eff.col(2));
  // arma::vec reg_data = co_eff.as_col();
  // uvec col_s = {0, 2};
  arma::vec reg_data = co_eff.cols(uvec {0, 2}).as_col();
  reg_data = join_cols(reg_data, z_score);
  
  return reg_data;
  
}  // end calc_lm_vec



////////////////////////////////////////////////////////////
// Define C++ enum type for the different methods for regularization,
// calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum meth_od {moment, least_squares, quantile, nonparametric, rank_sharpe, 
              max_sharpe, max_sharpe_median, min_var, min_varpca, rank, rankrob};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
meth_od calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return meth_od::moment;
  else if (method == "least_squares" || method == "l")
    return meth_od::least_squares;
  else if (method == "quantile" || method == "q")
    return meth_od::quantile;
  else if (method == "nonparametric" || method == "n")
    return meth_od::nonparametric;
  else if (method == "rank_sharpe")
    return meth_od::rank_sharpe;
  else if (method == "max_sharpe")
    return meth_od::max_sharpe;
  else if (method == "max_sharpe_median")
    return meth_od::max_sharpe_median;
  else if (method == "min_var")
    return meth_od::min_var;
  else if (method == "min_varpca")
    return meth_od::min_varpca;
  else if (method == "rank")
    return meth_od::rank;
  else if (method == "rankrob")
    return meth_od::rankrob;
  else 
    return meth_od::moment;
}  // end calc_method




////////////////////////////////////////////////////////////
//' Perform multivariate regression using different methods, and return a vector
//' of regression coefficients, t-values, and the last residual z-score.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{design} A \emph{time series} or a \emph{matrix} of design data
//'   (predictor or explanatory data).
//' 
//' @param \code{method} A \emph{string} specifying the type of regression model
//'   (see Details).  (The default is \code{method = "least_squares"})
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small eigenvalues in order to regularize the matrix inverse.  (The default
//'   is \code{0.001})
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of
//'   eigenvectors used for calculating the regularized inverse of the
//'   covariance \emph{matrix} (the default is the number of columns of
//'   \code{re_turns}).
//'   
//' @param \code{confi_level} The confidence level for calculating the
//'   quantiles. (the default is \code{confi_level = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @return A vector of regression coefficients, t-values, and the last
//'   residual z-score.
//'   For example, if the design matrix has \code{2} columns of data, then
//'   \code{calc_reg()} returns a vector with \code{7} elements: \code{3}
//'   regression coefficients (including the intercept coefficient), \code{3}
//'   corresponding t-values, and \code{1} z-score.
//'
//' @details 
//'   The function \code{calc_reg()} performs multivariate regression using
//'   different methods, and returns a vector of regression coefficients, their
//'   t-values, and the last residual z-score.
//' 
//'   If \code{method = "least_squares"} (the default) then it performs the
//'   standard least squares regression, the same as the function
//'   \code{calc_reg()}, and the function \code{lm()} from package \emph{stats}.
//'   It uses \code{RcppArmadillo} \code{C++} code so it's several times faster
//'   than \code{lm()}.
//'
//'   If \code{method = "quantile"} then it performs quantile regression (not
//'   implemented yet).
//'
//'   \code{calc_weights()} applies dimension regularization to calculate the
//'   inverse of the covariance \emph{matrix} of returns from its eigen
//'   decomposition, using the function \code{arma::eig_sym()}.
//'   
//'   In addition, it applies shrinkage to the \emph{vector} of mean column
//'   returns, by shrinking it to its common mean value.
//'   The shrinkage intensity \code{alpha} determines the amount of shrinkage 
//'   that is applied, with \code{alpha = 0} representing no shrinkage (with 
//'   the estimator of mean returns equal to the means of the columns of 
//'   \code{re_turns}), and \code{alpha = 1} representing complete shrinkage 
//'   (with the estimator of mean returns equal to the single mean of all the
//'   columns of \code{re_turns})
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("IEF", "VTI", "XLF")])
//' # Response equals IEF returns
//' res_ponse <- re_turns[, 1]
//' # Design matrix equals VTI and XLF returns
//' de_sign <- re_turns[, -1]
//' # Perform multivariate regression using lm()
//' reg_model <- lm(res_ponse ~ de_sign)
//' sum_mary <- summary(reg_model)
//' co_eff <- sum_mary$coefficients
//' # Perform multivariate regression using calc_reg()
//' reg_arma <- drop(HighFreq::calc_reg(response=res_ponse, design=de_sign))
//' # Compare the outputs of both functions
//' all.equal(reg_arma[1:(2*(1+NCOL(de_sign)))], 
//'   c(co_eff[, "Estimate"], co_eff[, "t value"]), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_reg(response=res_ponse, design=de_sign),
//'   Rcode=lm(res_ponse ~ de_sign),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::colvec calc_reg(arma::vec response, 
                      arma::mat design,
                      std::string method = "least_squares",
                      double eigen_thresh = 0.001,
                      int eigen_max = 0,
                      double confi_level = 0.1,
                      double alpha = 0.0) {
  
  // Add column for intercept to explanatory matrix
  int num_rows = design.n_rows;
  arma::mat design_p = join_rows(ones(num_rows), design);
  int num_cols = design_p.n_cols;
  int deg_free = (num_rows - num_cols);
  arma::colvec co_eff(num_cols, fill::zeros);
  arma::colvec reg_data(2*num_cols+1, fill::zeros);
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case meth_od::least_squares: {
    // Calculate alpha and beta coefficients for the model response ~ design
    co_eff = arma::solve(design_p, response);
    break;
  }  // end least_squares
  case meth_od::quantile: {
    // Not implemented yet
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return reg_data;
  }  // end default
  }  // end switch
  
  // Calculate residuals
  arma::colvec resid_uals = response - design_p*co_eff;
  
  // Calculate TSS, RSS, and ESS
  // double tot_sumsq = (num_rows-1)*arma::var(response);
  double res_sumsq = arma::dot(resid_uals, resid_uals);
  // double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate standard errors of beta coefficients
  arma::colvec std_err = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(design_p)*design_p)));
  // Calculate t-values and p-values of beta coefficients
  arma::colvec t_vals = co_eff/std_err;
  
  // Calculate z-score
  double z_score = resid_uals(num_rows-1)/arma::stddev(resid_uals);
  
  // Combine regression data
  reg_data.subvec(0, num_cols-1) = co_eff;
  reg_data.subvec(num_cols, 2*num_cols-1) = t_vals;
  reg_data(2*num_cols) = z_score;
  
  return reg_data;
  
}  // end calc_reg




////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of regression coefficients, t-values, and z-scores
//' at the end points of the design matrix.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{design} A \emph{time series} or a \emph{matrix} of design data
//'   (predictor or explanatory data).
//'   
//' @param \code{end_p} An \emph{integer} vector of end points.
//' 
//' @param \code{start_p} An \emph{integer} vector of start points.
//' 
//' @param \code{step} The number of time periods between the end points.
//'
//' @param \code{look_back} The number of end points in the look-back interval.
//'
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points.
//' 
//' @param \code{method} A \emph{string} specifying the type of regression model
//'   (see Details).  (The default is \code{method = "least_squares"})
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small eigenvalues in order to regularize the matrix inverse.  (The default
//'   is \code{0.001})
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of
//'   eigenvectors used for calculating the regularized inverse of the
//'   covariance \emph{matrix} (the default is the number of columns of
//'   \code{re_turns}).
//'   
//' @param \code{confi_level} The confidence level for calculating the
//'   quantiles. (the default is \code{confi_level = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @return A column \emph{vector} of the same length as the number of rows of
//'   \code{design}.
//'
//' @details 
//'   The function \code{roll_reg()} calculates a \emph{matrix} of regression
//'   coefficients, t-values, and z-scores at the end points of the design
//'   matrix.
//'   
//'   It first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{design}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//'   
//'   It then performs a loop over the end points, and at each end point it
//'   subsets the time series \code{design} over a look-back interval equal
//'   to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_reg()}, which
//'   calculates the regression data.
//'   
//'   For example, the rolling regression at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("IEF", "VTI", "XLF")])
//' # Response equals IEF returns
//' res_ponse <- re_turns[, 1]
//' # Design matrix equals VTI and XLF returns
//' de_sign <- re_turns[, -1]
//' # Calculate Z-scores from rolling time series regression using RcppArmadillo
//' look_back <- 11
//' z_scores <- HighFreq::roll_reg(response=res_ponse, design=de_sign, look_back=look_back)
//' # Calculate z-scores in R from rolling multivariate regression using lm()
//' z_scoresr <- sapply(1:NROW(de_sign), function(ro_w) {
//'   if (ro_w == 1) return(0)
//'   start_point <- max(1, ro_w-look_back+1)
//'   sub_response <- res_ponse[start_point:ro_w]
//'   sub_design <- de_sign[start_point:ro_w, ]
//'   reg_model <- lm(sub_response ~ sub_design)
//'   resid_uals <- reg_model$residuals
//'   resid_uals[NROW(resid_uals)]/sd(resid_uals)
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(z_scores[-(1:look_back)], z_scoresr[-(1:look_back)], 
//'   check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_reg(arma::vec response, 
                   arma::mat design, 
                   arma::uvec end_p = 0, 
                   arma::uvec start_p = 0, 
                   arma::uword step = 1, 
                   arma::uword look_back = 0,
                   std::string method = "least_squares",
                   double eigen_thresh = 0.001,
                   int eigen_max = 0,
                   double confi_level = 0.1,
                   double alpha = 0.0) {
  
  arma::uword num_rows = design.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(end_p) == 0) {
    end_pts = calc_endpoints(num_rows, step);
  } else {
    // Copy end points
    end_pts = end_p;
  }  // end if
  
  // Calculate start points if missing
  if (sum(start_p) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = start_p;
  }  // end if
  
  
  // Allocate regression matrix
  arma::vec sub_response;
  arma::mat sub_design;
  arma::colvec reg_data;
  arma::uword num_cols = design.n_cols;
  arma::uword num_points = end_pts.n_elem;
  arma::mat reg_stats(num_points, (2*(num_cols + 1) + 1), fill::zeros);
  
  
  // Perform loop over the end_pts
  for (arma::uword ep = 0; ep < num_points; ep++) {
    // Calculate regression coefficients
    if (end_pts(ep) > start_pts(ep)) {
      sub_response = response.subvec(start_pts(ep), end_pts(ep));
      sub_design = design.rows(start_pts(ep), end_pts(ep));
      reg_data = calc_reg(sub_response, sub_design, method=method, eigen_thresh=eigen_thresh, eigen_max=eigen_max, confi_level=confi_level, alpha=alpha);
      reg_stats.row(ep) = conv_to< rowvec >::from(reg_data);
    }  // end if
  }  // end for
  
  // Warmup period
  // reg_stats.rows(0, num_cols+1) = zeros(num_cols+2, (num_cols + 1));
  // for (arma::uword it = (num_cols+2); it < look_back; it++) {
  //   sub_response = response.subvec(0, it);
  //   sub_design = design.rows(0, it);
  //   reg_data = calc_reg(sub_response, sub_design);
  //   reg_stats.row(it) = conv_to< rowvec >::from(reg_data);
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = look_back; it < num_rows; it++) {
  //   sub_response = response.subvec(it-look_back+1, it);
  //   sub_design = design.rows(it-look_back+1, it);
  //   reg_data = calc_reg(sub_response, sub_design, method=method, eigen_thresh=eigen_thresh, eigen_max=eigen_max, confi_level=confi_level, alpha=alpha);
  //   reg_stats.row(it) = conv_to< rowvec >::from(reg_data);
  // }  // end for
  
  return reg_stats;
  
}  // end roll_reg




////////////////////////////////////////////////////////////
//' The function roll_wsum_mat() calculates the rolling weighted sum 
//' over a vector of data using \emph{RcppArmadillo}.
//' @export
// [[Rcpp::export]]
arma::mat roll_wsum_mat(arma::mat vec_tor, arma::mat weight_s) {
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


