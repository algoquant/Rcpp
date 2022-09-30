////////////////////////////
// Function roll_reg() for performing rolling PCA and robust regressions
// Add option to perform rolling predictions using the regressions.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/roll_reg.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @export
// [[Rcpp::export]]
arma::uvec calc_endpoints(arma::uword length, arma::uword step = 1, arma::uword stub = 0) {
  
  arma::uword extra = length % step;
  arma::uvec endp;
  
  if ((stub == 0) & (extra == 0)) {
    // No stub interval
    endp = arma::regspace<uvec>(step, step, length);
  } else if ((stub == 0) & (extra > 0)) {
    // Add stub interval at end
    endp = arma::regspace<uvec>(step, step, length + step);
    endp.back() = length;
  } else if ((stub > 0) & (extra == 0)) {
    // Add initial stub interval equal to stub
    endp = arma::regspace<uvec>(stub, step, length + step);
    endp.back() = length;
  } else if ((stub > 0) & (extra > 0) & (stub == extra)) {
    // Add initial stub interval equal to stub without stub at end
    endp = arma::regspace<uvec>(stub, step, length);
  } else {
    // Add initial stub interval equal to stub and with extra stub at end
    endp = arma::regspace<uvec>(stub, step, length + step);
    endp.back() = length;
  }  // end if
  
  // Subtract 1 from endp because C++ indexing starts at 0
  endp = endp - 1;
  return endp;
  
}  // end calc_endpoints


//' @export
// [[Rcpp::export]]
arma::uvec calc_startpoints(arma::uvec endp, arma::uword look_back) {
  
  arma::uword num_points = endp.n_elem;
  arma::uvec startp = arma::join_cols(arma::zeros<uvec>(look_back), 
                                       endp.subvec(0, num_points - look_back - 1) + 1);
  
  return startp;
  
}  // end calc_startpoints


// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword dimax = 0) {
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD
  arma::svd(svdu, svdval, svdv, tseries);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > eigen_thresh*arma::sum(svdval));
  
  if (dimax == 0) {
    // Set dimax
    dimax = svdnum - 1;
  } else {
    // Adjust dimax
    dimax = std::min(dimax - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, dimax);
  svdu = svdu.cols(0, dimax);
  svdv = svdv.cols(0, dimax);
  
  // Calculate the regularized inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_inv



// Performs the same regression calculations as the
// function lm() from package stats and returns a list.
//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(arma::vec response, arma::mat design) {
  
  // Add column for intercept to explanatory matrix
  int nrows = design.n_rows;
  arma::mat designp = join_rows(ones(nrows), design);
  int ncols = designp.n_cols;
  int deg_free = (nrows - ncols);
  
  // Calculate alpha and beta coefficients for the model response ~ design
  arma::colvec coeff = arma::solve(designp, response);
  // Calculate residuals
  arma::colvec residuals = response - designp*coeff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (nrows-1)*arma::var(response);
  double res_sumsq = arma::dot(residuals, residuals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate R-squared and F-statistic
  double rsquared = exp_sumsq/tot_sumsq;
  double fstat = (exp_sumsq*deg_free)/(res_sumsq*(ncols-1));
  // arma::rowvec stats=join_horiz(rsquared, fstat);
  Rcpp::NumericVector stats(2);
  stats(0) = rsquared;
  stats(1) = fstat;
  stats.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // Calculate standard errors of beta coefficients
  arma::colvec stderr = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(designp)*designp)));
  // Calculate t-values and p-values of beta coefficients
  arma::colvec t_vals = coeff/stderr;
  arma::colvec pvals = 2*Rcpp::pt(-abs(wrap(t_vals)), deg_free);
  NumericMatrix coeffmat = wrap(join_rows(join_rows(join_rows(coeff, stderr), t_vals), pvals));
  Rcpp::colnames(coeffmat) = Rcpp::CharacterVector::create("coeff", "stderr", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeffmat,
                            // Named("residuals") = residuals,
                            Named("z_score") = residuals(nrows-1)/arma::stddev(residuals),
                            Named("stats") = stats);
  
}  // end calc_lm



// Extract regression elements into a vector.
//' @export
// [[Rcpp::export]]
arma::vec calc_lm_vec(arma::vec response, arma::mat design) {
  
  Rcpp::List reg_list = calc_lm(response, design);
  arma::mat coeff = reg_list["coefficients"];
  arma::vec z_score = reg_list["z_score"];
  // arma::vec reg_data = coeff.col(0);
  // arma::vec reg_data = join_rows(coeff.col(0), coeff.col(2));
  // arma::vec reg_data = coeff.as_col();
  // uvec col_s = {0, 2};
  arma::vec reg_data = coeff.cols(uvec {0, 2}).as_col();
  reg_data = join_cols(reg_data, z_score);
  
  return reg_data;
  
}  // end calc_lm_vec



////////////////////////////////////////////////////////////
// Define C++ enum type for different methods of regularization,
// methodsfor calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum methodenum {moment, least_squares, quantile, nonparametric, regular, sharpem, 
                 maxsharpe, maxsharpemed, minvarlin, minvarquad, kellym, robustm, 
                 sumsq, sumone, voltarget, voleqw};

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
  else if (method == "sharpem")
    return methodenum::sharpem;
  else if (method == "maxsharpe")
    return methodenum::maxsharpe;
  else if (method == "maxsharpemed")
    return methodenum::maxsharpemed;
  else if (method == "minvarlin")
    return methodenum::minvarlin;
  else if (method == "minvarquad")
    return methodenum::minvarquad;
  else if (method == "kellym")
    return methodenum::kellym;
  else if (method == "robustm")
    return methodenum::robustm;
  else if (method == "sumsq")
    return methodenum::sumsq;
  else if (method == "voltarget")
    return methodenum::voltarget;
  else 
    return methodenum::moment;
}  // end calc_method



////////////////////////////////////////////////////////////
//' Perform multivariate regression using different methods, and return a vector
//' of regression coefficients, their t-values, and the last residual z-score.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predictor} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @param \code{intercept} A \emph{Boolean} specifying whether an intercept
//'   term should be added to the predictor (the default is \code{intercept =
//'   TRUE}).
//'
//' @param \code{method} A \emph{character string} specifying the type of the
//'   regression model the default is \code{method = "least_squares"} - see
//'   Details).
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{predictor} matrix (the default is \code{1e-5}).
//'   
//' @param \code{dimax} An \emph{integer} equal to the number of singular
//'   values used for calculating the regularized inverse of the \code{predictor}
//'   matrix (the default is \code{0} - equivalent to \code{dimax} equal to
//'   the number of columns of \code{predictor}).
//'   
//' @param \code{confl} The confidence level for calculating the
//'   quantiles of returns (the default is \code{confl = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @return A single-row matrix with
//' A vector with the regression coefficients, their t-values, and the
//'   last residual z-score.
//'
//' @details
//'   The function \code{calc_reg()} performs multivariate regression using
//'   different methods, and returns a vector of regression coefficients, their
//'   t-values, and the last residual z-score.
//'
//'   If \code{method = "least_squares"} (the default) then it performs the
//'   standard least squares regression, the same as the function
//'   \code{calc_lm()}, and the function \code{lm()} from the \code{R} package
//'   \emph{stats}.
//'   But it uses \code{RcppArmadillo} \code{C++} code so it's several times
//'   faster than \code{lm()}.
//'
//'   If \code{method = "regular"} then it performs shrinkage regression.  It
//'   calculates the regularized inverse of the \code{predictor} matrix from its
//'   singular value decomposition.  It performs regularization by selecting
//'   only the largest singular values equal in number to \code{dimax}.
//'   
//'   If \code{method = "quantile"} then it performs quantile regression (not
//'   implemented yet).
//' 
//'   If \code{intercept = TRUE} then an extra intercept column (unit column) is
//'   added to the predictor matrix (the default is \code{intercept = FALSE}).
//'   
//'   The length of the return vector depends on the number of columns of the
//'   \code{predictor} matrix (including the intercept column, if it's added).
//'   The length of the return vector is equal to the number of regression
//'   coefficients, plus their t-values, plus the z-score.
//'   The number of regression coefficients is equal to the number of columns of
//'   the \code{predictor} matrix (including the intercept column, if it's
//'   added).
//'   The number of t-values is equal to the number of coefficients.
//' 
//'   For example, if the number of columns of the \code{predictor} matrix is
//'   equal to \code{n}, and if \code{intercept = TRUE} (the default), then
//'   \code{calc_reg()} returns a vector with \code{2n+3} elements: \code{n+1}
//'   regression coefficients (including the intercept coefficient), \code{n+1}
//'   corresponding t-values, and \code{1} z-score value.
//'
//'   If \code{intercept = FALSE}, then \code{calc_reg()} returns a vector with
//'   \code{2n+1} elements: \code{n} regression coefficients (without the
//'   intercept coefficient), \code{n} corresponding t-values, and \code{1}
//'   z-score value.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' response <- returns[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predictor <- returns[, -1]
//' # Perform multivariate regression using lm()
//' lmod <- lm(response ~ predictor)
//' lmodsum <- summary(lmod)
//' coeff <- lmodsum$coefficients
//' # Perform multivariate regression using calc_reg()
//' reg_arma <- drop(HighFreq::calc_reg(response=response, predictor=predictor))
//' # Compare the outputs of both functions
//' all.equal(reg_arma[1:(2*(1+NCOL(predictor)))], 
//'   c(coeff[, "Estimate"], coeff[, "t value"]), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_reg(response=response, predictor=predictor),
//'   Rcode=lm(response ~ predictor),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_reg(const arma::mat& response, 
                   const arma::mat& predictor,
                   bool intercept = true,
                   std::string method = "least_squares",
                   double eigen_thresh = 1e-5, // Threshold level for discarding small singular values
                   arma::uword dimax = 0, // Regularization intensity
                   double confl = 0.1, // Confidence level for calculating the quantiles of returns
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
  
  // Apply different calculation methods for weights
  switch(calc_method(method)) {
  case methodenum::least_squares: {
    // Calculate regression coefficients for the model response ~ predictor
    coeff = arma::solve(predictori, response);
    break;
  }  // end least_squares
  case methodenum::regular: {
    // Calculate shrinkage regression coefficients
    coeff = calc_inv(predictori, eigen_thresh, dimax)*response;
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
  arma::mat stderrv = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(predictori)*predictori)));
  // Calculate t-values of the beta coefficients
  tvals = coeff/stderrv;
  
  // Calculate z-score
  mat zscore = residuals(nrows-1, 0)/arma::stddev(residuals);
  
  // Combine regression data
  reg_data.rows(0, ncols-1) = coeff;
  reg_data.rows(ncols, 2*ncols-1) = tvals;
  reg_data.row(2*ncols) = zscore;
  
  return reg_data.t();
  
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
//' @param \code{endp} An \emph{integer} vector of end points.
//' 
//' @param \code{startp} An \emph{integer} vector of start points.
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
//'   \code{returns}).
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
//' returns <- na.omit(rutils::etfenv$returns[, c("IEF", "VTI", "XLF")])
//' # Response equals IEF returns
//' response <- returns[, 1]
//' # Design matrix equals VTI and XLF returns
//' design <- returns[, -1]
//' # Calculate Z-scores from rolling time series regression using RcppArmadillo
//' look_back <- 11
//' z_scores <- HighFreq::roll_reg(response=response, design=design, look_back=look_back)
//' # Calculate z-scores in R from rolling multivariate regression using lm()
//' z_scoresr <- sapply(1:NROW(design), function(ro_w) {
//'   if (ro_w == 1) return(0)
//'   startpoint <- max(1, ro_w-look_back+1)
//'   sub_response <- response[startpoint:ro_w]
//'   sub_design <- design[startpoint:ro_w, ]
//'   reg_model <- lm(sub_response ~ sub_design)
//'   residuals <- reg_model$residuals
//'   residuals[NROW(residuals)]/sd(residuals)
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
                   arma::uvec endp = 0, 
                   arma::uvec startp = 0, 
                   arma::uword step = 1, 
                   arma::uword look_back = 0,
                   std::string method = "least_squares",
                   double eigen_thresh = 0.001,
                   int eigen_max = 0,
                   double confi_level = 0.1,
                   double alpha = 0.0) {
  
  arma::uword nrows = design.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate regression matrix
  arma::vec sub_response;
  arma::mat sub_design;
  arma::colvec reg_data;
  arma::uword ncols = design.n_cols;
  arma::uword num_points = endpts.n_elem;
  arma::mat reg_stats(num_points, (2*(ncols + 1) + 1), fill::zeros);
  
  
  // Perform loop over the endpts
  for (arma::uword ep = 0; ep < num_points; ep++) {
    // Calculate regression coefficients
    if (endpts(ep) > startpts(ep)) {
      sub_response = response.subvec(startpts(ep), endpts(ep));
      sub_design = design.rows(startpts(ep), endpts(ep));
      reg_data = calc_reg(sub_response, sub_design, method=method, eigen_thresh=eigen_thresh, eigen_max=eigen_max, confi_level=confi_level, alpha=alpha);
      reg_stats.row(ep) = conv_to< rowvec >::from(reg_data);
    }  // end if
  }  // end for
  
  // Warmup period
  // reg_stats.rows(0, ncols+1) = zeros(ncols+2, (ncols + 1));
  // for (arma::uword it = (ncols+2); it < look_back; it++) {
  //   sub_response = response.subvec(0, it);
  //   sub_design = design.rows(0, it);
  //   reg_data = calc_reg(sub_response, sub_design);
  //   reg_stats.row(it) = conv_to< rowvec >::from(reg_data);
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = look_back; it < nrows; it++) {
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
arma::mat roll_wsum_mat(arma::mat vectorv, arma::mat weights) {
  uword nrows = vectorv.n_rows;
  uword look_back = weights.n_rows;
  // arma::mat vectorv(nrows);
  // arma::vec rev_weights = arma::reverse(weights);
  // arma::vec rev_weights = weights;
  
  // warmup period
  // vectorv.n_rows(0, look_back-2) = vectorv.subvec(0, look_back-2);
  
  // remaining periods
  for (uword it = look_back; it < nrows; it++) {
    vectorv.row(it-1) = trans(weights)*vectorv.rows(it-look_back, it-1);
  }  // end for
  
  return vectorv;
}  // end roll_wsum_mat


