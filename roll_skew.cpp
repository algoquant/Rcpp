////////////////////////////
// Function roll_skew() for calculating the rolling skewness estimator 
// (third moment) using RcppArmadillo.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/roll_skew.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


////////////////////////////////////////////////////////////
// Switch statement in calc_skew() uses C++ enum type.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
enum skew_type {Pearson, Quantile, Nonparametric};
// Map string to C++ enum type for switch statement.
skew_type calc_skew_type(const std::string method) {
  if (method == "Pearson" || method == "pearson" || method == "p") 
    return skew_type::Pearson;
  else if (method == "Quantile" || method == "quantile" || method == "q")
    return skew_type::Quantile;
  else if (method == "Nonparametric" || method == "nonparametric" || method == "n")
    return skew_type::Nonparametric;
  else 
    return skew_type::Pearson;
}  // end calc_skew_type



////////////////////////////////////////////////////////////
//' Calculate the skewness of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{se_ries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{method} A \emph{string} specifying the type of skewness (see
//'   Details). (The default is the \code{method = "pearson"}.)
//'
//' @param \code{alpha} The confidence level for calculating the quantiles.
//'   (the default is \code{alpha = 0.25}).
//'
//' @return A single-row matrix with the skewness of the columns of
//'   \code{se_ries}.
//'
//' @details The function \code{calc_skew()} calculates the skewness of the
//'   columns of a \emph{time series} or a \emph{matrix} of data using
//'   \code{RcppArmadillo} \code{C++} code.
//'
//'   If \code{method = "pearson"} (the default) then \code{calc_skew()}
//'   calculates the Pearson skewness using the third moment of the data.
//'
//'   If \code{method = "quantile"} then it calculates the skewness using the
//'   differences between the quantiles of the data.
//'
//'   If \code{method = "nonparametric"} then it calculates the skewness as the
//'   difference between the mean of the data minus its median, divided by the
//'   standard deviation.
//'   
//'   The code examples below compare the function \code{calc_skew()} with the
//'   skewness calculated using \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Calculate VTI returns
//' returns <- na.omit(rutils::etfenv$returns[, "VTI", drop=FALSE])
//' # Calculate the Pearson skewness
//' HighFreq::calc_skew(returns)
//' # Compare HighFreq::calc_skew() with Pearson skewness
//' calc_skewr <- function(x) {
//'   x <- (x-mean(x)); nr <- NROW(x);
//'   nr*sum(x^3)/(var(x))^1.5/(nr-1)/(nr-2)
//' }  # end calc_skewr
//' all.equal(HighFreq::calc_skew(returns), 
//'   calc_skewr(returns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(returns),
//'   Rcode=calc_skewr(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile skewness
//' HighFreq::calc_skew(returns, method = "quantile", alpha = 0.1)
//' # Compare HighFreq::calc_skew() with quantile skewness
//' calc_skewq <- function(x) {
//'   	quantiles <- quantile(x, c(0.25, 0.5, 0.75), type=5)
//'   	(quantiles[3] + quantiles[1] - 2*quantiles[2])/(quantiles[3] - quantiles[1])
//' }  # end calc_skewq
//' all.equal(drop(HighFreq::calc_skew(returns, method = "quantile")), 
//'   calc_skewq(returns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(returns, method = "quantile"),
//'   Rcode=calc_skewq(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric skewness
//' HighFreq::calc_skew(returns, method = "nonparametric")
//' # Compare HighFreq::calc_skew() with R nonparametric skewness
//' all.equal(drop(HighFreq::calc_skew(returns, method = "nonparametric")), 
//'   (mean(returns)-median(returns))/sd(returns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(returns, method = "nonparametric"),
//'   Rcode=(mean(returns)-median(returns))/sd(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_skew(arma::mat se_ries,
                    std::string method = "pearson", 
                    double alpha = 0.25) {
  
  // Switch statement for all the different methods of skewness
  switch(calc_skew_type(method)) {
  case skew_type::Pearson: {  // Pearson
    double nrows = se_ries.n_rows;
    arma::mat mean_s = arma::mean(se_ries);
    arma::mat var_s = arma::var(se_ries);
    // De-mean the columns of se_ries
    se_ries.each_row() -= mean_s;
    return (nrows/(nrows-1)/(nrows-2))*arma::sum(arma::pow(se_ries, 3))/arma::pow(var_s, 1.5);
  }  // end pearson
  case skew_type::Quantile: {  // Quantile
    arma::vec probs = {alpha, 0.5, 1.0 - alpha};
    arma::mat quantiles = quantile(se_ries, probs);
    return (quantiles.row(2) + quantiles.row(0) - 2*quantiles.row(1))/(quantiles.row(2) - quantiles.row(0));
  }  // end quantile
  case skew_type::Nonparametric: {  // Nonparametric
    return (arma::mean(se_ries) - arma::median(se_ries))/arma::stddev(se_ries);
  }  // end nonparametric
  default : {
    cout << "Invalid method" << endl;
    return arma::zeros<rowvec>(se_ries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_skew



////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of skewness estimates over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{se_ries} A \emph{time series} or a \emph{matrix} of data.
//'    
//' @param \code{ste_p} The number of time periods between the end points.
//'
//' @param \code{look_back} The number of end points in the look-back interval.
//'
//' @param \code{method} A \emph{string} specifying the type of skewness.  (The
//'   default is the \code{method = "pearson"}.)
//'
//' @param \code{alpha} The confidence level for calculating the quantiles.
//'   (the default is \code{alpha = 0.25}).
//'
//' @return A \emph{matrix} with the same number of columns as the input time
//'   series \code{se_ries}, and the number of rows equal to the number of end
//'   points.
//'   
//' @details 
//'   The function \code{roll_skew()} calculates a \emph{matrix} of skewness
//'   estimates over rolling look-back intervals attached at the end points of
//'   the \emph{time series} \code{se_ries}.
//'   
//'   It first calculates a vector of end points separated by \code{ste_p} time
//'   periods. It calculates the end points along the rows of \code{se_ries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{ste_p} time periods.
//'   
//'   It then performs a loop over the end points, and at each end point it
//'   subsets the time series \code{se_ries} over a look-back interval equal
//'   to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_skew()}, which
//'   calculates the skewness.
//'   See the function \code{calc_skew()} for a description of the skewness
//'   methods.
//'   
//'   For example, the rolling skewness at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{ste_p = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_skew()} with the parameter \code{ste_p = 1}
//'   performs the same calculation as the function \code{roll_skew()} from
//'   package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo}
//'   \code{C++} code.
//'
//'   The function \code{roll_skew()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the rolling skewness at 25 day end points, with a 75 day look-back
//' variance <- HighFreq::roll_skew(returns, ste_p=25, look_back=3)
//' # Compare the skewness estimates over 11-period lookback intervals
//' all.equal(HighFreq::roll_skew(returns, look_back=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_skew(returns, n=11)), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with RcppRoll
//' library(microbenchmark)
//' summary(microbenchmark(
//'   RcppArmadillo=HighFreq::roll_skew(returns, look_back=11),
//'   RcppRoll=RcppRoll::roll_skew(returns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_skew(arma::mat se_ries, 
                    arma::uword ste_p = 1, 
                    arma::uword look_back = 11, 
                    std::string method = "pearson", 
                    double alpha = 0.25) {
  
  // Calculate end points
  arma::uword nrows = se_ries.n_rows;
  arma::uvec endp = calc_endpoints(nrows, ste_p);
  // Start points equal to end points lagged by look_back
  arma::uvec startp = calc_startpoints(endp, look_back);
  // Allocate skewness matrix
  arma::uword num_points = endp.n_elem;
  arma::mat skew_ness = arma::zeros<mat>(num_points, se_ries.n_cols);
  
  // Perform loop over the endp
  for (arma::uword ep = 0; ep < num_points; ep++) {
    // Calculate skewness
    if (endp(ep) > startp(ep)) {
      skew_ness.row(ep) = calc_skew(se_ries.rows(startp(ep), endp(ep)), method);
    }  // end if
  }  // end for
  
  return skew_ness;
  
}  // end roll_skew



//' @export
// [[Rcpp::export]]
arma::mat roll_skewo(arma::mat se_ries,
                     std::string method = "pearson", 
                     double alpha = 0.25, 
                     const arma::uword look_back = 11) {
  
  arma::uword nrows = se_ries.n_rows;
  arma::uword ncols = se_ries.n_cols;
  arma::mat skew_ness = arma::zeros<mat>(nrows, ncols);
  arma::mat sub_series;
  
  // Warmup period
  // skew_ness.rows(0, ncols+1) = arma::zeros(ncols+2, (ncols + 1));
  for (arma::uword it = 3; it < look_back; it++) {
    sub_series = se_ries.rows(0, it);
    skew_ness.row(it) = calc_skew(sub_series, method, alpha);
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < nrows; it++) {
    sub_series = se_ries.rows(it-look_back+1, it);
    skew_ness.row(it) = calc_skew(sub_series, method, alpha);
  }  // end for
  
  return skew_ness;
  
}  // end roll_skewo

