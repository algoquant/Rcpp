////////////////////////////
// Function roll_kurtosis() for calculating the rolling kurtosis estimator 
// (fourth moment) using RcppArmadillo.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/roll_kurtosis.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


////////////////////////////////////////////////////////////
// Switch statement in calc_kurtosis() uses C++ enum type.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
enum kurtosis_type {Pearson, Quantile, Nonparametric};
// Map string to C++ enum type for switch statement.
kurtosis_type calc_kurtosis_type(const std::string& method) {
  if (method == "Pearson" || method == "pearson" || method == "p") 
    return kurtosis_type::Pearson;
  else if (method == "Quantile" || method == "quantile" || method == "q")
    return kurtosis_type::Quantile;
  else if (method == "Nonparametric" || method == "nonparametric" || method == "n")
    return kurtosis_type::Nonparametric;
  else 
    return kurtosis_type::Pearson;
}  // end calc_kurtosis_type



////////////////////////////////////////////////////////////
//' Calculate the kurtosis of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{se_ries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{method} A \emph{string} specifying the type of kurtosis (see
//'   Details). (The default is the \code{method = "pearson"}.)
//'
//' @param \code{alpha} The confidence level for calculating the quantiles.
//'   (the default is \code{alpha = 0.25}).
//'
//' @return A single-row matrix with the kurtosis of the columns of
//'   \code{se_ries}.
//'
//' @details 
//'   The function \code{calc_kurtosis()} calculates the kurtosis of the columns of
//'   a \emph{time series} or a \emph{matrix} of data using \code{RcppArmadillo}
//'   \code{C++} code.
//'
//'   If \code{method = "pearson"} (the default) then \code{calc_kurtosis()}
//'   calculates the Pearson kurtosis using the third moment of the data.
//'
//'   If \code{method = "quantile"} then it calculates the kurtosis using the
//'   differences between the quantiles of the data.
//'
//'   If \code{method = "nonparametric"} then it calculates the kurtosis as the
//'   difference between the mean of the data minus its median, divided by the
//'   standard deviation.
//'   
//'   If the number of rows of \code{se_ries} is less than \code{3} then it
//'   returns zeros.
//'   
//'   The code examples below compare the function \code{calc_kurtosis()} with the
//'   kurtosis calculated using \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Calculate VTI returns
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the Pearson kurtosis
//' HighFreq::calc_kurtosis(returns)
//' # Compare HighFreq::calc_kurtosis() with Pearson kurtosis
//' calc_kurtosisr <- function(x) {
//'   x <- (x-mean(x)); nr <- NROW(x);
//'   nr*sum(x^3)/(var(x))^1.5/(nr-1)/(nr-2)
//' }  # end calc_kurtosisr
//' all.equal(HighFreq::calc_kurtosis(returns), 
//'   calc_kurtosisr(returns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(returns),
//'   Rcode=calc_kurtosisr(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile kurtosis
//' HighFreq::calc_kurtosis(returns, method = "quantile", alpha = 0.1)
//' # Compare HighFreq::calc_kurtosis() with quantile kurtosis
//' calc_kurtosisq <- function(x) {
//'   	quantiles <- quantile(x, c(0.25, 0.5, 0.75), type=5)
//'   	(quantiles[3] + quantiles[1] - 2*quantiles[2])/(quantiles[3] - quantiles[1])
//' }  # end calc_kurtosisq
//' all.equal(drop(HighFreq::calc_kurtosis(returns, method = "quantile")), 
//'   calc_kurtosisq(returns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(returns, method = "quantile"),
//'   Rcode=calc_kurtosisq(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric kurtosis
//' HighFreq::calc_kurtosis(returns, method = "nonparametric")
//' # Compare HighFreq::calc_kurtosis() with R nonparametric kurtosis
//' all.equal(drop(HighFreq::calc_kurtosis(returns, method = "nonparametric")), 
//'   (mean(returns)-median(returns))/sd(returns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(returns, method = "nonparametric"),
//'   Rcode=(mean(returns)-median(returns))/sd(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_kurtosis(arma::mat se_ries,
                        std::string method = "pearson", 
                        double alpha = 0.25) {
  // Return zeros if not enough data
  if (se_ries.n_rows < 3) {
    return arma::zeros<rowvec>(se_ries.n_cols);
  }  // end if
  
  // Switch statement for all the different methods of skew
  switch(calc_kurtosis_type(method)) {
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
  
}  // end calc_kurtosis


////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of kurtosis estimates over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{se_ries} A \emph{time series} or a \emph{matrix} of data.
//'    
//' @param \code{ste_p} The number of time periods between the end points.
//'
//' @param \code{look_back} The number of end points in the look-back interval.
//'
//' @param \code{method} A \emph{string} specifying the type of kurtosis.  (The
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
//'   The function \code{roll_kurtosis()} calculates a \emph{matrix} of kurtosis
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
//'   It passes the subset time series to the function \code{calc_kurtosis()}, which
//'   calculates the kurtosis.
//'   See the function \code{calc_kurtosis()} for a description of the kurtosis
//'   methods.
//'   
//'   For example, the rolling kurtosis at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{ste_p = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_kurtosis()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Define end points and start points
//' endp <- 1 + HighFreq::calc_endpoints(NROW(returns), ste_p)
//' startp <- HighFreq::calc_startpoints(endp, 3)
//' # Calculate the rolling kurtosis at 25 day end points, with a 75 day look-back
//' kurto_sis <- HighFreq::roll_kurtosis(returns, ste_p=25, look_back=3)
//' # Calculate the rolling kurtosis using R code
//' kurto_r <- sapply(1:NROW(endp), function(it) {
//'   HighFreq::calc_kurtosis(returns[startp[it]:endp[it], ])
//' })  # end sapply
//' # Compare the kurtosis estimates
//' all.equal(drop(kurto_sis), kurto_r, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_kurtosis(returns, ste_p=25, look_back=3),
//'   Rcode=sapply(1:NROW(endp), function(it) {
//'     HighFreq::calc_kurtosis(returns[startp[it]:endp[it], ])
//'   }),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_kurtosis(arma::mat se_ries, 
                        arma::uword ste_p = 1, 
                        arma::uword look_back = 11, 
                        std::string method = "pearson", 
                        double alpha = 0.25) {
  
  // Calculate end points
  arma::uword nrows = se_ries.n_rows;
  arma::uvec endp = calc_endpoints(nrows, ste_p);
  // Start points equal to end points lagged by look_back
  arma::uvec startp = calc_startpoints(endp, look_back);
  // Allocate kurtosis matrix
  arma::uword num_points = endp.n_elem;
  arma::mat kurto_sis = arma::zeros<mat>(num_points, se_ries.n_cols);
  
  // Perform loop over the endp
  for (arma::uword ep = 0; ep < num_points; ep++) {
    // Calculate kurtosis
    if (endp(ep) > startp(ep)) {
      kurto_sis.row(ep) = calc_kurtosis(se_ries.rows(startp(ep), endp(ep)), method);
    }  // end if
  }  // end for
  
  return kurto_sis;
  
}  // end roll_kurtosis


