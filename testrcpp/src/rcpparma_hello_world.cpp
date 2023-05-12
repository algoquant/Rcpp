// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);
	                     
    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}



////////////////////////////////////////////////////////////
//' Calculate the trailing weighted means of streaming \emph{time series} data
 //' using an online recursive formula.
 //' 
 //' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
 //' 
 //' @param \code{lambda} A decay factor which multiplies past estimates.
 //'   
 //' @param \code{weightv} A single-column \emph{matrix} of weights.
 //'
 //' @return A \emph{matrix} with the same dimensions as the input argument
 //'   \code{tseries}.
 //'
 //' @details
 //'   The function \code{run_mean()} calculates the trailing weighted means of
 //'   the streaming \emph{time series} data \eqn{p_t} by recursively weighting
 //'   present and past values using the decay factor \eqn{\lambda}. If the
 //'   \code{weightv} argument is equal to zero, then the function
 //'   \code{run_mean()} simply calculates the exponentially weighted moving
 //'   average value of the streaming \emph{time series} data \eqn{p_t}:
 //'   \deqn{
 //'     \bar{p}_t = \lambda \bar{p}_{t-1} + (1-\lambda) p_t = (1-\lambda) \sum_{j=0}^{n} \lambda^j p_{t-j}
 //'   }
 //'   
 //'   Some applications require applying additional weight factors, like for
 //'   example the volume-weighted average price indicator (VWAP).  Then the
 //'   streaming prices can be multiplied by the streaming trading volumes.
 //'   
 //'   If the argument \code{weightv} has the same number of rows as the argument
 //'   \code{tseries}, then the function \code{run_mean()} calculates the
 //'   trailing weighted means in two steps.
 //'   
 //'   First it calculates the trailing mean weights \eqn{\bar{w}_t}:
 //'   \deqn{
 //'     \bar{w}_t = \lambda \bar{w}_{t-1} + (1-\lambda) w_t
 //'   }
 //'   
 //'   Second it calculates the trailing mean products \eqn{\bar{w p}_t} of the
 //'   weights \eqn{w_t} and the data \eqn{p_t}:
 //'   \deqn{
 //'     \bar{w p}_t = \lambda \bar{w p}_{t-1} + (1-\lambda) w_t p_t
 //'   }
 //'   Where \eqn{p_t} is the streaming data, \eqn{w_t} are the streaming
 //'   weights, \eqn{\bar{w}_t} are the trailing mean weights, and \eqn{\bar{w p}_t}
 //'   are the trailing mean products of the data and the weights.
 //'   
 //'   The trailing mean weighted value \eqn{\bar{p}_t} is equal to the ratio of the
 //'   data and weights products, divided by the mean weights:
 //'   \deqn{
 //'     \bar{p}_t = \frac{\bar{w p}_t}{\bar{w}_t}
 //'   }
 //' 
 //'   The above online recursive formulas are convenient for processing live
 //'   streaming data because they don't require maintaining a buffer of past
 //'   data.
 //'   The formulas are equivalent to a convolution with exponentially decaying
 //'   weights, but they're much faster to calculate.
 //'   Using exponentially decaying weights is more natural than using a sliding
 //'   look-back interval, because it gradually "forgets" about the past data.
 //'   
 //'   The value of the decay factor \eqn{\lambda} must be in the range between
 //'   \code{0} and \code{1}.  
 //'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
 //'   values have a greater weight, and the trailing mean values have a stronger
 //'   dependence on past data.  This is equivalent to a long look-back
 //'   interval.
 //'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
 //'   past values have a smaller weight, and the trailing mean values have a
 //'   weaker dependence on past data.  This is equivalent to a short look-back
 //'   interval.
 //' 
 //'   The function \code{run_mean()} performs the same calculation as the
 //'   standard \code{R} function\cr\code{stats::filter(x=series, filter=lambda,
 //'   method="recursive")}, but it's several times faster.
 //' 
 //'   The function \code{run_mean()} returns a \emph{matrix} with the same
 //'   dimensions as the input argument \code{tseries}.
 //'   
 //' @examples
 //' \dontrun{
 //' # Calculate historical prices
 //' ohlc <- rutils::etfenv$VTI
 //' closep <- quantmod::Cl(ohlc)
 //' # Calculate the trailing means
 //' lambda <- 0.9
 //' meanv <- HighFreq::run_mean(closep, lambda=lambda, weightv=0)
 //' # Calculate trailing means using R code
 //' pricef <- (1-lambda)*filter(closep, 
 //'   filter=lambda, init=as.numeric(closep[1, 1])/(1-lambda), 
 //'   method="recursive")
 //' all.equal(drop(meanv), unclass(pricef), check.attributes=FALSE)
 //' 
 //' # Compare the speed of RcppArmadillo with R code
 //' library(microbenchmark)
 //' summary(microbenchmark(
 //'   Rcpp=HighFreq::run_mean(closep, lambda=lambda, weightv=0),
 //'   Rcode=filter(closep, filter=lambda, init=as.numeric(closep[1, 1])/(1-lambda), method="recursive"),
 //'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
 //'   
 //' # Calculate weights equal to the trading volumes
 //' weightv <- quantmod::Vo(ohlc)
 //' # Calculate the trailing weighted means
 //' meanw <- HighFreq::run_mean(closep, lambda=lambda, weightv=weightv)
 //' # Plot dygraph of the trailing weighted means
 //' datav <- xts(cbind(meanv, meanw), zoo::index(ohlc))
 //' colnames(datav) <- c("means trailing", "means weighted")
 //' dygraphs::dygraph(datav, main="Trailing Means") %>%
 //'   dyOptions(colors=c("blue", "red"), strokeWidth=2) %>%
 //'   dyLegend(show="always", width=500)
 //' }
 //' 
 //' @export
 // [[Rcpp::export]]
 arma::mat run_mean(const arma::mat& tseries, 
                    double lambda, // Decay factor which multiplies the past values 
                    const arma::colvec& weightv = 0) {
   
   arma::uword nrows = tseries.n_rows;
   arma::uword nweights = weightv.n_elem;
   arma::mat meanm(nrows, tseries.n_cols);
   double lambda1 = 1-lambda;
   
   if (nweights == 1) {
     meanm.row(0) = tseries.row(0);
     // Calculate means without weights
     for (arma::uword it = 1; it < nrows; it++) {
       // Calculate the means using the decay factor
       meanm.row(it) = lambda*meanm.row(it-1) + lambda1*tseries.row(it);
     }  // end for
   } else if (nweights == nrows) {
     // Calculate means with weights
     double meanw = weightv(0);
     double meanww;
     meanm.row(0) = meanw*tseries.row(0);
     for (arma::uword it = 1; it < nrows; it++) {
       // Calculate the means using the decay factor
       meanww = meanw;
       meanw = lambda*meanw + lambda1*weightv(it);
       meanm.row(it) = lambda*meanm.row(it-1) + lambda1*weightv(it)*tseries.row(it);
       // Divide by the mean weight
       meanm.row(it-1) = meanm.row(it-1)/meanww;
     }  // end for
     meanm.row(nrows-1) = meanm.row(nrows-1)/meanw;
   }  // end if
   
   return meanm;
   
 }  // end run_mean

