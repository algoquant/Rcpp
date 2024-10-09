////////////////////////////
// Tests of C++ Functions
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_temp.cpp")

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;


////////////////////////////////////////////////////////////
//' Scrub the bad data in a \emph{time series} or a \emph{matrix}. 
//' 
//' @param \code{tseries} A single-column \emph{time series} or a
//'   \emph{vector}.
//'
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag.
//'   (The default is \code{lagg = 1}.)
//'
//' @param \code{pad_zeros} \emph{Boolean} argument: Should the output be padded
//'   with zeros? (The default is \code{pad_zeros = TRUE}.)
//'
//' @return A column \emph{vector} with the same number of elements as the input
//'   time series.
//'
//' @details
//'   The function \code{lag_vec()} applies a lag to the input \emph{time
//'   series} \code{tseries} by shifting its elements by the number equal to the
//'   argument \code{lagg}.  For positive \code{lagg} values, the elements are
//'   shifted forward in time (down), and for negative \code{lagg} values they
//'   are shifted backward (up).
//'   
//'   The output \emph{vector} is padded with either zeros (the default), or
//'   with data from \code{tseries}, so that it has the same number of element
//'   as \code{tseries}.
//'   If the \code{lagg} is positive, then the first element is copied and added
//'   upfront.
//'   If the \code{lagg} is negative, then the last element is copied and added
//'   to the end.
//'   
//'   As a rule, if \code{tseries} contains returns data, then the output
//'   \emph{matrix} should be padded with zeros, to avoid data snooping.
//'   If \code{tseries} contains prices, then the output \emph{matrix} should
//'   be padded with the prices.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' retp <- rnorm(1e6)
//' # Compare lag_vec() with rutils::lagit()
//' all.equal(drop(HighFreq::lag_vec(retp)), 
//'   rutils::lagit(retp))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::lag_vec(retp),
//'   Rcode=rutils::lagit(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
void do_scrub(arma::mat& pricev, double tolv=0.1) {
  
  bool isvalid = true;
  int nrows = pricev.n_rows;
  int nscrub = 0;
  double diffv = 0;
  
  // Allocate output matrix
  arma::mat posv = arma::zeros(nrows, 1);
  
  // Perform loop over the end points
  for (int it = 1; it < nrows; it++) {
    diffv = abs(pricev(it) - pricev(it-1));
    if ((diffv > tolv) && isvalid) {
      pricev(it) = pricev(it-1);
      isvalid = false;
      nscrub++;
    } else {
      isvalid = true;
    }  // end if
  }  // end for
  
  cout << "nscrub: " << nscrub << endl;
  // return nscrub;
  
}  // end do_scrub




//' Find the positions of the closest values in a single-column \emph{matrix}
//' that match those in the original \emph{matrix}.
//' Given two vectors of dates, ts1 and ts2, find the closest 
//' dates in ts2 that match those in ts1.
//' The lengths of ts1 and ts2 don't have to be equal.
//' In practice it's better to select ts1 to be the shorter series.
//' @export
// [[Rcpp::export]]
arma::mat calc_nearest(const arma::mat& ts1, const arma::mat& ts2) {

  int nrows1 = ts1.n_rows;
  int nrows2 = ts2.n_rows;
  int npts = 5;
  int startp = 0;
  int endp = nrows2 - 1;
  int posi = 0;
  double diffv = 0;
  double diffm = 0;
  
  // Allocate output matrix
  arma::mat posv = arma::zeros(nrows1, 1);

  // Perform loop over the end points
  for (int it = 0; it < nrows1; it++) {
    startp = max(posi - npts, 0);
    // endp = min(posi + npts, nrows2);
    endp = nrows2;
    diffm = 10;
    for (int jt = startp; jt < endp; jt++) {
      diffv = abs(ts1(it) - ts2(jt));
      if (diffv < diffm) {
        diffm = diffv;
        posv(it) = jt + 1;
        posi = jt;
      }  // end if
    }  // end for
  }  // end for
  
  return posv;
  
}  // end calc_nearest


// [[Rcpp::export]]
arma::mat is_na(arma::mat& pricev) {
  
  // Allocate output matrix
  int nrows = pricev.n_rows;
  arma::mat posv = arma::zeros(nrows, 1);
  
  // cout << "pricev(0): " << pricev(0) << endl;
  // Perform loop over the end points
  for (int it = 0; it < nrows; it++) {
    if (std::isnan(pricev(it)) || std::isinf(pricev(it))) {
    // if (pricev(it).has_nan() || pricev(it).has_inf()) {
      posv(it) = 1;
    }  // end if
  }  // end for
  
  return posv;
  
}  // end is_na


// [[Rcpp::export]]
bool has_na(arma::mat& pricev) {
  
  if (pricev.has_nan() || pricev.has_inf()) {
    return true;
  } else {
    return false;
  }  // end if
    
}  // end has_na




