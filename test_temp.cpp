////////////////////////////
// Tests of C++ Functions
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_temp.cpp")

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

// using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;


////////////////////////////////////////////////////////////
//' Scrub the bad data in a \emph{time series} or a \emph{matrix}. 
//' 
//' @param \code{timeser} A single-column \emph{time series} or a
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
//'   series} \code{timeser} by shifting its elements by the number equal to the
//'   argument \code{lagg}.  For positive \code{lagg} values, the elements are
//'   shifted forward in time (down), and for negative \code{lagg} values they
//'   are shifted backward (up).
//'   
//'   The output \emph{vector} is padded with either zeros (the default), or
//'   with data from \code{timeser}, so that it has the same number of element
//'   as \code{timeser}.
//'   If the \code{lagg} is positive, then the first element is copied and added
//'   upfront.
//'   If the \code{lagg} is negative, then the last element is copied and added
//'   to the end.
//'   
//'   As a rule, if \code{timeser} contains returns data, then the output
//'   \emph{matrix} should be padded with zeros, to avoid data snooping.
//'   If \code{timeser} contains prices, then the output \emph{matrix} should
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
void do_scrub(arma::mat& timeser, double tolv=0.1) {
  
  bool isvalid = true;
  int nrows = timeser.n_rows;
  int nscrub = 0;
  double diffv = 0;
  
  // Allocate output matrix
  arma::mat posv = arma::zeros(nrows, 1);
  
  // Perform loop over the end points
  for (int it = 1; it < nrows; it++) {
    diffv = abs(timeser(it) - timeser(it-1));
    if ((diffv > tolv) && isvalid) {
      timeser(it) = timeser(it-1);
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


// Return TRUE (1) for the NA or Inf values in a matrix.
// [[Rcpp::export]]
arma::mat is_na(arma::mat& timeser) {
  
  // Allocate output matrix
  int nrows = timeser.n_rows;
  arma::mat posv = arma::zeros(nrows, 1);
  
  // cout << "timeser(0): " << timeser(0) << endl;
  // Perform loop over the end points
  for (int it = 0; it < nrows; it++) {
    if (std::isnan(timeser(it)) || std::isinf(timeser(it))) {
    // if (timeser(it).has_nan() || timeser(it).has_inf()) {
      posv(it) = 1;
    }  // end if
  }  // end for
  
  return posv;
  
}  // end is_na


// Return TRUE (1) if the matrix has NA or Inf values.
// [[Rcpp::export]]
bool has_na(arma::mat& timeser) {
  
  if (timeser.has_nan() || timeser.has_inf()) {
    return true;
  } else {
    return false;
  }  // end if
    
}  // end has_na



// Calculate the exponential moving average (EMA) of streaming \emph{time
// series} data using an online recursive formula.
// Handle the NA or Inf values.
// [[Rcpp::export]]
arma::mat run_mean(const arma::mat& timeser, 
                   double lambdaf, // Decay factor which multiplies the past values 
                   const arma::colvec& weightv = 0) {
  
  // Allocate output matrix
  int nrows = timeser.n_rows;
  arma::mat posv = arma::zeros(nrows, 1);
  
  // cout << "timeser(0): " << timeser(0) << endl;
  // Perform loop over the end points
  for (int it = 0; it < nrows; it++) {
    if (std::isnan(timeser(it)) || std::isinf(timeser(it))) {
      // if (timeser(it).has_nan() || timeser(it).has_inf()) {
      posv(it) = 1;
    }  // end if
  }  // end for
  
  return posv;
  
}  // end run_mean


// Calculate the weights of a portfolio by selecting the nth smallest drawdowns.
// Select the values of the nth smallest elements in a vector.
// [[Rcpp::export]]
arma::mat calc_weights(const arma::mat& drawv, // Vector of drawdowns
                       int nums = 3, // Number of stocks to select
                       double threshd = -0.1) { // Threshold for the drawdown
  
  arma::uword nrows = drawv.n_elem;
  arma::mat weightv = drawv;
  
  // Loop over the drawdowns if the number of stocks to select is less than the number of drawdowns
  if (nums < nrows) {

    // Calculate the ranks of the drawdowns
    arma::uvec rankk = (arma::sort_index(arma::sort_index(drawv)));
    
    // Loop over the drawdowns
    for (uword it = 0; it < nrows; it++) {
      
      if ((rankk[it] >= nums) || (drawv[it] > threshd))
        weightv[it] = 0;
      
    }  // end for
    
  }  // end if

  // Scale the weights to reduce leverage
  // double sumw = abs(arma::sum(arma::sum(weightv)));
  // if (sumw > 0) {
  //   weightv = weightv/sumw;
  // }  // end if
  
  return weightv;
  
}  // end calc_weights



// Calculate the weights of a portfolio over time.
// [[Rcpp::export]]
arma::mat run_weights(const arma::mat& retm, // Time series of returns
                      double lambdaf = 0.3, // Decay factor which multiplies the past returns 
                      int nums = 3, // Number of stocks to select
                      double threshd = -0.1, // Threshold for the drawdown
                      double threshvar = 0.0004) { // Threshold for the variance

  arma::uword nrows = retm.n_rows;
  arma::uword ncols = retm.n_cols;
  arma::mat meanm = arma::zeros(nrows, ncols); // Mean matrix
  arma::mat varm = arma::zeros(nrows, ncols); // Variance matrix
  arma::mat drawm = arma::zeros(nrows, ncols); // Drawdown matrix
  arma::mat weightm = arma::zeros(nrows, ncols);
  double drawd; // Drawdown value
  double lambda1 = 1-lambdaf;
  double lambda2 = pow(lambdaf, 2);
  double lambda21 = 1-lambda2;
  
  // Copy the first row of the returns to the mean matrix
  // meanm.row(0) = retm.row(0);
  double valuc = retm(0, 0);
  for (arma::uword cn = 0; cn < ncols; cn++) {
    valuc = retm(0, cn);
    if (std::isnan(valuc) || std::isinf(valuc)) {
      // Set the mean equal to zero
      meanm(0, cn) = 0;
      varm(0, cn) = 1;
      drawm(0, cn) = 0;
    } else {
      // Set the mean equal to the current return
      meanm(0, cn) = valuc;
      varm(0, cn) = max(pow(valuc, 2), 0.00001);
      // Calculate the weight as the Kelly ratio
      // drawm(0, cn) = meanm(0, cn)/varm(0, cn);
      drawm(0, cn) = meanm(0, cn);
    }  // end if
  }  // end column for
  
  // Loop over the rows of the time series
  for (arma::uword rn = 1; rn < nrows; rn++) {
    
    // Loop over the columns and check for any NA or Inf values
    for (arma::uword cn = 0; cn < ncols; cn++) {
      
      // Calculate the means and variances
      valuc = retm(rn, cn);
      if (std::isnan(valuc) || std::isinf(valuc)) {
        // Set the mean equal to zero
        meanm(rn, cn) = 0;
        varm(rn, cn) = 1;
        drawm(rn, cn) = 0;
      } else {
        // Calculate the mean using the decay factor
        meanm(rn, cn) = lambdaf*meanm(rn-1, cn) + lambda1*valuc;
        varm(rn, cn) = lambda2*varm(rn-1, cn) + lambda21*max(pow(valuc - meanm(rn, cn), 2), 0.00001);
      }  // end if
      
      // Calculate the drawdown if the variance is below the threshold
      if (varm(rn, cn) < threshvar) {
        // Calculate the drawdown equal to the Kelly ratio
        // drawd = meanm(rn, cn)/varm(rn, cn);
        // Calculate the drawdown equal to the mean
        drawd = meanm(rn, cn);
        // Reset the drawdown if it becomes more negative
        if (drawd < drawm(rn-1, cn)) {
          drawm(rn, cn) = drawd;
        } else if (drawd > 0) {
          // Reset the drawdown to zero if it's positive
          drawm(rn, cn) = 0;
        } else {
          // Carry forward the previous drawdown
          drawm(rn, cn) = drawm(rn-1, cn);
        }  // end if
      } else {
        // Set the drawdown zero if the variance is above the threshold
        drawm(rn, cn) = 0;
      }  // end if
      
    }  // end column for
    
    // Calculate the weights by selecting the nth smallest drawdowns and rescaling them
    weightm.row(rn) = calc_weights(drawm.row(rn), nums, threshd);
    // weightm.row(rn) = drawm.row(rn);
      
  }  // end row for
  
  return -weightm;
  
}  // end run_weights


// Calculate the downside variance of a \emph{time series} using an online recursive formula.
// [[Rcpp::export]]
arma::mat run_vard(const arma::mat& timeser, 
                   double lambdaf, // Decay factor which multiplies the past values
                   double retarg = 0.0) { // Target return for the downside variance
  
  arma::uword nrows = timeser.n_rows;
  arma::uword ncols = timeser.n_cols;
  // arma::mat meanm = arma::zeros(nrows, ncols);
  arma::mat vars = arma::zeros(nrows, ncols);
  // double lambda1 = 1-lambdaf;
  double lambda2 = pow(lambdaf, 2);
  double lambda21 = 1-lambda2;
  
  // Perform loop over the rows
  // meanm.row(0) = timeser.row(0);
  arma::mat varow = timeser.row(0) - retarg;
  varow.transform([](double x) {return min(x, 0.0);});
  vars.row(0) = arma::square(varow);
  // vars.row(0) = 0;
  if (!(timeser.has_nan() || timeser.has_inf())) {
    // No NA or Inf values
    for (arma::uword it = 1; it < nrows; it++) {
      // Calculate the means using the decay factor
      // meanm.row(it) = lambdaf*meanm.row(it-1) + lambda1*timeser.row(it);
      // Variance is the weighted sum of the past variance and the square of the data minus its mean
      varow = timeser.row(it) - retarg;
      varow.transform([](double x) {return min(x, 0.0);});
      vars.row(it) = lambda2*vars.row(it-1) + lambda21*arma::square(varow);
    }  // end for
  } else {
    // Loop over the columns because of the NA or Inf values.
    // Calculate variance with NA or Inf values
    double valuc = timeser(0, 0);
    // meanm.row(0) = timeser.row(0);
    for (arma::uword rn = 1; rn < nrows; rn++) {
      for (arma::uword cn = 0; cn < ncols; cn++) {
        // Calculate the negative (downside) returns
        valuc = min(timeser(rn, cn) - retarg, 0.0);
        valuc = pow(valuc, 2);
        if (std::isnan(vars(rn-1, cn)) || std::isinf(vars(rn-1, cn)) || std::isnan(valuc) || std::isinf(valuc)) {
          // Set the mean and variance equal to the current value
          // meanm(rn, cn) = valuc;
          vars(rn, cn) = valuc;
        } else {
          // Calculate the mean using the decay factor
          // meanm(rn, cn) = lambdaf*meanm(rn-1, cn) + lambda1*valuc;
          // Variance is the weighted sum of the past variance and the square of the data minus its mean
          vars(rn, cn) = lambda2*vars(rn-1, cn) + lambda21*valuc;
        }  // end if
      }  // end column for
    }  // end row for
    // meanm = timeser;
  }  // end if
  
  // return arma::join_rows(meanm, vars);
  return vars;
  
}  // end run_vard



// [[Rcpp::export]]
arma::mat get_pos(arma::mat& drawv) {

  return drawv.transform([](double x) {return max(x, 0.0);});

}  // end get_pos

