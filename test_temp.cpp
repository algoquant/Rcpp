////////////////////////////
// Functions to test C++ syntax with Armadillo
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_temp.cpp")

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


////////////////////////////////////////////////
// Tests misc




////////////////////////////////////////////////
// Included to facilitate Tests - remove later - don't migrate

////////////////////////////////////////////////////////////
// Define C++ enum type for the different methods for regularization,
// calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum meth_od {moment, least_squares, quantile, nonparametric, regular, rank_sharpe, 
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
  else if (method == "regular")
    return meth_od::regular;
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



arma::uvec calc_endpoints(arma::uword length, arma::uword step = 1, arma::uword stub = 0) {
  
  arma::uword remainder = length % step;
  arma::uvec endp;
  
  if ((stub == 0) & (remainder == 0)) {
    // No stub interval
    endp = arma::regspace<uvec>(step, step, length);
  } else if ((stub == 0) & (remainder > 0)) {
    // Add stub interval at end
    endp = arma::regspace<uvec>(step, step, length + step);
    endp.back() = length;
  } else if ((stub > 0) & (remainder == 0)) {
    // Add initial stub interval equal to stub
    endp = arma::regspace<uvec>(stub, step, length + step);
    endp.back() = length;
  } else if ((stub > 0) & (remainder > 0) & (stub == remainder)) {
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



arma::uvec calc_startpoints(arma::uvec endp, arma::uword look_back) {
  
  arma::uword num_pts = endp.n_elem;
  arma::uvec startp = arma::join_cols(arma::zeros<uvec>(look_back), 
                                      endp.subvec(0, num_pts - look_back - 1) + 1);
  
  return startp;
  
}  // end calc_startpoints



//' @export
// [[Rcpp::export]]
arma::mat diff_it(const arma::mat& tseries, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
  arma::uword num_rows = (tseries.n_rows-1);
  arma::uword num_cols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    // Matrix difference without padding
    arma::mat diff_mat = (tseries.rows(lagg, num_rows) - tseries.rows(0, num_rows - lagg));
    if (pad_zeros) {
      // Pad diff_mat with zeros at the front
      return arma::join_cols(arma::zeros<mat>(lagg, num_cols), diff_mat);
    } else {
      // Don't pad the output
      return diff_mat;
    }  // end if pad_zeros
  } else {
    // Negative lag
    // Matrix difference without padding
    arma::mat diff_mat = (tseries.rows(0, num_rows + lagg) - tseries.rows(-lagg, num_rows));
    if (pad_zeros) {
      // Pad diff_mat with zeros at the back
      return arma::join_cols(diff_mat, arma::zeros<mat>(-lagg, num_cols));
    } else {
      // Don't pad the output
      return diff_mat;
    }  // end if pad_zeros
  }  // end if lagg
  
}  // end diff_it





////////////////////////////////////////////////
// Test versions to be migrated to package HighFreq




////////////////////////////////////////////////
// Old stuff - can be deleted later

// Calculate the rolling maximum or minimum of streaming data using a lambda decay factor
//' @export
// [[Rcpp::export]]
arma::mat run_maxmin(arma::mat tseries, double lambda, bool calc_max = true) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat maxmin = tseries;
  arma::mat means = tseries;
  double lambda1 = 1-lambda;
  
  if (calc_max) {
    // Perform loop over rows
    for (arma::uword it = 1; it < num_rows; it++) {
      // Calculate the mean as a weighted sum
      means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
      // Calculate the max from a weighted sum
      maxmin.row(it) = arma::max(maxmin.row(it), means.row(it-1) + lambda*(maxmin.row(it-1) - means.row(it-1)));
    }  // end for
  } else {
    // Perform loop over rows
    for (arma::uword it = 1; it < num_rows; it++) {
      // Calculate the mean as a weighted sum
      means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
      // Calculate the max from a weighted sum
      maxmin.row(it) = arma::min(maxmin.row(it), means.row(it-1) + lambda*(maxmin.row(it-1) - means.row(it-1)));
    }  // end for
  }  // end if
  
  return maxmin;
  
}  // end run_maxmin


// Calculate the rolling maximum of streaming data using a lambda decay factor
//' @export
// [[Rcpp::export]]
arma::colvec armax(arma::colvec tseries, arma::colvec tseries2) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::colvec maxs = tseries;
  
  // Perform loop over rows
  for (arma::uword it = 0; it < num_rows; it++) {
    // Calculate the mean as a weighted sum
    // means.row(it) = lambda1*means.row(it) + lambda*means.row(it-1);
    // Calculate the max from a weighted sum
    // maxs.row(it) = arma::max(maxs.row(it), means.row(it-1) + lambda*(maxs.row(it-1) - means.row(it-1)));
    maxs.row(it) = arma::max(tseries.row(it), tseries2.row(it));
  }  // end for
  
  return maxs;
  
}  // end armax



// Demonstration of using default NULL arguments in Rcpp code.
// calc_endpoints_null() implements default NULL arguments in Rcpp code.
// Calculate the end points with a stub interval.
bool calc_endpoints_null(arma::mat& se_ries,
                         arma::uword ste_p = 1,
                         Rcpp::Nullable<int> stu_b = R_NilValue, 
                         Rcpp::Nullable<Rcpp::IntegerVector> end_points = R_NilValue) {
  
  arma::uword num_rows = se_ries.n_rows;
  arma::uvec end_p;
  // arma::uvec end_p = arma::ones<uvec>(1);
  // end_p.reset();
  // bool is_empty;
  
  if (end_points.isNotNull()) {
    // Simply copy end_points
    end_p = Rcpp::as<uvec>(end_points);
  } else if (stu_b.isNotNull()) {
    // Calculate end points with stu_b
    end_p = arma::regspace<uvec>(Rcpp::as<uword>(stu_b), ste_p, num_rows + ste_p);
    end_p = end_p.elem(find(end_p < num_rows));
  }  // end if
  
  // Return Boolean
  if (end_p.is_empty())
    // end_p is empty if arguments end_points and stu_b are both NULL
    return true;
  else  
    return false;
  
}  // end calc_endpoints_null



////////////////////////////////////////////////
// Experimental functions.


////////////////////////////////////////////////////////////
//' Calculate by reference the rolling convolutions (weighted sums) of a
//' \emph{time series} with a \emph{column vector} of weights.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//' 
//' @param \code{weights} A \emph{column vector} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{tseries}.
//'
//' @details 
//'   The function \code{roll_conv_ref()} is experimental and doesn't quite work yet.
//'   
//'   The function \code{roll_conv_ref()} calculates the convolutions of the
//'   \emph{matrix} columns with a \emph{column vector} of weights.  It performs a loop
//'   down over the \emph{matrix} rows and multiplies the past (higher) values
//'   by the weights.  It calculates the rolling weighted sums of the past
//'   values.
//'   
//'   The function \code{roll_conv_ref()} accepts a \emph{pointer} to the argument
//'   \code{tseries}, and replaces the old \emph{matrix} values with the
//'   weighted sums.
//'   It performs the calculation in place, without copying the \emph{matrix} in
//'   memory (which greatly increases the computation speed).
//'   
//'   The function \code{roll_conv_ref()} uses the \code{RcppArmadillo} function
//'   \code{arma::conv2()}. It performs a similar calculation to the standard
//'   \code{R} function \code{filter(x=re_turns, filter=weight_s,
//'   method="convolution", sides=1)}, but it's over \code{6} times faster, and
//'   it doesn't produce any leading \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Calculate a time series of prices
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("IEF", "VTI")])
//' # Create simple weights equal to a 1 value plus zeros
//' weight_s <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sums
//' weight_ed <- HighFreq::roll_conv_ref(re_turns, weight_s)
//' # Compare with original
//' all.equal(coredata(re_turns), weight_ed, check.attributes=FALSE)
//' # Second example
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*(1:11))
//' weight_s <- matrix(weight_s/sum(weight_s), nc=1)
//' # Calculate rolling weighted sums
//' weight_ed <- HighFreq::roll_conv_ref(re_turns, weight_s)
//' # Calculate rolling weighted sums using filter()
//' filter_ed <- filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filter_ed[-(1:11), ], weight_ed[-(1:11), ], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_conv_ref(arma::mat& tseries, const arma::mat& weights) {
  
  arma::uword look_back = weights.n_rows-1;
  arma::uword num_rows = tseries.n_rows-1;
  
  // Calculate the convolutions
  arma::mat convmat = arma::conv2(tseries, weights, "full");
  // Copy the convolutions
  tseries.rows(look_back, num_rows) = convmat.rows(look_back, num_rows);
  // tseries = arma::conv2(tseries, weights, "same");
  
  return weights;
  
}  // end roll_conv_ref





// double armax(double lambda1, double lambda2) {
//   
//   return std::max(lambda1, lambda2);
//   
// }  // end armax

