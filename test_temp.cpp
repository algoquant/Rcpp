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


//' @export
// [[Rcpp::export]]
arma::vec lag_vec(arma::vec& vec_tor, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
  arma::uword num_rows = (vec_tor.n_elem-1);
  
  if (lagg > 0) {
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros(lagg), 
                             vec_tor.subvec(0, num_rows-lagg));
    } else {
      // Pad front with first element of vec_tor
      return arma::join_cols(arma::repelem(vec_tor.subvec(0, 0), lagg, 1), 
                             vec_tor.subvec(0, num_rows-lagg));
    }  // end if
  } else {
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(vec_tor.subvec(-lagg, num_rows), 
                             arma::zeros(-lagg));
    } else {
      // Pad back with last element of vec_tor
      return arma::join_cols(vec_tor.subvec(-lagg, num_rows), 
                             arma::repelem(vec_tor.subvec(num_rows, num_rows), -lagg, 1));
    }  // end if
  }  // end if
  
}  // end lag_vec



//' @export
// [[Rcpp::export]]
arma::mat lag_it(arma::mat& t_series, 
                 arma::sword lagg = 1, 
                 bool pad_zeros = true) {
  
  arma::uword num_rows = (t_series.n_rows-1);
  arma::uword num_cols = t_series.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros(lagg, num_cols), 
                             t_series.rows(0, num_rows-lagg));
    } else {
      // Pad front with first element of t_series
      return arma::join_cols(arma::repmat(t_series.rows(0, 0), lagg, 1), 
                             t_series.rows(0, num_rows-lagg));
    }  // end if
  } else {
    // Negative lag
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(t_series.rows(-lagg, num_rows), 
                             arma::zeros(-lagg, num_cols));
    } else {
      // Pad back with last element of t_series
      return arma::join_cols(t_series.rows(-lagg, num_rows), 
                             arma::repmat(t_series.rows(num_rows, num_rows), -lagg, 1));
    }  // end if
  }  // end if
  
  // Old code below
  // if (lagg > 0)
  //   // Positive lag
  //   return arma::join_cols(arma::repelem(t_series.row(0), lagg, 1), 
  //                          t_series.rows(0, num_rows-lagg));
  // else
  //   // Negative lag
  //   return arma::join_cols(t_series.rows(-lagg, num_rows), 
  //                          arma::repelem(t_series.row(num_rows), -lagg, 1));
  
}  // end lag_it



//' @export
// [[Rcpp::export]]
arma::mat diff_it(arma::mat& t_series, 
                  arma::uword lagg = 1, 
                  bool padd = true) {
  
  arma::uword num_rows = (t_series.n_rows-1);
  // Matrix difference without padding
  arma::mat diff_mat = (t_series.rows(lagg, num_rows) - t_series.rows(0, num_rows - lagg));
  
  if (padd)
    // Pad diff_mat with warmup period at the beginning
    return arma::join_cols(t_series.rows(0, lagg - 1), diff_mat);
  else
    // Don't pad the output
    return diff_mat;
  
}  // end diff_it



//' @export
// [[Rcpp::export]]
arma::vec calc_var_vec_agg(arma::vec& re_turns, 
                           arma::uword lagg = 1) {
  // Cumulative re_turns
  // arma::vec cum_sum = arma::cumsum(re_turns);
  
  return (re_turns - lag_vec(re_turns, lagg));
  // return (cum_sum - lag_vec(cum_sum, lagg));
  // return arma::var(cum_sum);
}  // end calc_var_vec_agg



//' @export
// [[Rcpp::export]]
arma::uvec calc_endpoints(arma::uword len_gth, arma::uword ste_p, bool front=true) {
  
  // Calculate number of intervals that fit over len_gth
  arma::uword num_points = len_gth/ste_p;
  arma::uvec end_p;
  
  if (len_gth == ste_p*num_points) {
    // No stub interval
    // Include the first placeholder end point - legacy code
    // end_p = arma::cumsum(arma::ones<uvec>(num_points));
    // end_p = arma::regspace<uvec>(0, ste_p, len_gth);
    end_p = arma::regspace<uvec>(ste_p, ste_p, len_gth);
  } else {
    // Need to add stub interval
    // Include the first placeholder end point - legacy code
    // end_p = arma::cumsum(arma::ones<uvec>(num_points + 1));
    // end_p = arma::regspace<uvec>(0, ste_p, len_gth + ste_p);
    end_p = arma::regspace<uvec>(ste_p, ste_p, len_gth + ste_p);
    if (front) {
      // Stub interval at beginning
      end_p = end_p - ste_p + len_gth % ste_p;
    } else {
      // Stub interval at end
      // The last end point must be equal to len_gth
      // end_p(num_points + 1) = len_gth;
      end_p(num_points) = len_gth;
    }  // end if
  }  // end if
  
  // Set the first end point to zero - it's a placeholder
  // end_p(0) = 0;
  // Subtract 1 from end_p because indexing starts at 0
  end_p = end_p - 1;
  return end_p;
  
}  // end calc_endpoints


//' @export
// [[Rcpp::export]]
arma::uvec calc_startpoints(arma::uvec end_points, arma::uword look_back) {

  arma::uword num_points = end_points.n_elem;
  arma::uvec start_p = arma::join_cols(arma::zeros<uvec>(look_back), 
                                       end_points.subvec(0, num_points - look_back - 1) + 1);
  
  return start_p;
  
}  // end calc_startpoints


//' @export
// [[Rcpp::export]]
arma::uvec calc_endpoints_stub(arma::uword len_gth, arma::uword ste_p, arma::uword stu_b = 0) {
  
  // Calculate number of intervals that fit over len_gth
  // arma::uword num_points = len_gth/ste_p;
  arma::uvec endp = arma::regspace<uvec>(stu_b, ste_p, len_gth + ste_p);
  
  return endp.elem(find(endp <= len_gth));
  
}  // end calc_endpoints_stub




////////////////////////////////////////////////////////////
//' Calculate a vector of end points that divides a vector into equal intervals.
//'
//' @param \code{len_gth} An \emph{integer} equal to the length of the vector to
//'   be divide into equal intervals.
//'   
//' @param \code{ste_p} The number of elements in each interval.
//' 
//' @param \code{front} \emph{Boolean} argument: if \code{TRUE} then add a stub
//'   interval at the beginning, else add a stub interval at the end.  (default
//'   is \code{TRUE})
//'
//' @return An \emph{integer} vector of equally spaced end points (vector of
//'   integers).
//'
//' @details The end points are a vector of integers which divide the vector of
//'   length equal to \code{len_gth} into equally spaced intervals.
//'   If a whole number of intervals doesn't fit over the vector, then
//'   \code{calc_endpoints()} adds a stub interval either at the beginning (the
//'   default) or at the end.
//'   The end points are shifted by \code{-1} because indexing starts at
//'   \code{0} in \code{C++} code.
//'
//'   The function \code{calc_endpoints()} is similar to the function
//'   \code{rutils::calc_endpoints()} from package
//'   \href{https://github.com/algoquant/rutils}{rutils}.
//'   
//'   The end points produced by \code{calc_endpoints()} don't include the first
//'   placeholder end point, which is usually equal to zero.
//'   For example, consider the end points for a vector of length \code{20}
//'   divided into intervals of length \code{5}: \code{0, 5, 10, 15, 20}.
//'   In order for all the differences between neighboring end points to be
//'   equal to \code{5}, the first end point must be equal to \code{0}.
//'   The first end point is a placeholder and doesn't correspond to any vector
//'   element.
//'   
//'   This works in \code{R} code because the vector element corresponding to
//'   index \code{0} is empty.  For example, the \code{R} code: \code{(4:1)[c(0,
//'   1)]} produces \code{4}.  So in \code{R} we can select vector elements
//'   using the end points starting at zero.
//'   
//'   In \code{C++} the end points must be shifted by \code{-1} because indexing
//'   starts at \code{0}: \code{-1, 4, 9, 14, 19}.  But there is no vector
//'   element corresponding to index \code{-1}. So in \code{C++} we cannot
//'   select vector elements using the end points starting at \code{-1}. The
//'   solution is to drop the first placeholder end point.
//'   
//' @examples
//' # Calculate end points without a stub interval
//' HighFreq::calc_endpoints(25, 5)
//' # Calculate end points with initial stub interval
//' HighFreq::calc_endpoints(23, 5)
//' # Calculate end points with a stub interval at the end
//' HighFreq::calc_endpoints(23, 5, FALSE)
//'
//' @export
// [[Rcpp::export]]
arma::rowvec calc_var(arma::mat& re_turns, 
                      arma::uword ste_p = 1) {
  
  // arma::uword num_cols = re_turns.n_cols;
  
  if (ste_p == 1)
    // Calculate the variance without aggregations
    return arma::var(re_turns);
  else {
    arma::mat cum_sum = arma::cumsum(re_turns, 0);
    // Calculate the variance of aggregated returns
    // arma::uvec end_p = calc_endpoints(re_turns.n_rows, ste_p);
    // cum_sum = cum_sum.rows(end_p);
    // return arma::var(diff_it(cum_sum, 1, false));
    
    // Perform loop over the stubs
    arma::uword num_rows = re_turns.n_rows;
    arma::mat aggs;
    arma::uvec end_p;
    arma::mat var_s(ste_p, re_turns.n_cols);
    for (arma::uword stu_b = 0; stu_b < ste_p; stu_b++) {
      end_p = arma::regspace<uvec>(stu_b, ste_p, num_rows + ste_p);
      end_p = end_p.elem(find(end_p < num_rows));
      aggs = cum_sum.rows(end_p);
      var_s.row(stu_b) = arma::var(diff_it(aggs, 1, false));
    }  // end for
    return mean(var_s);
  }  // end if
  
}  // end calc_var




////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sum over a \emph{time series} or a
//' \emph{matrix} using \emph{Rcpp}.
//' 
//' @param \code{t_series} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of data points included in calculating the rolling sum (the default
//'   is \code{look_back = 1}).
//'   
//' @param \code{stu_b} An \emph{integer} value equal to the first stub interval
//'   for calculating the end points.
//' 
//' @param \code{end_points} An \emph{unsigned integer} vector of end
//' points.
//'   
//' @param \code{weight_s} A column \emph{vector} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{t_series}.
//'
//' @details The function \code{roll_sum()} calculates the rolling sums over the
//'   columns of the \code{t_series} data.  
//'   The sums are calculated over a number of data points equal to
//'   \code{look_back}.
//'   
//'   The function \code{roll_sum()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{t_series}.
//' 
//'   The arguments \code{stu_b}, \code{end_points}, and \code{weight_s} are
//'   optional.
//'   
//'   If either the arguments \code{stu_b} or \code{end_points} are supplied,
//'   then the rolling sums are calculated at the end points. 
//'   
//'   If only the argument \code{stu_b} is supplied, then the end points are
//'   calculated from the \code{stu_b} and \code{look_back} arguments. The first
//'   end point is equal to \code{stu_b} and the end points are spaced
//'   \code{look_back} periods apart.
//'   
//'   If the argument \code{weight_s} is supplied, then weighted sums are
//'   calculated.
//'   Then the function \code{roll_sum()} calculates the rolling weighted sums
//'   of the past values.
//'   
//'   The function \code{roll_sum()} calculates the rolling weighted sums as
//'   convolutions of the \code{t_series} columns with the \emph{vector} of
//'   weights using the \code{RcppArmadillo} function \code{arma::conv2()}.
//'   It performs a similar calculation to the standard \code{R} function
//'   \code{stats::filter(x=t_series, filter=weight_s, method="convolution",
//'   sides=1)}, but it's over \code{6} times faster, and it doesn't produce any
//'   leading \code{NA} values. using fast \emph{RcppArmadillo} \code{C++} code.
//'   The function \code{roll_sum()} is several times faster than
//'   \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Create series of historical returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "IEF")])
//' # Define parameters
//' look_back <- 22
//' stu_b <- 21
//' # Calculate rolling sums at each point
//' c_sum <- HighFreq::roll_sum(re_turns, look_back=look_back)
//' r_sum <- rutils::roll_sum(re_turns, look_back=look_back)
//' all.equal(c_sum, coredata(r_sum), check.attributes=FALSE)
//' r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
//' lag_sum <- rbind(matrix(numeric(2*look_back), nc=2), r_sum[1:(NROW(r_sum) - look_back), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points
//' c_sum <- HighFreq::roll_sum(re_turns, look_back=look_back, stu_b=stu_b)
//' end_p <- (stu_b + look_back*(0:(NROW(re_turns) %/% look_back)))
//' end_p <- end_p[end_p < NROW(re_turns)]
//' r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
//' r_sum <- r_sum[end_p+1, ]
//' lag_sum <- rbind(numeric(2), r_sum[1:(NROW(r_sum) - 1), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points - pass in end_points
//' c_sum <- HighFreq::roll_sum(re_turns, end_points=end_p)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*(1:11))
//' weight_s <- matrix(weight_s/sum(weight_s), nc=1)
//' # Calculate rolling weighted sum
//' c_sum <- HighFreq::roll_sum(re_turns, weight_s=weight_s)
//' # Calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
//' all.equal(c_sum[-(1:11), ], filter_ed[-(1:11), ], check.attributes=FALSE)
//' 
//' # Calculate rolling weighted sums at end points
//' c_sum <- HighFreq::roll_sum(re_turns, end_points=end_p, weight_s=weight_s)
//' all.equal(c_sum, filter_ed[end_p+1, ], check.attributes=FALSE)
//' 
//' # Create simple weights equal to a 1 value plus zeros
//' weight_s <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_sum(re_turns, weight_s)
//' # Compare with original
//' all.equal(coredata(re_turns), weight_ed, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sum(arma::mat& t_series,
                   arma::uword look_back = 1,
                   Rcpp::Nullable<int> stu_b = R_NilValue, 
                   Rcpp::Nullable<Rcpp::IntegerVector> end_points = R_NilValue, 
                   Rcpp::Nullable<Rcpp::NumericVector> weight_s = R_NilValue) {
  
  arma::uword num_rows = t_series.n_rows;
  arma::mat cum_sum;

  if (weight_s.isNotNull()) {
    // Copy weight_s
    arma::vec weights_vec = Rcpp::as<vec>(weight_s);
    arma::uword num_weights = weights_vec.n_elem;
    // Calculate the weighted averages as convolutions
    cum_sum = arma::conv2(t_series, weights_vec, "full");
    // Copy the warmup period
    // cout << "num_weights = " << num_weights << endl;
    cum_sum.rows(0, num_weights-2) = t_series.rows(0, num_weights-2);
    cum_sum = cum_sum.rows(0, num_rows-1);
    // cout << "cum_sum.n_rows = " << cum_sum.n_rows << endl;
  } else {
    // Calculate cumulative returns
    cum_sum = arma::cumsum(t_series, 0);
  }  // end if
  
  
  // Declare empty end points
  arma::uvec end_p;
  // Update end points
  if (end_points.isNotNull()) {
    // Copy end_points
    end_p = Rcpp::as<uvec>(end_points);
  } else if (stu_b.isNotNull()) {
    // Calculate end points with stu_b
    end_p = arma::regspace<uvec>(Rcpp::as<uword>(stu_b), look_back, num_rows + look_back);
    end_p = end_p.elem(find(end_p < num_rows));
  }  // end if
  
  
  // Calculate the rolling sums
  if (end_p.is_empty() && weight_s.isNotNull()) {
    // Do nothing
    // Return the weighted averages (convolutions) at each point
    // return cum_sum;
  } else if (end_p.is_empty() && !weight_s.isNotNull()) {
    // Return rolling sums at each point
    cum_sum = diff_it(cum_sum, look_back, true);
  } else if (!end_p.is_empty() && weight_s.isNotNull()) {
    // Return the weighted averages (convolutions) at end points
    cum_sum = cum_sum.rows(end_p);
  } else if (!end_p.is_empty() && !weight_s.isNotNull()) {
    // Return the rolling sums at end points
    cum_sum = cum_sum.rows(end_p);
    cum_sum = diff_it(cum_sum, 1, true);
  }  // end if
  
  return cum_sum;
  
}  // end roll_sum



////////////////////////////////////////////////
// Old versions

//' @export
// [[Rcpp::export]]
bool test_empty(arma::mat& re_turns,
                   arma::uword ste_p = 1,
                   Rcpp::Nullable<int> stu_b = R_NilValue, 
                   Rcpp::Nullable<Rcpp::IntegerVector> end_points = R_NilValue) {
  
  arma::uword num_rows = re_turns.n_rows;
  // Calculate cumulative returns
  // arma::mat cum_sum = arma::cumsum(re_turns, 0);
  arma::uvec end_p;
  // arma::uvec end_p = arma::ones<uvec>(1);
  // end_p.reset();
  
  // bool is_empty;
  
  if (end_points.isNotNull()) {
    // Calculate rolling sums at end_points
    end_p = Rcpp::as<uvec>(end_points);
  } else if (stu_b.isNotNull()) {
    // Calculate rolling sums at end points with stu_b
    end_p = arma::regspace<uvec>(Rcpp::as<uword>(stu_b), ste_p, num_rows + ste_p);
    end_p = end_p.elem(find(end_p < num_rows));
  }  // end if
  
  if (end_p.is_empty())
    return true;
  else  
    return false;
    
}  // end test_empty


////////////////////////////////////////////////

//' @export
// [[Rcpp::export]]
arma::mat roll_sumo(arma::mat re_turns, 
                    arma::uword look_back) {
  
  // cout << "end_p = " << end_p << endl;
  
  // Calculate cumulative returns at end points
  arma::mat cum_sum = arma::cumsum(re_turns, 0);
  
  // Return the aggregated returns
  return diff_it(cum_sum, look_back, true);
  
}  // end roll_sumo



//' @export
// [[Rcpp::export]]
arma::mat roll_retso(arma::mat& re_turns, 
                   arma::uword ste_p = 1,
                   arma::uword stu_b = 0) {

  // cout << "stu_b = " << stu_b << endl;
  // Calculate end points
  arma::uword num_rows = re_turns.n_rows;
  arma::uvec end_p = arma::regspace<uvec>(stu_b, ste_p, num_rows + ste_p);
  end_p = end_p.elem(find(end_p < num_rows));
  // cout << "end_p = " << end_p << endl;
  
  // Calculate cumulative returns at end points
  arma::mat cum_sum = arma::cumsum(re_turns, 0);
  cum_sum = cum_sum.rows(end_p);

  // Return the aggregated returns
  return diff_it(cum_sum, 1, true);
  
  //// Attempt at alternative code
  // Attempt to add to argument list: default value is empty vector - doesn't work 
  // arma::uvec end_points = arma::ones<uvec>(0), 
  // if (end_points.is_empty())
  // Calculate the variance without aggregations
  // end_points = calc_endpoints(re_turns.n_rows, ste_p);
  // This R code works with is_empty() but it requires assigning integer(0):
  // foo <- roll_retso(re_turns, end_points=integer(0), ste_p=22)
  //// End alternative code
  
}  // end roll_retso

