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



////////////////////////////////////////////////
// Test versions to be migrated to package HighFreq


////////////////////////////////////////////////////////////
//' Calculate the variance for overlapping aggregated returns: the variance of
//' returns aggregated over k time periods. a vector of end points that divides
//' a vector into equal intervals.
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
arma::rowvec calc_var(arma::mat& se_ries, 
                      arma::uword ste_p = 1) {
  
  // arma::uword num_cols = se_ries.n_cols;
  
  if (ste_p == 1)
    // Calculate the variance without aggregations
    return arma::var(se_ries);
  else {
    arma::mat cum_sum = arma::cumsum(se_ries, 0);
    // Calculate the variance of aggregated returns
    // arma::uvec end_p = calc_endpoints(se_ries.n_rows, ste_p);
    // cum_sum = cum_sum.rows(end_p);
    // return arma::var(diff_it(cum_sum, 1, false));
    
    // Perform loop over the stubs
    arma::uword num_rows = se_ries.n_rows;
    arma::mat aggs;
    arma::uvec end_p;
    arma::mat var_s(ste_p, se_ries.n_cols);
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
// Old versions

//' Calculate the end points with a stub interval passed in as an argument.
//' @export
// [[Rcpp::export]]
arma::uvec calc_endpoints_stub(arma::uword len_gth, arma::uword ste_p, arma::uword stu_b = 0) {
  
  // Calculate number of intervals that fit over len_gth
  // arma::uword num_points = len_gth/ste_p;
  arma::uvec endp = arma::regspace<uvec>(stu_b, ste_p, len_gth + ste_p);
  
  return endp.elem(find(endp <= len_gth));
  
}  // end calc_endpoints_stub


// Demonstration of using default NULL arguments in Rcpp code.
// calc_endpoints_null() implements default NULL arguments in Rcpp code.
// Calculate the end points with a stub interval.
//' @export
// [[Rcpp::export]]
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
// New version with stub argument
//' @param \code{stu_b} An \emph{integer} value equal to the first stub interval
//'   for calculating the end points.
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

////////////////////////////////////////////////



// Calculate the cumulative returns (rolling sums) at end points.
//' @export
// [[Rcpp::export]]
arma::mat roll_rets(arma::mat& se_ries, 
                   arma::uword ste_p = 1,
                   arma::uword stu_b = 0) {

  // cout << "stu_b = " << stu_b << endl;
  // Calculate end points
  arma::uword num_rows = se_ries.n_rows;
  arma::uvec end_p = arma::regspace<uvec>(stu_b, ste_p, num_rows + ste_p);
  end_p = end_p.elem(find(end_p < num_rows));
  // cout << "end_p = " << end_p << endl;
  
  // Calculate cumulative returns at end points.
  arma::mat cum_sum = arma::cumsum(se_ries, 0);
  cum_sum = cum_sum.rows(end_p);

  // Return the differences of the cumulative returns
  return diff_it(cum_sum, 1, true);
  
}  // end roll_rets

