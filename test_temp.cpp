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
arma::vec lag_vec(arma::vec& vec_tor, int lagg=1, bool pad_zeros=true) {
  
  int len_gth = (vec_tor.n_elem-1);
  
  if (lagg > 0) {
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros(lagg), 
                             vec_tor.subvec(0, len_gth-lagg));
    } else {
      // Pad front with first element of vec_tor
      return arma::join_cols(arma::repelem(vec_tor.subvec(0, 0), lagg, 1), 
                             vec_tor.subvec(0, len_gth-lagg));
    }  // end if
  } else {
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(vec_tor.subvec(-lagg, len_gth), 
                             arma::zeros(-lagg));
    } else {
      // Pad back with last element of vec_tor
      return arma::join_cols(vec_tor.subvec(-lagg, len_gth), 
                             arma::repelem(vec_tor.subvec(len_gth, len_gth), -lagg, 1));
    }  // end if
  }  // end if
  
}  // end lag_vec


//' @export
// [[Rcpp::export]]
arma::mat lag_it(arma::mat& t_series, int lagg=1) {
  int num_rows = (t_series.n_rows-1);
  
  if (lagg > 0)
    return arma::join_cols(arma::repelem(t_series.row(0), lagg, 1), 
                           t_series.rows(0, num_rows-lagg));
  else
    return arma::join_cols(t_series.rows(-lagg, num_rows), 
                           arma::repelem(t_series.row(num_rows), -lagg, 1));
  
}  // end lag_it



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
arma::uvec calc_endpoints_stub(arma::uword len_gth, arma::uword ste_p, arma::uword stub = 0) {
  
  // Calculate number of intervals that fit over len_gth
  // arma::uword num_points = len_gth/ste_p;
  arma::uvec endp = arma::regspace<uvec>(stub, ste_p, len_gth + ste_p);
  
  return endp.elem(find(endp <= len_gth));
  
}  // end calc_endpoints_stub



//' @export
// [[Rcpp::export]]
arma::mat diff_it(arma::mat& mat_rix, int lagg=1, bool padd=true) {
  
  int num_rows = (mat_rix.n_rows-1);
  // Matrix difference without padding
  arma::mat diff_mat = (mat_rix.rows(lagg, num_rows) - mat_rix.rows(0, num_rows - lagg));
  
  if (padd)
    // Pad diff_mat with zeros at the beginning
    return arma::join_cols(arma::zeros(lagg, mat_rix.n_cols), diff_mat);
  else
    // Don't pad the output
    return diff_mat;
  
}  // end diff_it




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
  
  // int num_cols = re_turns.n_cols;
  
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
    for (arma::uword stub=1; stub <= ste_p; stub++) {
      end_p = arma::regspace<uvec>(stub, ste_p, num_rows + ste_p);
      end_p = end_p.elem(find(end_p <= num_rows));
      aggs = cum_sum.rows(end_p-1);
      var_s.row(stub-1) = arma::var(diff_it(aggs, 1, false));
    }  // end for
    return mean(var_s);
  }  // end if
  
}  // end calc_var



//' @export
// [[Rcpp::export]]
arma::mat roll_var(arma::mat& t_series, arma::uword look_back = 1, arma::uword ste_p = 1) {
  
  // Calculate end points
  arma::uword num_rows = t_series.n_rows;
  arma::uvec end_p = calc_endpoints(num_rows, ste_p);
  // Allocate variance matrix
  arma::uword num_points = end_p.n_elem;
  arma::mat vari_ance = arma::zeros(num_points, t_series.n_cols);
  // Number of time periods in lookback interval
  // arma::uword num_steps = look_back*ste_p;
  // Calculate start points equal to end points lagged by num_steps
  // arma::uvec start_p = end_p - num_steps + 1;
  // Another better way:
  // arma::uvec start_p = arma::join_cols(arma::zeros<uvec>(look_back), 
  //                                      end_p.subvec(0, num_points - look_back - 1) + 1);
  arma::uvec start_p = calc_startpoints(end_p, look_back);
  // Make start points zero if less than zero - start_p is unsigned integer!
  // start_p.elem(find(start_p > num_rows)).zeros();

  // Perform loop over the end_points
  for (arma::uword it = 0; it < num_points; it++) {
    // Calculate variance
    vari_ance.row(it) = arma::var(t_series.rows(start_p(it), end_p(it)));
  }  // end for
  
  return vari_ance;

}  // end roll_var


//' @export
// [[Rcpp::export]]
double calc_var_ohlc(arma::mat& oh_lc, 
                     const std::string& calc_method="yang_zhang", 
                     arma::colvec lag_close=0, 
                     arma::colvec in_dex=0, 
                     const bool& scal_e=true) {
  
  int num_rows = oh_lc.n_rows;
  double co_eff = 0.34/(1.34 + (num_rows+1)/(num_rows-1));
  
  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
    // cout << "oh_lc.n_rows = " << num_rows << endl;
    // cout << "in_dex.n_rows = " << in_dex.n_rows << endl;
  }  // end if
  
  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec open_close(clo_se.n_rows);
  if (lag_close.n_rows == 1) {
    open_close = arma::join_cols(clo_se.subvec(0, 0), clo_se.subvec(0, clo_se.n_elem-2));
    open_close = (oh_lc.col(0) - open_close)/in_dex;
  } else {
    open_close = (oh_lc.col(0) - lag_close)/in_dex;
  }  // end if
  arma::colvec close_open = (clo_se - oh_lc.col(0))/in_dex;
  arma::colvec close_high = (clo_se - oh_lc.col(1))/in_dex;
  arma::colvec close_low = (clo_se - oh_lc.col(2))/in_dex;
  arma::colvec high_low = (oh_lc.col(1) - oh_lc.col(2))/in_dex;
  arma::colvec high_open = (oh_lc.col(1) - oh_lc.col(0))/in_dex;
  arma::colvec low_open = (oh_lc.col(2) - oh_lc.col(0))/in_dex;
  
  if (calc_method == "close") {
    // cout << "Calc method is Close" << endl;
    return arma::var(arma::diff(clo_se));
  } else if (calc_method == "rogers_satchell") {
    // cout << "Calc method is Rogers-Satchell" << endl;
    return -(arma::dot(close_high, high_open) +
             arma::dot(close_low, low_open))/num_rows;
  } else if (calc_method == "garman_klass") {
    // cout << "Calc method is Garman-Klass" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows;
  } else if (calc_method == "garman_klass_yz") {
    // cout << "Calc method is Garman-Klass-YZ" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows + 
            arma::var(open_close);
  } else if (calc_method == "yang_zhang") {
    // cout << "Calc method is Yang-Zhang" << endl;
    return arma::var(open_close) + co_eff*arma::var(close_open) +
      (co_eff-1)*(arma::dot(close_high, high_open) + 
      arma::dot(close_low, low_open))/num_rows;
  } else {
    cout << "Wrong calc method!" << endl;
    return 1;
  }  // end if
  
  // cout << "Calc method is " << calc_method << endl;
  
}  // end calc_var_ohlc



//' @export
// [[Rcpp::export]]
arma::vec roll_var_ohlc(arma::mat& oh_lc, 
                        arma::uword ste_p = 1, 
                        arma::uword look_back = 1, 
                        const std::string& calc_method = "yang_zhang", 
                        arma::colvec in_dex = 0, 
                        const bool& scal_e = true) {
  
  // Calculate end points
  arma::uword num_rows = oh_lc.n_rows;
  arma::uvec end_p = calc_endpoints(num_rows, ste_p);
  // Start points equal to end points lagged by look_back
  arma::uvec start_p = calc_startpoints(end_p, look_back);
  // Allocate variance matrix
  arma::uword num_points = end_p.n_elem;
  arma::vec vari_ance = arma::zeros(num_points);
  
  // Extract OHLC close prices
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec lag_close = lag_it(clo_se);
  
  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
  }  // end if
  
  // Define data subsets
  arma::mat sub_ohlc;
  arma::colvec sub_close;
  arma::colvec sub_index;
  
  // cout << "Before subset start_p(0) = " << start_p(0) << endl;
  // cout << "Before subset end_p(0) = " << end_p(0) << endl;
  // sub_ohlc = oh_lc.rows(start_p(it), end_p(it));
  // sub_close = lag_close.rows(start_p(it), end_p(it));
  // sub_index = in_dex.rows(start_p(it), end_p(it));
  
  // cout << "Before loop num_rows = " << num_rows << endl;
  // Perform loop over the end_points
  for (arma::uword it = 0; it < num_points; it++) {
    if (end_p(it) > start_p(it)) {
      // cout << "Inside loop it = " << it << endl;
      sub_ohlc = oh_lc.rows(start_p(it), end_p(it));
      sub_close = lag_close.rows(start_p(it), end_p(it));
      sub_index = in_dex.rows(start_p(it), end_p(it));
      // Calculate variance
      vari_ance(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
    }  // end if
  }  // end for
  
  // Old code below
  
  // Warmup period
  // for (arma::uword it = 1; it < look_back; it++) {
  //   arma::mat sub_ohlc = oh_lc.rows(0, it);
  //   arma::colvec sub_close = lag_close.rows(0, it);
  //   arma::colvec sub_index = in_dex.subvec(0, it);
  //   vari_ance(it, 1) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = look_back; it < num_rows; it++) {
  //   arma::mat sub_ohlc = oh_lc.rows(it-look_back+1, it);
  //   arma::colvec sub_close = lag_close.rows(it-look_back+1, it);
  //   arma::colvec sub_index = in_dex.subvec(it-look_back+1, it);
  //   vari_ance(it, 1) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
  // }  // end for
  
  return vari_ance;
  
}  // end roll_var_ohlc


////////////////////////////////////////

//' @export
// [[Rcpp::export]]
arma::mat roll_var_ohlcp(arma::mat& oh_lc, 
                        const std::string& calc_method="yang_zhang", 
                        arma::colvec in_dex=0, 
                        const bool& scal_e=true, 
                        arma::uword look_back=11) {
  
  arma::uword num_rows = oh_lc.n_rows;
  arma::mat var_vec = arma::zeros(num_rows, 2);
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec lag_close = lag_it(clo_se);
  
  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
  }  // end if
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    arma::mat sub_ohlc = oh_lc.rows(0, it);
    arma::colvec sub_close = lag_close.rows(0, it);
    arma::colvec sub_index = in_dex.subvec(0, it);
    var_vec(it, 0) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < num_rows; it++) {
    arma::mat sub_ohlc = oh_lc.rows(it-look_back+1, it);
    arma::colvec sub_close = lag_close.rows(it-look_back+1, it);
    arma::colvec sub_index = in_dex.subvec(it-look_back+1, it);
    var_vec(it, 0) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
  }  // end for

  // New code
  arma::uvec end_p = calc_endpoints(num_rows, 1);
  arma::uvec start_p = calc_startpoints(end_p, look_back);
  arma::uword num_points = end_p.n_elem;
  arma::mat sub_ohlc;
  arma::colvec sub_close;
  arma::colvec sub_index;
  for (arma::uword it = 0; it < num_points; it++) {
    if (end_p(it) > start_p(it)) {
      // cout << "Inside loop it = " << it << endl;
      sub_ohlc = oh_lc.rows(start_p(it), end_p(it));
      sub_close = lag_close.rows(start_p(it), end_p(it));
      sub_index = in_dex.rows(start_p(it), end_p(it));
      // Calculate variance
      var_vec(it, 1) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
    }
  }  // end for
  
  return var_vec;
  
}  // end roll_var_ohlcp

////////////////////////////////////////


//' @export
// [[Rcpp::export]]
arma::mat roll_points(arma::mat& oh_lc, 
                        arma::uword ste_p = 1, 
                        arma::uword look_back = 1, 
                        const std::string& calc_method = "yang_zhang", 
                        arma::colvec in_dex = 0, 
                        const bool& scal_e = true) {
  
  // Calculate end points
  arma::uword num_rows = oh_lc.n_rows;
  arma::uvec end_p = calc_endpoints(num_rows, ste_p);
  // Start points equal to end points lagged by look_back
  arma::uvec start_p = calc_startpoints(end_p, look_back);
  // Allocate variance matrix
  arma::uword num_points = end_p.n_elem;
  arma::mat vari_ance = arma::zeros(num_points, 2);
  
  // cout << "Before subset start_p(0) = " << start_p(0) << endl;
  // cout << "Before subset end_p(0) = " << end_p(0) << endl;
  // sub_ohlc = oh_lc.rows(start_p(it), end_p(it));
  // sub_close = lag_close.rows(start_p(it), end_p(it));
  // sub_index = in_dex.rows(start_p(it), end_p(it));
  
  // cout << "Before loop num_rows = " << num_rows << endl;
  // Perform loop over the end_points
  for (arma::uword it = 0; it < num_points; it++) {
    if (end_p(it) > start_p(it)) {
      // cout << "Inside loop it = " << it << endl;
      // Calculate variance
      // vari_ance(it, 0) = start_p(it);
      // vari_ance(it, 1) = end_p(it);
      vari_ance(it, 0) = mean(mean(oh_lc.rows(start_p(it), end_p(it))));
    }
  }  // end for
  
  // Old code below
  
  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    // vari_ance(it, 2) = 0;
    // vari_ance(it, 3) = it;
    vari_ance(it, 1) = mean(mean(oh_lc.rows(0, it)));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < num_rows; it++) {
    // vari_ance(it, 2) = it-look_back+1;
    // vari_ance(it, 3) = it;
    vari_ance(it, 1) = mean(mean(oh_lc.rows(it-look_back+1, it)));
  }  // end for
  
  return vari_ance;
  
}  // end roll_points



