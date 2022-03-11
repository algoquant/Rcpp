// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/roll_mad.cpp")

arma::uvec calc_endpoints(arma::uword length, arma::uword step = 1, bool front = true) {
  
  // Calculate number of intervals that fit over length
  arma::uword num_points = length/step;
  arma::uvec endp;
  
  if (length == step*num_points) {
    // No stub interval
    // Include the first placeholder end point - legacy code
    // endp = arma::cumsum(arma::ones<uvec>(num_points));
    // endp = arma::regspace<uvec>(0, step, length);
    endp = arma::regspace<uvec>(step, step, length);
  } else {
    // Need to add stub interval
    // Include the first placeholder end point - legacy code
    // endp = arma::cumsum(arma::ones<uvec>(num_points + 1));
    // endp = arma::regspace<uvec>(0, step, length + step);
    endp = arma::regspace<uvec>(step, step, length + step);
    if (front) {
      // Stub interval at beginning
      endp = endp - step + length % step;
    } else {
      // Stub interval at end
      // The last end point must be equal to length
      endp(num_points) = length;
    }  // end if
  }  // end if
  
  // Set the first end point to zero - it's a placeholder
  // endp(0) = 0;
  // Subtract 1 from endp because indexing starts at 0
  endp = endp - 1;
  return endp;
  
}  // end calc_endpoints


arma::uvec calc_startpoints(arma::uvec endp, arma::uword look_back) {
  
  arma::uword num_points = endp.n_elem;
  arma::uvec startp = arma::join_cols(arma::zeros<uvec>(look_back), 
                                       endp.subvec(0, num_points - look_back - 1) + 1);
  
  return startp;
  
}  // end calc_startpoints


//' @export
// [[Rcpp::export]]
arma::mat roll_mad(arma::mat se_ries, 
                   arma::uword step = 1, 
                   arma::uword look_back = 11, 
                   std::string method="moment") {
  
  // Calculate end points
  arma::uword nrows = se_ries.n_rows;
  arma::uvec endp = calc_endpoints(nrows, step);
  // Start points equal to end points lagged by look_back
  arma::uvec startp = calc_startpoints(endp, look_back);
  // Allocate variance matrix
  arma::uword num_points = endp.n_elem;
  arma::mat variance = arma::zeros<mat>(num_points, se_ries.n_cols);
  arma::mat sub_series;
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_points; ep++) {
    // Calculate variance
    if (endp(ep) > startp(ep)) {
      sub_series = se_ries.rows(startp(ep), endp(ep));
      sub_series.each_row() -= arma::median(sub_series, 0);
      variance.row(ep) = arma::median(arma::abs(sub_series), 0);
    }  // end if
  }  // end for
  
  return variance;
  
}  // end roll_mad

