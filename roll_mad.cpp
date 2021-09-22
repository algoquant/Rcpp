// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/roll_mad.cpp")

arma::uvec calc_endpoints(arma::uword length, arma::uword step = 1, bool front = true) {
  
  // Calculate number of intervals that fit over length
  arma::uword num_points = length/step;
  arma::uvec end_p;
  
  if (length == step*num_points) {
    // No stub interval
    // Include the first placeholder end point - legacy code
    // end_p = arma::cumsum(arma::ones<uvec>(num_points));
    // end_p = arma::regspace<uvec>(0, step, length);
    end_p = arma::regspace<uvec>(step, step, length);
  } else {
    // Need to add stub interval
    // Include the first placeholder end point - legacy code
    // end_p = arma::cumsum(arma::ones<uvec>(num_points + 1));
    // end_p = arma::regspace<uvec>(0, step, length + step);
    end_p = arma::regspace<uvec>(step, step, length + step);
    if (front) {
      // Stub interval at beginning
      end_p = end_p - step + length % step;
    } else {
      // Stub interval at end
      // The last end point must be equal to length
      end_p(num_points) = length;
    }  // end if
  }  // end if
  
  // Set the first end point to zero - it's a placeholder
  // end_p(0) = 0;
  // Subtract 1 from end_p because indexing starts at 0
  end_p = end_p - 1;
  return end_p;
  
}  // end calc_endpoints


arma::uvec calc_startpoints(arma::uvec end_p, arma::uword look_back) {
  
  arma::uword num_points = end_p.n_elem;
  arma::uvec start_p = arma::join_cols(arma::zeros<uvec>(look_back), 
                                       end_p.subvec(0, num_points - look_back - 1) + 1);
  
  return start_p;
  
}  // end calc_startpoints


//' @export
// [[Rcpp::export]]
arma::mat roll_mad(arma::mat se_ries, 
                   arma::uword step = 1, 
                   arma::uword look_back = 11, 
                   std::string method="moment") {
  
  // Calculate end points
  arma::uword num_rows = se_ries.n_rows;
  arma::uvec end_p = calc_endpoints(num_rows, step);
  // Start points equal to end points lagged by look_back
  arma::uvec start_p = calc_startpoints(end_p, look_back);
  // Allocate variance matrix
  arma::uword num_points = end_p.n_elem;
  arma::mat vari_ance = arma::zeros<mat>(num_points, se_ries.n_cols);
  arma::mat sub_series;
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_points; ep++) {
    // Calculate variance
    if (end_p(ep) > start_p(ep)) {
      sub_series = se_ries.rows(start_p(ep), end_p(ep));
      sub_series.each_row() -= arma::median(sub_series, 0);
      vari_ance.row(ep) = arma::median(arma::abs(sub_series), 0);
    }  // end if
  }  // end for
  
  return vari_ance;
  
}  // end roll_mad

