#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/sim_boll.cpp")

// The function sim_boll() calculates the positions of 
// a Bollinger strategy, using the z-scores.
// It doubles down and adds to the position when the 
// z-score increases in the same direction.
// It holds its position as long as the z-score doesn't 
// change sign.
// It flips the position when the z-score changes sign.
// The strategy is path dependent and requires performing 
// a loop over time.

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec sim_boll(const arma::vec& zscores) {
    int nrows = zscores.n_elem;
    arma::vec posv(nrows, fill::zeros);  // Stock positions, initialized to 0
    
    for (int it = 1; it < nrows; ++it) {
      // If the z-score increased then increase the position size
      // If the z-score changed sign then flip the position
      if ((std::abs(zscores(it-1)) > std::abs(posv(it-1))) || 
            ((zscores(it-1) * posv(it-1)) > 0)) {
            // Increase position size
            posv(it) = -zscores(it-1);
        } else {
            // Do nothing - maintain previous position
            posv(it) = posv(it-1);
        }
    }
    
    return posv;
}
