////////////////////////////
// Copies of functions HighFreq::calc_weights() and 
// HighFreq::back_test() for testing.
// 
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/back_test_run.cpp")

#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword eigen_max = 0) {
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD
  arma::svd(svdu, svdval, svdv, tseries);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > eigen_thresh*arma::sum(svdval));
  
  if (eigen_max == 0) {
    // Set eigen_max
    eigen_max = svdnum - 1;
  } else {
    // Adjust eigen_max
    eigen_max = std::min(eigen_max - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, eigen_max);
  svdu = svdu.cols(0, eigen_max);
  svdv = svdv.cols(0, eigen_max);
  
  // Calculate the regularized inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_inv



//' @export
// [[Rcpp::export]]
arma::mat back_test(const arma::mat& returns, // Asset returns
                    double lambda = 0.0, // Decay factor for averaging the portfolio weights
                    double eigen_thresh = 1e-5, // Threshold level for discarding small singular values
                    arma::uword eigen_max = 0, // Regularization intensity
                    bool scale = true, // Scale the weights
                    double vol_target = 0.01,
                    double coeff = 1.0,
                    double bid_offer = 0.0) {
  
  double lambda1 = 1-lambda;
  arma::uword nweights = returns.n_cols;
  arma::uword nrows = returns.n_rows;
  arma::mat retsrow = returns.row(0);
  arma::mat retsrowt = returns.row(0).t();
  arma::mat means = returns.row(0).t();
  arma::mat rets2 = arma::square(returns);
  arma::mat vars = rets2.row(0);
  arma::vec weights(nweights, fill::zeros);
  // arma::vec weights_past = ones(nweights)/sqrt(nweights);
  arma::mat covmat = returns.row(0).t()*returns.row(0);
  arma::mat pnls = zeros(nrows, 1);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < nrows; it++) {
    retsrow = returns.row(it);
    retsrowt = returns.row(it).t();
    // cout << "it: " << it << endl;
    // arma::mat invmat = calc_inv(cov(returns), eigen_max);
    // Calculate out-of-sample returns
    pnls.row(it) = retsrow*weights;
    // Add transaction costs
    // pnls.row(endp(it-1)+1) -= bid_offer*sum(abs(weights - weights_past))/2;
    // Update the means using the decay factor
    means = lambda1*retsrowt + lambda*means;
    // Calculate the weights
    // Portfolio optimization
    // Update the covariance matrix using the decay factor
    covmat = lambda1*retsrowt*retsrow + lambda*covmat;
    weights = calc_inv(covmat, eigen_thresh, eigen_max)*means;
    // Momentum
    // vars = lambda1*(rets2.row(it) - arma::square(means.t())) + lambda*vars;
    // vars = lambda1*rets2.row(it) + lambda*vars;
    // vars.replace(0, 1);
    // weights = conv_to<vec>::from(means.t()/vars);
    // weights = (weights - arma::mean(weights));
    weights = weights/std::sqrt(sum(arma::square(weights)));
  }  // end for
  
  // Return the strategy pnls
  return pnls;
  
}  // end back_test

