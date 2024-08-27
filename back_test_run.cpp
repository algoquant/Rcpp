////////////////////////////
// Functions to backtest trading strategies.
// Contains a copy of the function HighFreq::calc_inv().
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
arma::mat calc_inv(const arma::mat& matrixv, 
                   arma::uword dimax = 0, // Number of eigen vectors for dimension reduction
                   double singmin = 0.0) { // Threshold for discarding small eigen values
  
  // Allocate Eigen variables
  arma::uword ncols = matrixv.n_cols;
  arma::vec eigenval(ncols, fill::zeros); // Eigen values
  arma::mat eigenvec(ncols, ncols, fill::zeros); // Eigen vectors
  // Calculate the eigen decomposition - the eigenvalues are in ascending order
  arma::eig_sym(eigenval, eigenvec, matrixv);
  // Calculate the number of non-small singular values
  arma::uword neigen = arma::sum(eigenval > singmin*arma::sum(eigenval));
  
  // If no regularization then set dimax to neigen
  if (dimax == 0) {
    // Set dimax
    dimax = neigen;
  } else {
    // Adjust dimax
    dimax = std::min(dimax, neigen);
  }  // end if
  
  // Remove all small singular values
  // cout << "neigen: " << neigen << endl;
  // cout << "eigenval: " << eigenval << endl;
  // cout << "ncols: " << ncols << endl;
  // cout << "ncols-dimax: " << (ncols-dimax) << endl;
  eigenval = eigenval.subvec(ncols-dimax, ncols-1);
  eigenvec = eigenvec.cols(ncols-dimax, ncols-1);
  
  // Calculate the reduced inverse from the eigen decomposition
  return eigenvec*arma::diagmat(1/eigenval)*eigenvec.t();
  
}  // end calc_inv



//' Backtest a trading strategy using the exponential 
//' moving average of the covariance matrix.
//' @export
// [[Rcpp::export]]
arma::mat back_testp(const arma::mat& retp, // Asset returns
                    double lambda = 0.0, // Decay factor for averaging the portfolio weights
                    double singmin = 1e-5, // Threshold level for discarding small singular values
                    arma::uword dimax = 0, // Regularization intensity
                    bool scale = true, // Scale the weights
                    double voltarget = 0.01,
                    double coeff = 1.0,
                    double bidask = 0.0) {
  
  double lambda1 = 1-lambda;
  arma::uword nweights = retp.n_cols;
  arma::uword nrows = retp.n_rows;
  arma::mat reti = retp.row(0);
  arma::mat retit = reti.t();
  arma::mat meanv = retit;
  // arma::mat rets2 = arma::square(retp);
  // arma::mat varv = rets2.row(0);
  arma::vec weightv(nweights, fill::zeros);
  // arma::vec weightm(nweights, fill::zeros);
  // arma::vec weights_past = ones(nweights)/sqrt(nweights);
  arma::mat covmat = retp.row(0).t()*retp.row(0);
  arma::mat pnls = zeros(nrows, 1);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < nrows; it++) {
    reti = retp.row(it);
    retit = reti.t();
    // cout << "it: " << it << endl;
    // arma::mat invmat = calc_inv(cov(retp), dimax);
    // Calculate out-of-sample PnLs
    // weightm = lambda1*weightv + lambda*weightm;
    pnls.row(it) = reti*weightv;
    // Add transaction costs
    // pnls.row(endp(it-1)+1) -= bidask*sum(abs(weightv - weights_past))/2;
    // Update the means using the decay factor
    meanv = lambda1*retit + lambda*meanv;
    // Calculate the weights
    // Portfolio optimization
    // Update the covariance matrix using the decay factor
    covmat = lambda1*retit*reti + lambda*covmat;
    // cout << "covmat: " << covmat << endl;
    weightv = calc_inv(covmat, dimax, singmin)*meanv;
    // Momentum
    // varv = lambda1*(rets2.row(it) - arma::square(meanv.t())) + lambda*varv;
    // varv = lambda1*rets2.row(it) + lambda*varv;
    // varv.replace(0, 1);
    // weightv = conv_to<vec>::from(meanv.t()/varv);
    // weightv = (weightv - arma::mean(weightv));
    // weightv = weightv/std::sqrt(sum(arma::square(weightv)));
  }  // end for
  
  // Return the strategy pnls
  return pnls;
  
}  // end back_testp



//' @export
// [[Rcpp::export]]
arma::mat back_testx(const arma::mat& retp, // Asset returns
                     const arma::rowvec& betav, // Asset betas
                     double lambda = 0.0, // Decay factor for averaging the portfolio weights
                     double singmin = 1e-5, // Threshold level for discarding small singular values
                     arma::uword dimax = 0, // Regularization intensity
                     bool scale = true, // Scale the weights
                     double voltarget = 0.01,
                     double coeff = 1.0,
                     double bidask = 0.0) {
  
  double lambda1 = 1-lambda;
  arma::uword nweights = retp.n_cols;
  arma::uword nrows = retp.n_rows;
  arma::mat reti = retp.row(0);
  arma::mat retit = reti.t();
  arma::mat meanv = retit;
  // arma::mat rets2 = arma::square(retp);
  // arma::mat varv = rets2.row(0);
  arma::vec weightv(nweights, fill::zeros);
  // arma::vec weightm(nweights, fill::zeros);
  // arma::vec weights_past = ones(nweights)/sqrt(nweights);
  arma::mat covmat = retp.row(0).t()*retp.row(0);
  arma::mat pnls = zeros(nrows, 1);
  arma::mat betas = zeros(nrows, 1);
  
  // cout << "reti: " << reti << endl;
  // cout << "betav: " << betav << endl;
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < nrows; it++) {
    reti = retp.row(it);
    retit = reti.t();
    // cout << "it: " << it << endl;
    // arma::mat invmat = calc_inv(cov(retp), dimax);
    // Calculate out-of-sample PnLs
    // weightm = lambda1*weightv + lambda*weightm;
    pnls.row(it) = reti*weightv;
    betas.row(it) = betav*weightv;
    // Add transaction costs
    // pnls.row(endp(it-1)+1) -= bidask*sum(abs(weightv - weights_past))/2;
    // Update the means using the decay factor
    meanv = lambda1*retit + lambda*meanv;
    // Calculate the weights
    // Portfolio optimization
    // Update the covariance matrix using the decay factor
    covmat = lambda1*retit*reti + lambda*covmat;
    // cout << "covmat: " << covmat << endl;
    weightv = calc_inv(covmat, dimax, singmin)*meanv;
    // Momentum
    // varv = lambda1*(rets2.row(it) - arma::square(meanv.t())) + lambda*varv;
    // varv = lambda1*rets2.row(it) + lambda*varv;
    // varv.replace(0, 1);
    // weightv = conv_to<vec>::from(meanv.t()/varv);
    // weightv = (weightv - arma::mean(weightv));
    weightv = weightv/std::sqrt(sum(arma::square(weightv)));
  }  // end for
  
  // Return the strategy pnls
  // return pnls;
  return arma::join_rows(pnls, betas);
  
}  // end back_testx

