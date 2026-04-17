// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
////////////////////////////
// RcppArmadillo functions for backtesting portfolio momentum strategies.
// 
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/run_portf.cpp")
// 
// Copyright: 2026 Jerzy Pawlowski
////////////////////////////

#include "RcppArmadillo.h"
// #include <vector>
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
// Use STL
using namespace std;
// For eigen solver SymEigsSolver
// using namespace arma::newarp;


void push_covar(const arma::rowvec& rets, // Row of new asset returns
                arma::mat& covmat,  // Covariance matrix
                arma::rowvec& retmcov, // Trailing means of the returns
                const double& lambdacov) { // Covariance decay factor
  
  double lambdac1 = 1-lambdacov;
  
  // Update the means of the returns
  retmcov = lambdacov*retmcov + lambdac1*rets;
  // Calculate the centered returns
  arma::rowvec retc = (rets - retmcov);
  
  // Update the covariance of the returns
  covmat = lambdacov*covmat + lambdac1*arma::trans(retc)*retc;
  
}  // end push_covar



arma::mat calc_invred(const arma::mat& matrixv, 
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
  eigenval = eigenval.subvec(ncols-dimax, ncols-1);
  eigenvec = eigenvec.cols(ncols-dimax, ncols-1);
  
  // Calculate the reduced inverse from the eigen decomposition
  return eigenvec*arma::diagmat(1/eigenval)*eigenvec.t();
  
}  // end calc_invred


// [[Rcpp::export]]
arma::mat sim_portfoptim(const arma::mat& rets, // Asset returns
                         const arma::uword dimax, // Number of eigen vectors for dimension reduction
                         const double lambdaf, // Returns decay factor
                         const double lambdacov, // Covariance decay factor
                         const double lambdaw, // Weight decay factor
                         const double volt = 0.01) { // Volatility target for scaling the weights
  
  arma::uword nrows = rets.n_rows;
  arma::uword ncols = rets.n_cols;
  arma::rowvec retsc; // Row of scaled asset returns
  arma::rowvec retm(ncols, fill::zeros); // EMA of asset returns
  arma::rowvec varv; // Variance of asset returns
  arma::rowvec weightv; // Row of portfolio weights
  double pnlm = 0; // EMA of strategy PnLs
  double pnlvar = 0; // Variance of strategy PnLs
  arma::colvec varportf(nrows); // Portfolio variance
  arma::mat invmat(ncols, ncols, fill::ones); // Inverse covariance matrix
  arma::mat weightm(nrows, ncols, fill::zeros); // Portfolio weights
  arma::colvec pnls(nrows, fill::zeros); // Strategy PnLs
  double weightd; // Sum of squared weights
  double lambda1 = 1-lambdaf;
  double lambdac1 = 1-lambdacov;
  double lambdaw1 = 1-lambdaw;
  
  // Initialize EMA of returns with first row of returns
  arma::rowvec retmcov = rets.row(0); // EMA of returns for covariance
  // Calculate the centered returns
  arma::rowvec retc = retmcov;
  // Initialize covariance matrix to unit diagonal matrix - doesn't perform well because takes too long to warm up
  // arma::mat covmat = arma::eye(ncols, ncols);
  // Initialize covariance matrix to outer product of the returns - initial matrix is singular but it warms up faster than the unit matrix
  // arma::mat covmat = arma::trans(retmcov)*retmcov;
  // Initialize covariance matrix to outer product of the returns plus unit diagonal matrix to prevent singularity
  arma::mat covmat = arma::trans(retmcov)*retmcov + arma::dot(retmcov, retmcov)*arma::eye(ncols, ncols);
  
  // Perform loop over the rows (time)
  for (arma::uword it = 1; it < nrows; it++) {
    
    // retsc = rets.row(it);
    // Calculate the variance of the returns using the diagonal of the covariance matrix
    varv = arma::trans(covmat.diag());
    varv = varv + 1e-8; // Add small number to prevent division by zero
    // varv.replace(0, 1);
    // Scale the returns by the trailing volatility
    retsc = rets.row(it)/arma::sqrt(varv);
    // weightm.row(it) = weightm.row(it)*volt/weightd;
    // Calculate the strategy PnLs - the products of the lagged weights times the asset returns
    pnls(it) = arma::dot(retsc, weightm.row(it-1));
    // pnlvar = lambdaf*pnlvar + lambda1*pow(pnls(it) - pnlm, 2);
    // pnlm = lambdaf*pnlm + lambda1*pnls(it);
    // Scale the weights using the volatility target
    // weightd = std::sqrt(pnlvar);
    // if (weightd == 0) weightd = 1;
    // pnls(it) = pnls(it)*volt/weightd;
    // Update the covariance means of the returns
    retmcov = lambdacov*retmcov + lambdac1*retsc;
    // Calculate the centered returns
    retc = (retsc - retmcov);
    // Update the covariance of the returns
    covmat = lambdacov*covmat + lambdac1*arma::trans(retc)*retc;
    // Calculate the trailing EMA of asset returns
    retm = lambdaf*retm + lambda1*rets.row(it);
    // Calculate the reduced inverse of the covariance matrix using its eigen decomposition
    invmat = calc_invred(covmat, dimax);
    // Calculate the weights using recursive inverse
    // calc_invrec(covmat, invmat);
    // weightm.row(it) = retm*invmat;
    // Calculate the EMA weights
    weightv = retm*invmat;
    // Scale the weights to the volatility target
    // varportf(it) = arma::as_scalar(weightv*covmat*weightv.t());
    // weightv = weightv*volt/std::sqrt(varportf(it));
    // Scale the weights so their sum of squares is equal to one
    // weightd = std::sqrt(arma::sum(arma::square(weightv)));
    // weightv = weightv/weightd;
    // cout << "it: " << it << endl << "weightv: " << endl << weightv << endl;
    // Calculate the portfolio variance using the weights
    weightm.row(it) = lambdaw*weightm.row(it-1) + lambdaw1*weightv;
    // Scale the weights so their sum of squares is equal to one
    weightd = std::sqrt(arma::sum(arma::square(weightm.row(it))));
    // if (weightd == 0) weightd = 1;
    weightm.row(it) = weightm.row(it)/weightd;
    // Scale the weights to the volatility target
    // varportf(it) = arma::as_scalar(weightm.row(it)*covmat*weightm.row(it).t());
    // weightm.row(it) = weightm.row(it)*volt/std::sqrt(varportf(it));
  }  // end for
  
  // Return the strategy PnLs and the portfolio weights
  // return arma::join_rows(pnls, weightm, varportf);
  return arma::join_rows(pnls, weightm);
  
}  // end sim_portfoptim




