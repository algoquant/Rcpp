// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
////////////////////////////
// RcppArmadillo functions for matrix inverse and pseudo-inverse.
// 
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/calc_inv.cpp")
// 
// Copyright: 2026 Jerzy Pawlowski
////////////////////////////

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

// using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;


////////////////////////////////////////////////////////////
//' Calculate the inverse of a square positive definite \emph{matrix} using
//' \code{RcppArmadillo}.
//' 
//' @param \code{matrixv} A square positive definite \emph{matrix}.
//'
//' @return A \emph{matrix} equal to the inverse of \code{matrixv}.
//'
//' @details
//'   The function \code{calc_inv()} calculates the inverse of the matrix
//'   \code{matrixv} using the \code{RcppArmadillo} function \code{arma::inv()}.
//'   
//'   The product of a matrix \eqn{\strong{C}} times its inverse
//'   \eqn{\strong{C}^{-1}} is equal to the identity matrix \eqn{\strong{I}}:
//'   \deqn{
//'     \strong{C} \, \strong{C}^{-1} = \strong{I}
//'   }
//'   
//'   The matrix \code{matrixv} must be square and positive definite (non-singular).
//'   For singular or near-singular matrices, consider using \code{calc_ginv()}
//'   which computes the Moore-Penrose pseudo-inverse.
//'   But \code{calc_ginv()} can be several times slower than \code{calc_inv()}
//'   for large matrices.
//'   For symmetric positive definite matrices, consider using 
//'   \code{calc_inv_sympd()} which is slightly faster.
//'   
//' @examples
//' \dontrun{
//' # Create a random positive definite matrix
//' matv <- matrix(rnorm(25), nc=5)
//' matv <- t(matv) %*% matv
//' # Add small identity matrix to ensure positive definiteness
//' matv <- matv + 1e-2*diag(5)
//' # Calculate matrix inverse using RcppArmadillo
//' invmat <- HighFreq::calc_inv(matv)
//' # Calculate matrix inverse in R
//' invr <- solve(matv)
//' all.equal(invmat, invr, check.attributes=FALSE)
//' }  # end dontrun
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& matrixv) {
  
  // Calculate the reduced inverse from the eigen decomposition
  return arma::inv(matrixv);
  
}  // end calc_inv

 
 ////////////////////////////////////////////////////////////
 //' Calculate the generalized inverse (Moore-Penrose pseudo-inverse) of a 
 //' \emph{matrix} using \code{RcppArmadillo}.
 //' 
 //' @param \code{matrixv} A \emph{matrix} (can be rectangular or singular).
 //'
 //' @return A \emph{matrix} equal to the pseudo-inverse of \code{matrixv}.
 //'
 //' @details
 //'   The function \code{calc_ginv()} calculates the generalized inverse
 //'   (Moore-Penrose pseudo-inverse) of the matrix \code{matrixv} using the
 //'   \code{RcppArmadillo} function \code{arma::pinv()}.
 //'   
 //'   The generalized inverse \eqn{\strong{C}^{+}} of a matrix \eqn{\strong{C}} 
 //'   satisfies the Moore-Penrose condition:
 //'   \deqn{
 //'     \strong{C} \, \strong{C}^{+} \, \strong{C} = \strong{C}
 //'   }
 //'   Which is a generalization of the standard inverse condition: 
 //'   \deqn{
 //'     \strong{C} \, \strong{C}^{-1} = \strong{I}
 //'   }
 //'   
 //'   Unlike the standard inverse, the generalized inverse exists for any matrix,
 //'   including rectangular and singular matrices. For non-singular square 
 //'   matrices, the generalized inverse equals the standard inverse.
 //'   
 //'   The generalized inverse is computed using singular value decomposition (SVD),
 //'   which makes it more robust but several times slower than \code{calc_inv()}.
 //'   For large non-singular matrices, \code{calc_inv()} is preferred for speed.
 //'   
 //' @examples
 //' \dontrun{
 //' # Create a singular matrix (rank deficient)
 //' matv <- matrix(rnorm(25), nc=5)
 //' matv[, 5] <- matv[, 1]  # Make it singular
 //' # Calculate generalized inverse using RcppArmadillo
 //' invmat <- HighFreq::calc_ginv(matv)
 //' # Calculate generalized inverse using MASS package in R
 //' invr <- MASS::ginv(matv)
 //' all.equal(invmat, invr, check.attributes=FALSE)
 //' # Verify Moore-Penrose condition
 //' all.equal(matv %*% invmat %*% matv, matv, check.attributes=FALSE)
 //' }  # end dontrun
 //' 
 //' @export
 // [[Rcpp::export]]
 arma::mat calc_ginv(const arma::mat& matrixv) {
   
   // Calculate the generalized inverse (Moore-Penrose pseudo-inverse)
   return arma::pinv(matrixv);
   
 }  // end calc_ginv
 
 
 
 // [[Rcpp::export]]
 arma::mat calc_inv_sympd(const arma::mat& matrixv) { // Threshold for discarding small eigen values
   
   // Calculate the reduced inverse from the eigen decomposition
   return arma::inv_sympd(matrixv);
   
 }  // end calc_inv_sympd
 
 
 
 