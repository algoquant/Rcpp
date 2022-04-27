////////////////////////////
// Functions to test C++ syntax with Armadillo
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_temp.cpp")

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
void mult_rows(arma::vec& weightv, arma::mat& matrixv) {
  
  matrixv.each_row() %= weightv.t();
  
  // matrixv.each_row() %= weightv.t();
  // return matrixv;
  
}  // end mult_rows



////////////////////////////////////////////////////////////
//' Multiply the rows or columns of a \emph{matrix} times a \emph{vector},
//' element-wise and in place (without copying).
//' 
//' @param \code{vector} A \emph{numeric} \emph{vector}.
//' 
//' @param \code{matrix} A \emph{numeric} \emph{matrix}.
//' 
//' @param \code{byrow} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the rows of \code{matrix} by \code{vector}, otherwise multiply the columns
//'   (the default is \code{byrow = TRUE}.)
//' 
//' @return Void (no return value).
//' 
//' @details
//'   The function \code{mult_mat_ref()} multiplies the rows or columns of a
//'   \emph{matrix} times a \emph{vector}, element-wise and in place (without
//'   copying).
//'
//'   It accepts a \emph{pointer} to the argument \code{matrix}, and replaces
//'   the old \code{matrix} values with the new values. It performs the
//'   calculation in place, without copying the \emph{matrix} in memory, which
//'   can significantly increase the computation speed for large matrices.
//'
//'   If \code{byrow = TRUE} (the default), then function \code{mult_mat_ref()}
//'   multiplies the rows of the argument \code{matrix} times the argument
//'   \code{vector}.
//'   Otherwise it multiplies the columns of \code{matrix}.
//' 
//'   In \code{R}, \emph{matrix} multiplication is performed by columns.
//'   Performing multiplication by rows is often required, for example when
//'   multiplying stock returns by portfolio weights.
//'   But performing multiplication by rows requires explicit loops in \code{R},
//'   or it requires \emph{matrix} transpose.  And both are slow.
//'
//'   The function \code{mult_mat_ref()} uses \code{RcppArmadillo} \code{C++}
//'   code, so when multiplying large \emph{matrix} columns it's several times
//'   faster than vectorized \code{R} code, and it's even much faster compared
//'   to \code{R} when multiplying the \emph{matrix} rows.
//' 
//'   The function \code{mult_mat_ref()} performs loops over the \emph{matrix} rows
//'   and columns using the \emph{Armadillo} operators \code{each_row()} and
//'   \code{each_col()}, instead of performing explicit \code{for()} loops (both
//'   methods are equally fast).
//'   
//' @examples
//' \dontrun{
//' # Create vector and matrix data
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' 
//' # Multiply the matrix rows using R
//' matrixr <- t(vectorv*t(matrixv))
//' # Multiply the matrix rows using C++
//' HighFreq::mult_mat_ref(vectorv, matrixv, byrow=TRUE)
//' all.equal(matrixr, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat_ref(vectorv, matrixv, byrow=TRUE),
//'     Rcode=t(vectorv*t(matrixv)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//'     
//' # Multiply the matrix columns using R
//' matrixr <- vectorv*matrixv
//' # Multiply the matrix columns using C++
//' HighFreq::mult_mat_ref(vectorv, matrixv, byrow=FALSE)
//' all.equal(matrixr, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat_ref(vectorv, matrixv, byrow=FALSE),
//'     Rcode=vectorv*matrixv,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
void mult_mat_ref(arma::vec vector,
                  arma::mat matrix,
                  bool byrow = true) {
  
  arma::uword nelem = vector.n_elem;
  arma::uword nrows = matrix.n_rows;
  arma::uword ncols = matrix.n_cols;
  
  if (byrow && (nelem == ncols)) {
    // Multiply every row of matrix by vector
    matrix.each_row() %= vector.t();
  } else if (!byrow && (nelem == nrows)) {
    // Multiply every column of matrix by vector
    matrix.each_col() %= vector;
  } else {
    // Do nothing
    cout << "Nothing done: Vector length is neither equal to the number of columns nor to the rows of the matrix!" << std::endl;
  }  // end if
  
  // return matrix;
  
  }  // end mult_mat_ref
  
  