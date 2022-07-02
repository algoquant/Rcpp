////////////////////////////
// Matrix functions using Armadillo and the Standard Template Library (STL)
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/HighFreq/src/matrix_fun.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// Use STL
// using namespace std;
// #include <functional>
// #include <string>
// #include <map>
// #include <iostream>
// [[Rcpp::plugins(cpp11)]]


////////////////////////////////////////////////////////////
// Functions for matrix algebra
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Multiply the rows or columns of a \emph{matrix} times a \emph{vector},
//' element-wise.
//' 
//' @param \code{vector} A \emph{numeric} \emph{vector}.
//' 
//' @param \code{matrix} A \emph{numeric} \emph{matrix}.
//' 
//' @param \code{byrow} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the rows of \code{matrix} by \code{vector}, otherwise multiply the columns
//'   (the default is \code{byrow = TRUE}.)
//' 
//' @return A \emph{matrix} equal to the product of the arguments \code{matrix}
//'   times \code{vector}, with the same dimensions as the argument
//'   \code{matrix}.
//' 
//' @details
//'   The function \code{mult_mat()} multiplies the rows or columns of a
//'   \emph{matrix} times a \emph{vector}, element-wise.
//'
//'   If \code{byrow = TRUE} (the default), then function \code{mult_mat()}
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
//'   The function \code{mult_mat()} uses \code{RcppArmadillo} \code{C++}
//'   code, so when multiplying large \emph{matrix} columns it's several times
//'   faster than vectorized \code{R} code, and it's even much faster compared
//'   to \code{R} when multiplying the \emph{matrix} rows.
//' 
//'   The function \code{mult_mat()} performs loops over the \emph{matrix} rows
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
//' matrixp <- HighFreq::mult_mat(vectorv, matrixv, byrow=TRUE)
//' all.equal(matrixr, matrixp)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat(vectorv, matrixv, byrow=TRUE),
//'     Rcode=t(vectorv*t(matrixv)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//'     
//' # Multiply the matrix columns using R
//' matrixr <- vectorv*matrixv
//' # Multiply the matrix columns using C++
//' matrixp <- HighFreq::mult_mat(vectorv, matrixv, byrow=FALSE)
//' all.equal(matrixr, matrixp)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat(vectorv, matrixv, byrow=FALSE),
//'     Rcode=vectorv*matrixv,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat mult_mat(arma::vec vector,
                   arma::mat matrix,
                   bool byrow = true) {
  
  arma::uword ndata = vector.n_elem;
  arma::uword nrows = matrix.n_rows;
  arma::uword ncols = matrix.n_cols;
  
  if (byrow && (ndata == ncols)) {
    // Multiply every row of matrix by vector
    matrix = matrix.each_row() % vector.t();
  } else if (!byrow && (ndata == nrows)) {
    // Multiply every column of matrix by vector
    matrix = matrix.each_col() % vector;
  } else {
    // Do nothing
    cout << "Nothing done: Vector length is neither equal to the number of columns nor to the rows of the matrix!" << std::endl;
  }  // end if
  
  return matrix;
  
}  // end mult_mat



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
  
  arma::uword ndata = vector.n_elem;
  arma::uword nrows = matrix.n_rows;
  arma::uword ncols = matrix.n_cols;
  
  if (byrow && (ndata == ncols)) {
    // Multiply every row of matrix by vector
    matrix.each_row() %= vector.t();
  } else if (!byrow && (ndata == nrows)) {
    // Multiply every column of matrix by vector
    matrix.each_col() %= vector;
  } else {
    // Do nothing
    cout << "Nothing done: Vector length is neither equal to the number of columns nor to the rows of the matrix!" << std::endl;
  }  // end if
  
  // return matrix;
  
}  // end mult_mat_ref



////////////////////////////////////////////////////////////
//' Calculate the eigen decomposition of the \emph{covariance matrix} of returns
//' data using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of returns
//'   data.
//'
//' @return A list with two elements: a \emph{vector} of eigenvalues 
//'   (named "values"), and a \emph{matrix} of eigenvectors (named
//'   "vectors").
//'
//' @details
//'   The function \code{calc_eigen()} first calculates the \emph{covariance
//'   matrix} of \code{tseries}, and then calculates the eigen
//'   decomposition of the \emph{covariance matrix}.
//'
//' @examples
//' \dontrun{
//' # Create matrix of random data
//' datav <- matrix(rnorm(5e6), nc=5)
//' # Calculate eigen decomposition
//' eigend <- HighFreq::calc_eigen(scale(datav, scale=FALSE))
//' # Calculate PCA
//' pcad <- prcomp(datav)
//' # Compare PCA with eigen decomposition
//' all.equal(pcad$sdev^2, drop(eigend$values))
//' all.equal(abs(unname(pcad$rotation)), abs(eigend$vectors))
//' # Compare the speed of Rcpp with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_eigen(datav),
//'   Rcode=prcomp(datav),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_eigen(const arma::mat& tseries) {
  
  arma::mat eigenvec;
  arma::vec eigenval;
  arma::eig_sym(eigenval, eigenvec, arma::cov(tseries));
  
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  return Rcpp::List::create(Named("values") = arma::flipud(eigenval),
                            Named("vectors") = arma::fliplr(eigenvec));
  
}  // end calc_eigen




////////////////////////////////////////////////////////////
//' Calculate the regularized inverse of a \emph{matrix} of data using Singular
//' Value Decomposition (\emph{SVD}).
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of data.
//' 
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   matrix \code{tseries} (the default is \code{0.01}).
//'   
//' @param \code{dimax} An \emph{integer} equal to the number of singular
//'   values used for calculating the regularized inverse of the matrix
//'   \code{tseries} (the default is \code{dimax = 0} - equivalent to
//'   \code{dimax} equal to the number of columns of \code{tseries}).
//'
//' @return A \emph{matrix} equal to the regularized inverse of the matrix
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_inv()} calculates the regularized inverse of the
//'   matrix \code{tseries} using Singular Value Decomposition (\emph{SVD}).
//'   
//'   The function \code{calc_inv()} first performs Singular Value Decomposition
//'   (\emph{SVD}) of the matrix \code{tseries}.  
//'   The \emph{SVD} of a matrix \eqn{\strong{C}} is defined as the
//'   factorization:
//'   \deqn{
//'     \strong{C} = \strong{U}  \, \Sigma  \, \strong{V}^T
//'   }
//'   Where \eqn{\strong{U}} and \eqn{\strong{V}} are the left and right
//'   \emph{singular matrices}, and \eqn{\Sigma} is a diagonal matrix of
//'   \emph{singular values} \eqn{\Sigma = \{\sigma_i\}}.
//'   
//'   The inverse \eqn{\strong{C}^{-1}} of the matrix \eqn{\strong{C}} can be
//'   calculated from the \emph{SVD} matrices as:
//'   \deqn{
//'     \strong{C}^{-1} = \strong{V} \, \Sigma^{-1} \, \strong{U}^T
//'   }
//'   
//'   The \emph{regularized inverse} of the matrix \eqn{\strong{C}} is given
//'   by:
//'   \deqn{
//'     \strong{C}^{-1} = \strong{V}_n \, \Sigma_n^{-1} \, \strong{U}_n^T
//'   }
//'   Where \eqn{\strong{U}_n}, \eqn{\strong{V}_n} and \eqn{\Sigma_n} are the
//'   \emph{SVD} matrices with the rows and columns corresponding to zero
//'   \emph{singular values} removed.
//'   
//'   The function \code{calc_inv()} performs regularization by discarding the
//'   smallest singular values \eqn{\sigma_i} that are less than the threshold
//'   level \code{eigen_thresh} times the sum of all the singular values:
//'   \deqn{\sigma_i < eigen\_thresh \cdot (\sum{\sigma_i})}
//'   
//'   It then discards additional singular values so that only the largest
//'   \code{dimax} singular values remain.  
//'   It calculates the regularized inverse from the \emph{SVD} matrices using
//'   only the largest singular values up to \code{dimax}.  For example, if
//'   \code{dimax = 3} then it only uses the \code{3} largest singular
//'   values. This has the effect of dimension reduction.
//'   
//'   If the matrix \code{tseries} has a large number of small singular values,
//'   then the number of remaining singular values may be less than
//'   \code{dimax}.
//'   
//' @examples
//' \dontrun{
//' # Calculate ETF returns
//' returns <- na.omit(rutils::etfenv$returns)
//' # Calculate covariance matrix
//' covmat <- cov(returns)
//' # Calculate regularized inverse using RcppArmadillo
//' invmat <- HighFreq::calc_inv(covmat, dimax=3)
//' # Calculate regularized inverse from SVD in R
//' svdec <- svd(covmat)
//' dimax <- 1:3
//' invsvd <- svdec$v[, dimax] %*% (t(svdec$u[, dimax]) / svdec$d[dimax])
//' # Compare RcppArmadillo with R
//' all.equal(invmat, invsvd)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword dimax = 0) {
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD
  arma::svd(svdu, svdval, svdv, tseries);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > eigen_thresh*arma::sum(svdval));
  
  if (dimax == 0) {
    // Set dimax
    dimax = svdnum - 1;
  } else {
    // Adjust dimax
    dimax = std::min(dimax - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, dimax);
  svdu = svdu.cols(0, dimax);
  svdv = svdv.cols(0, dimax);
  
  // Calculate the regularized inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_inv



////////////////////////////////////////////////////////////
//' Scale (standardize) the columns of a \emph{matrix} of data using
//' \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of data.
//' 
//' @param \code{use_median} A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}).
//'   If \code{use_median = FALSE} then the centrality is calculated as the
//'   \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation} (the default is \code{FALSE})
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{tseries}.
//'
//' @details
//'   The function \code{calc_scaled()} scales (standardizes) the columns of the
//'   \code{tseries} argument using \code{RcppArmadillo}.
//'
//'   If the argument \code{use_median} is \code{FALSE} (the default), then it
//'   performs the same calculation as the standard \code{R} function
//'   \code{scale()}, and it calculates the centrality (central tendency) as the
//'   \emph{mean} and the dispersion as the \emph{standard deviation}.
//'
//'   If the argument \code{use_median} is \code{TRUE}, then it calculates the
//'   centrality as the \emph{median} and the dispersion as the \emph{median
//'   absolute deviation} (\emph{MAD}).
//'
//'   If the number of rows of \code{tseries} is less than \code{3} then it
//'   returns \code{tseries} unscaled.
//'   
//'   The function \code{calc_scaled()} uses \code{RcppArmadillo} \code{C++}
//'   code and is about \emph{5} times faster than function \code{scale()}, for
//'   a \emph{matrix} with \emph{1,000} rows and \emph{20} columns.
//'   
//' @examples
//' \dontrun{
//' # Create a matrix of random data
//' returns <- matrix(rnorm(20000), nc=20)
//' scaled <- calc_scaled(tseries=returns, use_median=FALSE)
//' scaled2 <- scale(returns)
//' all.equal(scaled, scaled2, check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=calc_scaled(tseries=returns, use_median=FALSE),
//'   Rcode=scale(returns),
//'   times=100))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_scaled(const arma::mat& tseries, bool use_median = false) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  arma::mat scaledmat(nrows, ncols, fill::zeros);
  arma::vec scaled(nrows, fill::zeros);
  double center;
  
  // Return zeros if not enough data
  if (nrows < 3) {
    return tseries;
  }  // end if
  
  // Perform a loop over the columns
  for (arma::uword it=0; it < ncols; it++) {
    if (use_median) {
      center = arma::median(tseries.col(it));
      scaled = (tseries.col(it) - center);
      scaled = scaled/arma::median(arma::abs(scaled));
      scaledmat.col(it) = scaled;
    } else {
      center = arma::mean(tseries.col(it));
      scaled = (tseries.col(it) - center);
      scaled = scaled/arma::stddev(scaled);
      scaledmat.col(it) = scaled;
    }  // end if
  }  // end for
  
  return scaledmat;
  
}  // end calc_scaled




