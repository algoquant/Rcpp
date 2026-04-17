#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/mult_mat.cpp")

// The function mult_mat() calculates the positions of 
// a Bollinger strategy, using the z-scores.
// It doubles down and adds to the position when the 
// z-score increases in the same direction.
// It holds its position as long as the z-score doesn't 
// change sign.
// It flips the position when the z-score changes sign.
// The strategy is path dependent and requires performing 
// a loop over time.

// benchmark in R
// weightv <- invred %*% colmeanv
// weightvv <- mult_mat(colmeanv, invred)
// weightvvv <- mult_matin(colmeanv, invred)
// all.equal(weightvvv, weightvv)
// library(microbenchmark)
// summary(microbenchmark(
//     innerr={invred %*% colmeanv},
//     innerarma=mult_mat(colmeanv, invred),
//     innercpp=mult_matin(colmeanv, invred),
//     times=10))[, c(1, 4, 5)]


using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat mult_mat(arma::vec vectorv,
                   arma::mat matrixv) {
  
  // Print the dimensions of the input vector and matrix
  // Rcpp::Rcout << "Dimensions of vectorv: " << vectorv.n_rows
  //             << " x " << vectorv.n_cols << std::endl;
  // Rcpp::Rcout << "Dimensions of matrixv: " << matrixv.n_rows
  //             << " x " << matrixv.n_cols << std::endl;
  
  // mat A(2, 3, fill::randu);
  // vec B(2,    fill::randu);
  // 
  // vec C = arma::affmul(A,B);
  // return C;
  // return arma::affmul(matrixv, vectorv);
  return (matrixv*vectorv);
  
}  // end mult_mat


// [[Rcpp::export]]
arma::mat mult_matin(arma::vec vectorv,
                   arma::mat matrixv) {
  
  // Multiply the matrix by the vector using loops
  // Get the number of rows and columns in the matrix
  int nrows = matrixv.n_rows;
  int ncols = matrixv.n_cols;
  // Create the output vector to store the output
  arma::vec outv(nrows, fill::zeros);
  // Multiply every row of matrixv by vectorv
  matrixv = matrixv.each_row() % vectorv.t();
  // Sum the rows of the resulting matrix to get the output vector
  outv = sum(matrixv, 1);
  
  // This is faster but it conflicts with mclapply because it uses OpenMP and causes a crash when used in parallel
  // return (matrixv*vectorv);
  
  // The standard C++ loops below are many times slower than the vectorized operations above
  // Loop through each row of the matrix
  // for (int i = 0; i < nrows; i++) {
  //   // Loop through each column of the matrix
  //   for (int j = 0; j < ncols; j++) {
  //     // Multiply the matrix element by the corresponding vector element
  //     outv(i) += matrixv(i, j) * vectorv(j);
  //   }
  // }

  return (outv);
  
}  // end mult_matin


