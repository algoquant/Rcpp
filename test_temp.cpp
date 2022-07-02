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


// The function outer_vec() calculates the outer product of two vectors.
// It accepts pointers to the two vectors and returns a matrix.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::mat outer_vec(const arma::vec& vec1, const arma::vec& vec2) {
  return vec1*vec2.t();
}  // end outer_vec


//' @export
// [[Rcpp::export]]
void run_covmat(arma::mat& covmat, arma::vec& returns, double lambda) {
  
  covmat = (1-lambda)*returns.t()*returns + lambda*covmat;

}  // end run_covmat


//' @export
// [[Rcpp::export]]
arma::mat calc_covmat(arma::mat& returns) {
  
  return returns.t()*returns;
  
}  // end calc_covmat


//' @export
// [[Rcpp::export]]
arma::mat matmeans(arma::mat returns, arma::uword dim) {
  
  return arma::mean(returns, dim);
  
}  // end calc_covmat


//' @export
// [[Rcpp::export]]
double vectorSum(NumericVector x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}  // end vectorSum


// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Sum : public Worker
{
  // source vector
  const RVector<double> input;
  
  // accumulated value
  double value;
  
  // constructors
  Sum(const NumericVector input) : input(input), value(0) {}
  Sum(const Sum& sum, Split) : input(sum.input), value(0) {}
  
  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    value += std::accumulate(input.begin() + begin, input.begin() + end, 0.0);
  }
  
  // join my value with that of another Sum
  void join(const Sum& rhs) {
    value += rhs.value;
  }
};

// [[Rcpp::export]]
double parallelVectorSum(NumericVector x) {
  
  // declare the SumBody instance
  Sum sum(x);
  
  // call parallel_reduce to start the work
  parallelReduce(0, x.length(), sum);
  
  // return the computed sum
  return sum.value;
}


