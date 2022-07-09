////////////////////////////
// Functions to test C++ syntax using the Standard Template Library (STL)
// Test C++ STL code for function pointers and functors.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_rcpp_stl_functor.cpp")

#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]

// Function for multiplying a number by 2.
// [[Rcpp::export]]
double double_it(double x) {
  return (2 * x);
}  // end mult_two


// Function for multiplying a number by 2.
// [[Rcpp::export]]
double mult_two(double factorv, double x) {
  return (factorv * x);
}  // end mult_two


// Define a functor for multiplying a number by a factorv.
class mult_class {
  
private:
  double factorv;
  
public:
  // Constructor
  mult_class(double input) : factorv(input) {}
  
  // Overloaded operator - the actual function
  double operator()(double x) {return (factorv*x);}
  
};  // end mult_class


// Call function inside transform.
// [[Rcpp::export]]
std::vector<double> double_vec(std::vector<double> vectorv) {
  std::vector<double> output(vectorv.size());
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), double_it);
  return output;
}  // end double_vec


// Call functor.
// [[Rcpp::export]]
std::vector<double> mult_vec(std::vector<double> vectorv, double factorv) {
  
  // Create the instance mult_it of the functor class mult_class
  mult_class mult_it(factorv);
  // Define output vector
  std::vector<double> output(vectorv.size());
  
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), 
                 // Call functor
                 mult_it
                   // Or pass a lambda function - similar speed
                   // [&factorv](double x) {return (factorv*x);}
  );
  // Or simply
  // std::transform(vectorv.begin(), vectorv.end(), output.begin(), mult_class(factorv));
  return output;
}  // end mult_vec


// Call lambda function.
// [[Rcpp::export]]
std::vector<double> mult_vec_lambda(std::vector<double> vectorv, double factorv) {
  std::vector<double> output(vectorv.size());
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), 
                 [&factorv](double x) {return (factorv*x);}
  );
  return output;
}  // end mult_vec_lambda



// Call function pointer.
// std::vector<double> double_vec_funcpt(std::vector<double> vectorv, double (*func_p)(double)) {
//   std::vector<double> output(vectorv.size());
//   std::transform(vectorv.begin(), vectorv.end(), output.begin(), func_p);
//   return output;
// }  // end double_vec_funcpt

//////////////////////////////

// Some functions below
// [[Rcpp::export]]
arma::mat run_mean(const arma::mat tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  means.row(0) = tseries.row(0);
  // Calculate means without weights
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the means using the decay factor
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
  }  // end for
  
  return means;
  
}  // end run_mean

// The functional run_fun2() accepts a pointer to a function of two variables (numeric and vector), and executes it.



