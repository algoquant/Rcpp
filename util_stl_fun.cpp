////////////////////////////
// Utility functions using the Standard Template Library (STL)
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/util_fun.cpp")

// Use STL
using namespace std;
#include <vector>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>


////////////////////////////////////////////////////////////
// Functions miscellaneous
////////////////////////////////////////////////////////////


// Add two numbers.
double add2(double x1, double x2) {
  return (x1 + x2);
}  // end add2

// Multiply two numbers.
double mult2(double x1, double x2) {
  return (x1 * x2);
}  // end mult2

// Multiply a vector by a number.
std::vector<double> mult_vec_lambda(double factorv, std::vector<double> vectorv) {
  std::vector<double> output(vectorv.size());
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), 
                 // Pass a lambda function
                 [&factorv](double x) {return (factorv*x);}
  );
  return output;
}  // end mult_vec_lambda


// Calculate the running weighted means of streaming time series data.
std::vector<double> run_mean(double lambda, std::vector<double> vectorv) {
  
  std::size_t nrows = vectorv.size();
  std::vector<double> means(nrows);
  double lambda1 = 1-lambda;
  
  means[0] = vectorv[0];
  // Calculate the running means using the decay factor lambda
  for (std::size_t it = 1; it < nrows; it++) {
    means[it] = lambda1*vectorv[it] + lambda*means[it-1];
  }  // end for
  
  return means;
  
}  // end run_mean

std::vector<double> run_int(double lambda, std::vector<int> vectorv) {
  
  std::size_t nrows = vectorv.size();
  std::vector<double> means(nrows);
  double lambda1 = 1-lambda;
  
  means[0] = vectorv[0];
  // Calculate the running means using the decay factor lambda
  for (std::size_t it = 1; it < nrows; it++) {
    means[it] = lambda1*vectorv[it] + lambda*means[it-1];
  }  // end for
  
  return means;
  
}  // end run_int





