////////////////////////////
// Test C++ Rcpp code for function pointers (not Rcpp::XPtr).
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_funptr.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
// Use STL
using namespace std;

// Below code is inspired by these:
// https://en.wikipedia.org/wiki/Function_pointer#Alternate_C_and_C++_Syntax
// https://en.wikipedia.org/wiki/Function_pointer#Example_in_C

// Define a few functions

// Add two numbers.
// [[Rcpp::export]]
double add2(double x1, double x2) {
  return (x1 + x2);
}  // end add2

// Multiply two numbers.
// [[Rcpp::export]]
double mult2(double x1, double x2) {
  return (x1*x2);
}  // end mult2

// Subtract 2 from vector
// [[Rcpp::export]]
arma::vec fun1(const arma::vec& tseries) {
  arma::vec y = tseries - 2;
  return y;
}  // end fun1

// Multiply the vector by 2
// [[Rcpp::export]]
arma::vec fun2(const arma::vec& tseries) {
  arma::vec y = 2*tseries;
  return y;
}  // end fun2


//////////////////
// Run a function of two double numbers based on a string.
// This method uses a functional without using a function pointer.

// The functional run2d() accepts a pointer to a function of two numeric variables, and executes it.
double run2d(double fun(double, double), double x1, double x2) {
  return fun(x1, x2);
}  // end run2d

// Run a function of two double numbers based on a string
// [[Rcpp::export]]
double runfun2(std::string funname, double x1, double x2) {
  if (funname == "add2")
    return run2d(add2, x1, x2);
  else if (funname == "mult2")
    return run2d(mult2, x1, x2);
  else
    throw std::invalid_argument("No such function!");
}  // end runfun2



//////////////////
// Run a function of two double numbers based on a string.
// This method uses a function pointer instead of a functional.

// Define type for pointer to function of two double numbers
typedef double (*funptr)(double, double);


// Return a pointer to a function of two double numbers based on a string
funptr getfunptr(std::string funname) {
  if (funname == "add2")
    return (&add2);
  else if (funname == "mult2")
    return (&mult2);
  else
    throw std::invalid_argument("No such function!");
}  // end getfunptr

// Run a function of two double numbers based on a string
// [[Rcpp::export]]
double runfun(std::string funname, double x1, double x2) {
  funptr funptrv = getfunptr(funname);
  // funptr fun = *funptrv;
  double y = funptrv(x1, x2);
  return y;
}  // end runfun


//////////////////
// Run a function of a vector based on a string.
// This method uses a function pointer.

// Define type for pointer to function of a vector
typedef arma::vec (*funptrvec)(const arma::vec& tseries);
// Define type for pointer to function of a matrix and a numeric
// typedef arma::mat (*funptr_mat_num)(const arma::mat& tseries, double lambda);

// Return a pointer to a function of a vector based on a string
funptrvec getfunptrvec(std::string funname) {
  if (funname == "fun1")
    return (&fun1);
  else if (funname == "fun2")
    return (&fun2);
  else
    throw std::invalid_argument("No such function!");
}  // end getfunptrvec

// Run a function of a vector based on a string
// [[Rcpp::export]]
arma::vec runfunvec(std::string funname, arma::vec& tseries) {
  funptrvec funptrv = getfunptrvec(funname);
  // funptrvec fun = *funptrv;
  arma::vec y = funptrv(tseries);
  return y;
}  // end runfunvec

