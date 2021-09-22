////////////////////////////
// Functions to test C++ syntax using the Standard Template Library (STL)
// Test C++ STL code for function pointers and functors.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_xptr.cpp")

#include <Rcpp.h>
using namespace Rcpp;
// Use STL
using namespace std;
// [[Rcpp::plugins(cpp11)]]

// Below code is inspired by this:
// https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
// https://en.wikipedia.org/wiki/Function_pointer#Alternate_C_and_C++_Syntax
// https://en.wikipedia.org/wiki/Function_pointer#Example_in_C

// Define template of function pointer
typedef double (*func_ptr)(double, double);

// Define function that accepts a function pointer
// [[Rcpp::export]]
NumericVector run_cpp(SEXP func_, double fac_tor, NumericVector vec_) {
  func_ptr func = *XPtr<func_ptr>(func_);
  
  // The lines below produce C++ warnings but they compile fine anyway
  Rcpp::NumericVector out_put(vec_.size());
  
  for (int i = 0; i < out_put.size(); i++){
    out_put[i] = func(fac_tor, vec_[i]);
  }
  // Rcpp::NumericVector out_put = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
  return out_put;
}  // end run_cpp



//////////////////////////////

// Below is copied from:
// https://gallery.rcpp.org/articles/passing-cpp-function-pointers-rcppxptrutils/

template <typename T>
NumericVector core_processing(T func, double l) {
  double accum = 0;
  for (int i = 0; i < 10; i++)
    accum += sum(as<NumericVector>(func(3, l)));
  return NumericVector(1, accum);
}

// [[Rcpp::export]]
NumericVector execute_r(Function func, double l) {
  return core_processing<Function>(func, l);
}

typedef SEXP (*func_ptr_sexp)(int, double);

// [[Rcpp::export]]
NumericVector execute_cpp(SEXP func_, double l) {
  func_ptr_sexp func = *XPtr<func_ptr_sexp>(func_);
  return core_processing<func_ptr_sexp>(func, l);
}

