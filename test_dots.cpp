////////////////////////////
// Functions to test C++ dots syntax
// To pass parameters to worker function
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_dots.cpp")
// Run in R:
// paramv <- list(first=1, second=2, third="blah")
// back_test(paramv)


// Compiler instructions
// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// Included for cpp11
// [[Rcpp::plugins(cpp11)]]


// Below is an example of passing parameters to 
// a worker function through a list.

// calc_weights is the worker function 
// calc_weights is not exported to R
void calc_weights(Rcpp::List paramv) {
  // if (!paramv.inherits("lm")) stop("Input must be a linear model");
  int first = paramv["first"];
  int second = paramv["second"];
  std::string third = paramv["third"];
  
  // NumericVector resid = as<NumericVector>(paramv["first"]);
  // NumericVector fitted = as<NumericVector>(paramv["second"]);
  std::cout << "Parameters passed into calc_weights" << endl;
  std::cout << "first = " << first << endl;
  std::cout << "second = " << second << endl;
  std::cout << "third = " << third << endl;
  
}  // end calc_weights


// back_test passes the list of parameters to calc_weights
// back_test is exported to R
//' @export
// [[Rcpp::export]]
void back_test(Rcpp::List paramv) {
  
  // Perform back_test loop
  for (arma::uword it=1; it < 3; it++) {
    std::cout << "Loop number: " << it << endl;
    // Calculate portfolio weights
    calc_weights(paramv);
  }  // end for
  
}  // end back_test


// Below is code that doesn't work with Rcpp

// Error: undefined reference to `addOne()'
// template <int ...Ns>
// void addOne() {
//   for (int i : {Ns...})
//     std::cout << i+1 << endl;
// }

// Error: Invalid parameter: '...' for Rcpp::export attribute at test_dots.cpp:17 
// void varargs_fn(int num, ...) {
//   va_list valst;
//   va_start(valst, num);
//   int n = 0;
//   
//   for (int i = 0 ; i < num ; i++) {
//     n = va_arg(valst, int);
//     // printf("%d\n", n);
//     std::cout << "Warning: Incorrect typev argument: " << n << endl;
//   }
//   va_end(valst);
// }  // end varargs_fn


// The Dots object is not recognized by compiler
// https://stackoverflow.com/questions/24590946/passing-many-argumentes-by-ellipsis-in-rcpp
// List force_dots(const Dots& dots) {
//   List out(n);
//   for(int i = 0; i < n; i++){
//     out[i] = Rcpp_eval(dots.promise(i), dots.environment(i));
//   }
//   return out ;
// }  // end force_dots()
// 
// List dots_example(NumericVector x, Dots dots) {
//   int n = dots.size();
//   List args = force_dots(dots);
//   return args;
// }
// 
