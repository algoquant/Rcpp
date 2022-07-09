////////////////////////////
// Test C++ Rcpp code for Rcpp::XPtr function pointers.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_xptr.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
// Use STL
using namespace std;

// Below code is inspired by these:
// https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
// https://stackoverflow.com/questions/14428687/rcpparmadillo-pass-user-defined-function/14428758#14428758
// https://en.wikipedia.org/wiki/Function_pointer#Alternate_C_and_C++_Syntax
// https://en.wikipedia.org/wiki/Function_pointer#Example_in_C


// Define a few regular functions

// [[Rcpp::export]]
double mult2(double factorv, double x) {
  return factorv*x;
}  // end mult2

// Double the vector
arma::vec fun1(const arma::vec& tseries) {
  arma::vec y = tseries + tseries;
  return (y);
}  // end fun1

// Multiply the vector by 10
arma::vec fun2(const arma::vec& tseries) {
  arma::vec y = 10*tseries;
  return (y);
}  // end fun2

// [[Rcpp::export]]
arma::mat run_mean(const arma::mat& tseries, double lambda) {
  
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


// [[Rcpp::export]]
arma::mat run_max(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat maxv = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  maxv.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the max using the decay factor
    maxv.row(it) = lambda*arma::max(tseries.row(it), maxv.row(it-1)) + lambda1*tseries.row(it);
    // Alternative formula for the same:
    // maxv.row(it) = arma::max(tseries.row(it), lambda*maxv.row(it-1) + lambda1*tseries.row(it));
  }  // end for
  
  return maxv;
  
}  // end run_max


// Define type for pointer to function of a vector
// typedef double (*funptr_vec)(double, double);
typedef arma::vec (*funptr_vec)(const arma::vec& tseries);
// Define type for pointer to function of a matrix and a numeric
typedef arma::mat (*funptr_mat_num)(const arma::mat& tseries, double lambda);


// Return a pointer to a function of a vector based on a string
// [[Rcpp::export]]
Rcpp::XPtr<funptr_vec> getfunptr_vec(std::string funname) {
  if (funname == "fun1")
    return (Rcpp::XPtr<funptr_vec>(new funptr_vec(&fun1)));
  else if (funname == "fun2")
    return (Rcpp::XPtr<funptr_vec>(new funptr_vec(&fun2)));
  else
    return Rcpp::XPtr<funptr_vec>(R_NilValue); // runtime error as NULL no Rcpp::XPtr
}  // end getfunptr_vec

// Run a function of a vector based on a string
// [[Rcpp::export]]
arma::vec runfunvec(std::string funname, const arma::vec& tseries) {
  Rcpp::XPtr<funptr_vec> xpfun = getfunptr_vec(funname);
  funptr_vec fun = *xpfun;
  arma::vec y = fun(tseries);
  return (y);
}  // end runfunvec



// Return a pointer to a function of a matrix based on a string
// [[Rcpp::export]]
Rcpp::XPtr<funptr_mat_num> getfunptr_mat_num(std::string funname) {
  if (funname == "run_mean")
    return (Rcpp::XPtr<funptr_mat_num>(new funptr_mat_num(&run_mean)));
  else if (funname == "run_max")
    return (Rcpp::XPtr<funptr_mat_num>(new funptr_mat_num(&run_max)));
  else
    return Rcpp::XPtr<funptr_mat_num>(R_NilValue); // runtime error as NULL no Rcpp::XPtr
}  // end getfunptr_mat_num

// Run a function of a matrix based on a string
// [[Rcpp::export]]
arma::mat runfunmat(std::string funname, arma::mat& tseries, double lambda) {
  Rcpp::XPtr<funptr_mat_num> xpfun = getfunptr_mat_num(funname);
  funptr_mat_num fun = *xpfun;
  arma::mat outp = fun(tseries, lambda);
  return (outp);
}  // end runfunmat


// Define function that accepts a function pointer
// Rcpp::NumericVector run_cpp(SEXP func_, double ratio, Rcpp::NumericVector vec_) {
//   funptr_vec func = *Rcpp::XPtr<funptr_vec>(func_);
//   
//   // The lines below produce C++ warnings but they compile fine anyway
//   Rcpp::NumericVector output(vec_.size());
//   
//   for (int i = 0; i < output.size(); i++){
//     output[i] = func(ratio, vec_[i]);
//   }
//   // Rcpp::NumericVector output = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
//   return output;
// }  // end run_cpp



//////////////////////////////

// Below is copied from:
// https://gallery.rcpp.org/articles/passing-cpp-function-pointers-rcppxptrutils/

// Define variable template
template <typename T>
// Define Rcpp functional
Rcpp::NumericVector core_processing(T func, double l) {
  double accum = 0;
  for (int i = 0; i < 10; i++)
    accum += sum(as<Rcpp::NumericVector>(func(3, l)));
  return Rcpp::NumericVector(1, accum);
}

// Use Rcpp::Function object to pass an R function into Rcpp
// https://teuder.github.io/rcpp4everyone_en/230_R_function.html
// [[Rcpp::export]]
Rcpp::NumericVector execute_r(Rcpp::Function func, double l) {
  return core_processing<Rcpp::Function>(func, l);
}

// Define R pointer to function of integer and numeric arguments
typedef SEXP (*funptr_vec_intdouble)(int, double);

// Functional accepts an R SEXP pointer to a function
// [[Rcpp::export]]
Rcpp::NumericVector execute_cpp(SEXP func_, double l) {
  funptr_vec_intdouble func = *XPtr<funptr_vec_intdouble>(func_);
  return core_processing<funptr_vec_intdouble>(func, l);
}

// ptr2fun <- cppXPtr(double_it)

// Define type for pointer to function of two numeric arguments
typedef double (*ptr_func2arg)(double, double);
// Define an Rcpp functional to run an Rcpp function of two numeric arguments in a loop over a vector
// The syntax below doesn't require if-else statements because it uses a generic function type
// [[Rcpp::export]]
Rcpp::NumericVector run_cpp(SEXP func_, double factorv, Rcpp::NumericVector vec_) {
  // Coerce function pointer to type ptr_func2arg
  ptr_func2arg func = *XPtr<ptr_func2arg>(func_);
  
  Rcpp::NumericVector output(vec_.size());
  // These two lines produce C++ warnings but they compile fine anyway
  for (int i = 0; i < output.size(); i++){
    output[i] = func(factorv, vec_[i]);
  }
  // Rcpp::NumericVector output = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
  return output;
}  // end run_cpp

