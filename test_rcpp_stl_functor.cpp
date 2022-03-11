////////////////////////////
// Functions to test C++ syntax using the Standard Template Library (STL)
// Test C++ STL code for function pointers and functors.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_rcpp_stl_functor.cpp")

#include <Rcpp.h>
using namespace Rcpp;
// Use STL
using namespace std;
// [[Rcpp::plugins(cpp11)]]


// Function for multiplying a number by 2.
// [[Rcpp::export]]
double double_it(double x) {
  return (2 * x);
}  // end mult_two


// Function for multiplying a number by 2.
double mult_two(double ratio, double x) {
  return (ratio * x);
}  // end mult_two


// Define a functor for multiplying a number by a ratio.
class mult_class {
  
  double ratio = 2.0;
  
public:
  // Constructor
  mult_class(double input) : ratio(input) {}
  
  // Overloaded operator - actual function
  double operator()(double x) {return (ratio*x);}
  
};  // end mult_class


// Call function pointer.
// [[Rcpp::export]]
std::vector<double> double_vec(std::vector<double> vectorv) {
  std::vector<double> output(vectorv.size());
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), double_it);
  return output;
}  // end double_vec


// Call functor.
// [[Rcpp::export]]
std::vector<double> mult_vec2(std::vector<double> vectorv, double ratio) {
  
  // Create the instance mult_it of the functor class mult_class
  mult_class mult_it(ratio);
  // Define output vector
  std::vector<double> output(vectorv.size());
  
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), 
                 mult_it
                   // Or pass a lambda function - similar speed
                   // [&ratio](double x) {return (ratio*x);}
  );
  // Or simply
  // std::transform(vectorv.begin(), vectorv.end(), output.begin(), mult_class(ratio));
  return output;
}  // end mult_vec2


// Call lambda function.
// [[Rcpp::export]]
std::vector<double> mult_vec_lambda(std::vector<double> vectorv, double ratio) {
  std::vector<double> output(vectorv.size());
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), 
                 [&ratio](double x) {return (ratio*x);}
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

typedef SEXP (*ptr_func2arg)(int, double);

// [[Rcpp::export]]
NumericVector execute_cpp(SEXP func_, double l) {
  ptr_func2arg func = *XPtr<ptr_func2arg>(func_);
  return core_processing<ptr_func2arg>(func, l);
}

typedef double (*ptr_func1arg)(double, double);
// [[Rcpp::export]]
NumericVector run_cpp(SEXP func_, double ratio, NumericVector vec_) {
  ptr_func1arg func = *XPtr<ptr_func1arg>(func_);
  
  Rcpp::NumericVector output(vec_.size());
  // These two lines produce C++ warnings but they compile fine anyway
  for (int i = 0; i < output.size(); i++){
    output[i] = func(ratio, vec_[i]);
  }
  // Rcpp::NumericVector output = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
  return output;
}  // end run_cpp

