////////////////////////////
// Functions to test C++ syntax using the Standard Template Library (STL)
// Test C++ STL code for function pointers and functors.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_rcpp_stl_functor.cpp")

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
double mult_two(double fac_tor, double x) {
  return (fac_tor * x);
}  // end mult_two


// Define a functor for multiplying a number by a fac_tor.
class mult_class {
  
  double fac_tor = 2.0;
  
public:
  // Constructor
  mult_class(double in_put) : fac_tor(in_put) {}
  
  // Overloaded operator - actual function
  double operator()(double x) {return (fac_tor*x);}
  
};  // end mult_class


// Call function pointer.
// [[Rcpp::export]]
std::vector<double> double_vec(std::vector<double> vec_tor) {
  std::vector<double> out_put(vec_tor.size());
  std::transform(vec_tor.begin(), vec_tor.end(), out_put.begin(), double_it);
  return out_put;
}  // end double_vec


// Call functor.
// [[Rcpp::export]]
std::vector<double> mult_vec2(std::vector<double> vec_tor, double fac_tor) {
  
  // Create the instance mult_it of the functor class mult_class
  mult_class mult_it(fac_tor);
  // Define output vector
  std::vector<double> out_put(vec_tor.size());
  
  std::transform(vec_tor.begin(), vec_tor.end(), out_put.begin(), 
                 mult_it
                   // Or pass a lambda function - similar speed
                   // [&fac_tor](double x) {return (fac_tor*x);}
  );
  // Or simply
  // std::transform(vec_tor.begin(), vec_tor.end(), out_put.begin(), mult_class(fac_tor));
  return out_put;
}  // end mult_vec2


// Call lambda function.
// [[Rcpp::export]]
std::vector<double> mult_vec_lambda(std::vector<double> vec_tor, double fac_tor) {
  std::vector<double> out_put(vec_tor.size());
  std::transform(vec_tor.begin(), vec_tor.end(), out_put.begin(), 
                 [&fac_tor](double x) {return (fac_tor*x);}
  );
  return out_put;
}  // end mult_vec_lambda



// Call function pointer.
// std::vector<double> double_vec_funcpt(std::vector<double> vec_tor, double (*func_p)(double)) {
//   std::vector<double> out_put(vec_tor.size());
//   std::transform(vec_tor.begin(), vec_tor.end(), out_put.begin(), func_p);
//   return out_put;
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

typedef SEXP (*ptr_func_2arg)(int, double);

// [[Rcpp::export]]
NumericVector execute_cpp(SEXP func_, double l) {
  ptr_func_2arg func = *XPtr<ptr_func_2arg>(func_);
  return core_processing<ptr_func_2arg>(func, l);
}

typedef double (*ptr_func_1arg)(double, double);
// [[Rcpp::export]]
NumericVector run_cpp(SEXP func_, double fac_tor, NumericVector vec_) {
  ptr_func_1arg func = *XPtr<ptr_func_1arg>(func_);
  
  Rcpp::NumericVector out_put(vec_.size());
  // These two lines produce C++ warnings but they compile fine anyway
  for (int i = 0; i < out_put.size(); i++){
    out_put[i] = func(fac_tor, vec_[i]);
  }
  // Rcpp::NumericVector out_put = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
  return out_put;
}  // end run_cpp

