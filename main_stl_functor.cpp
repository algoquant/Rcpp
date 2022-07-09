////////////////////////////
// Test C++ STL code for functors, functionals, and function pointers.
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  cd /Users/jerzy/Develop/Rcpp/
//  g++ -std=c++11 main_stl_functor.cpp util_stl_fun.cpp -o test
// Compile this file in R by running this command:
//  Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/main_stl_functor.cpp")

// Read more:
// https://www.cprogramming.com/tutorial/functors-function-objects-in-c++.html
// https://www.quantstart.com/articles/Function-Objects-Functors-in-C-Part-1/

using namespace std;
#include <vector>
#include "util_stl_fun.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
// #include <assert.h>

// The functional run_fun1() accepts a pointer to a function of two numeric variables, and executes it.
double run_fun1(double fun(double, double), double x1, double x2) {
  return fun(x1, x2);
}  // end run_fun1


// The functional run_fun2() accepts a pointer to a function of two variables (numeric and vector), and executes it.
std::vector<double> run_fun2(std::vector<double> fun(double, std::vector<double>), double factorv, std::vector<double> vectorv) {
  return fun(factorv, vectorv);
}  // end run_fun2


// The function is_greater() returns a Boolean depending on the order of the two elements of the inputs vector.
bool is_greater(const std::vector<double> inputs, int i1, int i2) {
  return (inputs[i1] < inputs[i2]);
}  // end is_greater


// A functor is a C++ class that acts as a function.

// Define a comparison functor.
class is_greater_functor {
public:
  // Overloaded operator - actual function
  bool operator()(std::vector<double> inputs, int i1, int i2) {
    return (inputs[i1] < inputs[i2]);
  }
};  // end is_greater_functor


// Functor for multiplying a number by 2
class double_class {
public:
  // Overloaded operator - actual function
  double operator()(double x) { return (2*x);}
};  // end double_class



// Functor class for multiplying a number by a factor.
class mult_class {
  
private:
  double factorv;
  
public:
  // Constructor
  mult_class(double input) : factorv(input) {}
  
  // Overloaded operator - the actual function
  double operator()(double x) {return (factorv*x);}
  
};  // end mult_class


// Uncomment the code below that you want to run
int main() {
  
	// int input;
  // std::vector<double> datav;
  // std::vector<int> datav = {7, 5, 16, 8, 16, 8};
  
  std::cout << "Execute functions using a functional and function pointers." << std::endl;
  std::cout << "Enter two numbers." << std::endl;
  double x1 = 5.0;
  double x2 = 10.0;
  std::cout << "Enter x1: ";
  std::cin >> x1;
  std::cout << "Enter x2: ";
  std::cin >> x2;
  
  // The functional run_fun1() accepts a pointer to a function of two numeric variables, and executes it.
  std::cout << "Add (x1 + x2): " << run_fun1(add2, x1, x2) << std::endl;
  std::cout << "Multiply (x1 * x2): " << run_fun1(mult2, x1, x2) << std::endl;
  std::cout << std::endl;
  
  std::cout << "Multiply a number using a functor." << std::endl;
  // Use the functor class mult_class
  double factorv;
  std::cout << "Enter multiplier factor: ";
  std::cin >> factorv;
  // Create the instance mult_it of the functor class mult_class
  mult_class mult_it(factorv);
  
  double numv;
  std::cout << "Enter a number: ";
  std::cin >> numv;
  // Call the functor mult_it
  double x = mult_it(numv);
  // assert(x == (factorv*numv));
  std::cout << "multiplier * number = " << x << std::endl;
  std::cout << std::endl;
  
  std::cout << "Multiply a vector using a functor and transform." << std::endl;
  std::size_t ndata;
  std::cout << "Enter the number of vector elements: ";
  std::cin >> ndata;
  // Use transform() and the functor class mult_class
  std::vector<double> datav(ndata);
  // Fill datav with a sequence of consecutive integers.
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(datav.begin(), datav.end(), 0);
  
  // Fill datav with an integers using pointer
  // for (auto &it: datav) {
  //   it = 4;
  // }  // end for
  // Old-style loop
  // for (auto it = 0; it < ndata; it++) {
  //   datav[it] = it;
  // }  // end for
  
  std::cout << "These are the vector elements before multiplication:" << std::endl;
  for (auto it: datav) {
    std::cout << it << " ";
  }  // end for
  std::cout << std::endl;
  // std::iota(datav.begin(), datav.begin() + ndata, 0);
  std::vector<double> output(ndata);
  std::transform(datav.begin(), datav.end(), output.begin(), mult_it);
  std::cout << "These are the vector elements after first multiplication:" << std::endl;
  for (auto it: output) {
    std::cout << it << " ";
  }  // end for
  std::cout << std::endl;
  
  
  output = run_fun2(mult_vec_lambda, factorv, datav);
  std::cout << "These are the vector elements after second multiplication:" << std::endl;
  for (auto it: output) {
    std::cout << it << " ";
  }  // end for
  std::cout << std::endl;
  
  // Decay factor
  double lambda = 0.9;
  std::fill(datav.begin(), datav.end(), 1.0);
  output = run_fun2(run_mean, lambda, datav);
  std::cout << "These are the running weighted means of the vector:" << std::endl;
  for (auto it: output) {
    std::cout << it << " ";
  }  // end for
  std::cout << std::endl;
  
  // Test for is_greater()
  // Create the instance is_greater2 of the functor class is_greater_functor
  is_greater_functor is_greater2;
  std::cout << "is_greater function: " << std::boolalpha << is_greater(output, 2, 3) << std::endl;
  std::cout << "is_greater functor: " << std::boolalpha << is_greater2(output, 2, 3) << std::endl;


}  // end main

