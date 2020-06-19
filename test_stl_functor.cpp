////////////////////////////
// Test C++ STL code for function pointers and functors.
////////////////////////////

// Compile this C++ file using MinGW:
// C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\test_stl_functor.cpp -o test


// Read more:
// https://www.quantstart.com/articles/Function-Objects-Functors-in-C-Part-1/

#include <iostream>
#include <vector>
#include <algorithm>
// #include <assert.h>

// A function that adds two numbers.
double calc_add(double x1, double x2) {
  return (x1 + x2);
}  // end calc_add

// A function that multiplies two numbers.
double calc_mult(double x1, double x2) {
  return (x1 * x2);
}  // end calc_mult

// The functional do_run() accepts a function pointer and executes it.
double do_run(double (func_tion)(double, double), double x1, double x2) {
  return (func_tion)(x1, x2);
}  // end do_run



// The function is_greater() returns a Boolean depending on the order of the two elements of in_puts.
bool is_greater(const std::vector<double> in_puts, int i1, int i2) {
  return (in_puts[i1] < in_puts[i2]);
}  // end is_greater


// A functor is a C++ class that acts as a function.

// Define a comparison functor as a struct.
struct is_greater_functor {
  // Overloaded operator - actual function
  bool operator()(std::vector<double> in_puts, int i1, int i2) {
    return (in_puts[i1] < in_puts[i2]);
  }
};  // end is_greater_functor


// Functor for multiplying a number by 2
class double_class {
public:
  // Overloaded operator - actual function
  double operator()(double x) { return (2*x);}
};  // end double_class



// Functor for multiplying a number by a fac_tor.
class mult_class {
  
  double fac_tor;
  
public:
  // Constructor
  mult_class(double in_put) : fac_tor(in_put) {}
  
  // Overloaded operator - actual function
  double operator()(double x) {return (fac_tor*x);}
  
};  // end mult_class



int main() {
  
	// int in_put;
	size_t num_el = 7;
  // std::string stri_ng;
  // std::vector<double> vec_tor;
  // std::vector<int> vec_tor = {7, 5, 16, 8, 16, 8};
  
  // double x1 = 5.0;
  // double x2 = 10.0;
  
  // std::cout << "Add: " << do_run(calc_add, x1, x2) << std::endl;
  // std::cout << "Multiply: " << do_run(calc_mult, x1, x2) << std::endl;
  
  double fac_tor;
  std::cout << "Enter multiplier: ";
  std::cin >> fac_tor;
  // Create the instance mult_it of the functor class mult_class
  mult_class mult_it(fac_tor);
  
  double va_r;
  std::cout << "Enter number: ";
  std::cin >> va_r;
  // Call the functor mult_it
  double x = mult_it(va_r);
  // assert(x == (fac_tor*va_r));
  std::cout << "mult_it = " << x << std::endl;
  
  std::vector<double> vec_tor(num_el);
  // Fill vec_tor with a sequence of consecutive integers.
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(vec_tor.begin(), vec_tor.end(), 0);
  // Fill vec_tor with an integers using pointer
  // for (auto &ele_ment: vec_tor) {
  //   ele_ment = 4;
  // }  // end for
  // Old-style loop
  // for (auto ele_ment = 0; ele_ment < num_el; ele_ment++) {
  //   vec_tor[ele_ment] = ele_ment;
  // }  // end for
  std::cout << "These are the vector elements before mult:" << std::endl;
  for (auto ele_ment: vec_tor) {
    std::cout << ele_ment << " ";
  }  // end for
  std::cout << std::endl;
  // std::iota(vec_tor.begin(), vec_tor.begin() + num_el, 0);
  std::vector<double> out_put(num_el);
  std::transform(vec_tor.begin(), vec_tor.end(), out_put.begin(), mult_it);
  std::cout << "These are the vector elements after mult:" << std::endl;
  for (auto ele_ment: out_put) {
    std::cout << ele_ment << " ";
  }  // end for
  std::cout << std::endl;
  
  // Test for is_greater()
  // Create the instance is_greater2 of the functor class is_greater_functor
  is_greater_functor is_greater2;
  std::cout << "is_greater function: " << std::boolalpha << is_greater(out_put, 2, 3) << std::endl;
  std::cout << "is_greater functor: " << std::boolalpha << is_greater2(out_put, 2, 3) << std::endl;


}  // end main

