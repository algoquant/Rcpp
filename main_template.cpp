////////////////////////////
// Test C++ templates
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  cd /Users/jerzy/Develop/Rcpp/
//  g++ -std=c++11 main_template.cpp util_stl_fun.cpp -o test
// Run it in Mac terminal: ./test
// 
// Compile this C++ file under Windows using MinGW:
//  C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\main_template.cpp fun_test.cpp -o test

// Read more:
// https://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file
// https://www.codeproject.com/Articles/48575/How-to-Define-a-Template-Class-in-a-h-File-and-Imp


using namespace std;
#include <vector>
#include "util_stl_fun.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

int main() {
  
  std::cout << "Calculate the maximum of two integers using a template function." << std::endl;
  std::cout << "Enter two numbers." << std::endl;
  int int1;
  int int2;
  std::cout << "Enter int1: ";
  std::cin >> int1;
  std::cout << "Enter int2: ";
  std::cin >> int2;
  // Call template function for integer arguments
  std::cout << calc_max(int1, int2) << std::endl;
  
  std::cout << "Calculate the maximum of two float numbers using a template function." << std::endl;
  std::cout << "Enter two float numbers." << std::endl;
  double num1;
  double num2;
  std::cout << "Enter num1: ";
  std::cin >> num1;
  std::cout << "Enter num2: ";
  std::cin >> num2;
  // Call template function for float arguments
  std::cout << calc_max(num1, num2) << std::endl;
  
  
  std::cout << "Multiply a vector using a template functional." << std::endl;
  std::cout << "Enter the number of vector elements: ";
  std::size_t ndata;
  std::cin >> ndata;
  // Decay factor
  double lambda = 0.9;
  std::cout << "Enter the decay factor: ";
  std::cin >> lambda;
  std::vector<int> intv(ndata);
  std::fill(intv.begin(), intv.end(), 2);
  std::vector<double> output(ndata);
  output = run_fun3(run_int, lambda, intv);
  // output = run_fun3(run_mean, lambda, datav);
  std::cout << "(template) These are the running weighted means of the vector:" << std::endl;
  for (auto it: output) {
    std::cout << it << " ";
  }  // end for
  std::cout << std::endl;

  return 0;
}

