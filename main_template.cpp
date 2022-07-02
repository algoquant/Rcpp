////////////////////////////
// Test C++ templates
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  cd /Users/jerzy/Develop/Rcpp/
//  g++ -std=c++11 main_test.cpp fun_test.cpp -o test
// Run it in Mac terminal: ./test
// 
// Compile this C++ file under Windows using MinGW:
//  C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\main_test.cpp fun_test.cpp -o test

// using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <string>
#include <numeric>

// Define template function
template <class T>
T calc_max(T a, T b) {
  T result;
  result = (a>b)? a : b;
  return (result);
}

int main() {
  int i = 5, j = 6, k;
  double l = 10.2, m = 5.1, n;
  // Call template function for integer argument
  k = calc_max(i,j);
  std::cout << k << std::endl;
  // Call template function for float argument
  n = calc_max(l,m);
  std::cout << n << std::endl;
  return 0;
}

