////////////////////////////
// Test C++ functions
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  g++ -std=c++11 /Users/jerzy/Develop/Rcpp/fun_test.cpp -o test
// Run it in Mac terminal: ./test
// 
// Compile this C++ file under Windows using MinGW:
//  C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\fun_test.cpp -o test

#include "fun_test.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <string>
#include <numeric>


// The function calc_unique() returns a vector with the unique
// elements of the integer input vector.
// It uses the STL set structure, without any explicit loop.
// It extracts the unique elements in two steps.
// First it calls a constructor for an unordered_set and copies 
// the input vector into a set.
// The STL set structure contains only unique elements, so any 
// duplicate elements are ignored.
// Then it calls a constructor for the output vector and copies 
// the set structure elements into it.
// The STL set can't be exported to Rcpp but can be used internally.
// https://www.techiedelight.com/convert-set-vector-cpp/
// https://codeyarns.com/2010/07/16/c-stl-inserting-vector-into-set/
// [[Rcpp::export]]
std::vector<int> calc_unique(std::vector<int> vectorv) {
  
  // Define unordered_set and copy the input vector into it.
  // The set contains only unique elements.
  std::unordered_set<int> setv(vectorv.begin(), vectorv.end());
  // Define output vector and copy the set into it.
  std::vector<int> output(setv.begin(), setv.end());
  
  return output;
  
}  // end calc_unique


// The function count_calls() counts the number of times it was called.
// It creates a static integer variable itv, and advances it every 
// time it is called.
// The variable itv is static so it remains alive outside the scope of 
// count_calls(), between the calls to count_calls().
// The function count_calls() is defined as a static function as 
// an illustration - it doesn't have to be static.
// A static function in C is only visible to those functions in the same source file. 
// Making a function static limits its scope to functions from the same source file. 
// https://www.cprogramming.com/tutorial/statickeyword.html
// [[Rcpp::export]]
void count_calls(unsigned int seedv) {
  
  // itv is static so it's initialized only once the first time count_calls() is called.
  static int itv = seedv;
  
  std::cout << "The function count_calls() was called " << itv << " times." << std::endl;
  
  itv++;
  
}  // end count_calls

