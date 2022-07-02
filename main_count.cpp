////////////////////////////
// Test C++ STL code for while loop and getline()
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  cd /Users/jerzy/Develop/Rcpp/
//  g++ -std=c++11 main_test.cpp fun_test.cpp -o test
// Run it in Mac terminal: ./test
// 
// Compile this C++ file under Windows using MinGW:
//  C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\main_test.cpp fun_test.cpp -o test

using namespace std;
#include "fun_test.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <string>
#include <numeric>

int main() {
  
  std::string stringv;
  // int input;
  // std::vector<int> inputs;
  // std::vector<int> inputs = {7, 5, 16, 8, 16, 8};
  
  // std::cout << "Enter string: ";
  // std::getline (std::cin, stringv);
  // std::cout << "String is " << stringv << std::endl;
  
  // std::cout << "Enter integer: ";
  // std::getline (std::cin, stringv);
  // 
  // if (stringv.empty()) {
  //   std::cout << "Input was empty" << std::endl;
  // } else {
  //   input = std::stoi(stringv);
  //   std::cout << "Integer is " << input << std::endl;
  // }  // end if
  
  stringv = "0";
  while (!stringv.empty()) {
    std::cout << "Enter integer: ";
    std::getline (std::cin, stringv);
    if (stringv.empty()) {
      std::cout << "Input was empty - bye!" << std::endl;
    } else {
      count_calls();
    }  // end if
  }  // end while
  
  // std::stringstream ss (s);
  // ss >> i;
  // std::cout << "Integer is " << i << std::endl;

  std::cout << std::endl;
  
}  // end main


