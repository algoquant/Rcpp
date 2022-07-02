////////////////////////////
// Test C++ try() and catch()
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  cd /Users/jerzy/Develop/Rcpp/
//  g++ -std=c++11 main_test.cpp fun_test.cpp -o test
// Run it in Mac terminal: ./test
// 
// Compile this C++ file under Windows using MinGW:
//  C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\main_test.cpp fun_test.cpp -o test

using namespace std;
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
  
  try {
    int age = 15;
    std::cout << "Enter your age: ";
    cin >> age;
    if (age >= 18) {
      std::cout << "Access granted - you are old enough." << std::endl;
    } else {
      throw (age);
    }
  }  // end try
  
  catch (int agerror) {
    std::cout << "Access denied - You must be at least 18 years old." << std::endl;
    std::cout << "Age is: " << agerror << std::endl;
  }  // end catch
  
  
  // std::stringstream ss (s);
  // ss >> i;
  // std::cout << "Integer is " << i << std::endl;

  // std::cout << std::endl;
  
}  // end main


