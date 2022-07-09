////////////////////////////
// Code for testing C++ structures
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
// g++ -std=c++11 /Users/jerzy/Develop/Rcpp/test_struct.cpp -o test
// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_struct.cpp")

////////////////////////////
// C++ Structures (struct)
// https://www.w3schools.com/cpp/cpp_structs.asp

// Use STL
#include<string>
#include <iostream>

// Declare a structure named "car"
struct car {
  std::string brand;
  std::string model;
  int year;
};


// main Functions a car structure and store it in myCar1;
int main() {
  // Create a car structure and store it in myCar1;
  car myCar1; 
  myCar1.brand = "BMW";
  myCar1.model = "X5";
  myCar1.year = 1999;
  
  // Create another car structure and store it in myCar2;
  car myCar2;
  myCar2.brand = "Ford";
  myCar2.model = "Mustang";
  myCar2.year = 1969;
  
  // Print the structure members
  std::cout << myCar1.brand << " " << myCar1.model << " " << myCar1.year << "\n";
  std::cout << myCar2.brand << " " << myCar2.model << " " << myCar2.year << "\n";
  
  return 0;
}

