////////////////////////////
// Update the class state variable in a loop.
// The state variable is persistent between function calls.
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  cd /Users/jerzy/Develop/Rcpp/
//  g++ -std=c++11 main_class.cpp -o test

using namespace std;
#include <iostream>

// Functor class for multiplying a number by a factor.
class update_state {
  
  // The private state variable is persistent between function calls.
private:
  double privstate;
  
public:
  // The public state variable is also persistent between function calls.
  double pubstate = 0;
  
  // Constructor
  update_state(double input) : privstate(input) {}
  
  // Overloaded operator - the actual function
  double operator()(double x) {
    privstate = privstate + x;
    return privstate;
    }  // end operator
  
  // Add to private state variable
  double add_to_private(double x) {
    privstate = privstate + x;
    return privstate;
  }  // end add_to_private
  
  // Add to public state variable
  void add_to_public(double x) {
    pubstate = pubstate + x;
  }  // end add_to_public
  
};  // end update_state


int main() {
  
  double inputv;
  double privstate;
  double pubv;
  std::cout << "Perform loop and add numbers to state variable." << std::endl;
  std::cout << "Enter the initial private state variable: ";
  std::cin >> inputv;
  // Create the instance update_it of the functor class update_state
  update_state update_it(inputv);
  std::cout << "The initial private state variable = " << inputv << std::endl;

  while (inputv > 0) {
    std::cout << "Enter a number to add to the private state variable (terminate with zero): ";
    std::cin >> inputv;
    std::cout << "Enter a number to add to the public state variable: ";
    std::cin >> pubv;
    update_it.add_to_public(pubv);
    // update_it.pubstate = update_it.pubstate + pubv;
    if (inputv > 0) {
      // privstate = update_it(inputv);
      privstate = update_it.add_to_private(inputv);
      std::cout << "The private state variable = " << privstate << std::endl;
      std::cout << "The public state variable = " << update_it.pubstate << std::endl;
    } else {
      std::cout << "Input was zero - bye!" << std::endl;
    }  // end if
  }  // end while

}  // end main

