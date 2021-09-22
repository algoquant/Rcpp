////////////////////////////
// Test C++ STL code for while loop and getline()
////////////////////////////

// Compile this C++ file using MinGW:
// C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\test_stl_main.cpp -o test

#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <string>


std::vector<int> calc_unique(std::vector<int> vec_tor) {
  
  // Define unordered_set and copy the input vector into it.
  // The set contains only unique elements.
  std::unordered_set<int> se_t(vec_tor.begin(), vec_tor.end());
  // Define output vector and copy the set into it.
  std::vector<int> out_put(se_t.begin(), se_t.end());
  
  return out_put;
  
}  // end calc_unique


int main() {
  
	int in_put;
  std::string stri_ng;
  std::vector<int> in_puts;
  // std::vector<int> in_puts = {7, 5, 16, 8, 16, 8};
  
  // std::cout << "Enter string: ";
  // std::getline (std::cin, stri_ng);
  // std::cout << "String is " << stri_ng << std::endl;
  
  // std::cout << "Enter integer: ";
  // std::getline (std::cin, stri_ng);
  // 
  // if (stri_ng.empty()) {
  //   std::cout << "Input was empty" << std::endl;
  // } else {
  //   in_put = std::stoi(stri_ng);
  //   std::cout << "Integer is " << in_put << std::endl;
  // }  // end if
  
  stri_ng = "0";
  while (!stri_ng.empty()) {
    std::cout << "Enter integer: ";
    std::getline (std::cin, stri_ng);
    if (stri_ng.empty()) {
      std::cout << "Input was empty" << std::endl;
    } else {
      in_put = std::stoi(stri_ng);
      in_puts.push_back(in_put);
    }  // end if
  }  // end while
  
  // std::stringstream ss (s);
  // ss >> i;
  // std::cout << "Integer is " << i << std::endl;
  
  
	std::vector<int> uniqu_e = calc_unique(in_puts);

	std::cout << "The sum of the inputs = " << std::accumulate(in_puts.begin(), in_puts.end(), 0) << std::endl;
	// std::cout << "The sum of the inputs = " << sum_it(in_puts) << std::endl;
	
	std::cout << "These are the unique elements:" << std::endl;
	for (auto ele_ment: uniqu_e) {
		std::cout << ele_ment << " ";
	}  // end for
	

}  // end main


