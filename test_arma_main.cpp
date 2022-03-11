// Compile this C++ file using MinGW:
// C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\test_arma_main.cpp -o test

#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <string>


arma::uvec calc_ranks(const arma::vec& vectorv) {
  return (arma::sort_index(arma::sort_index(vectorv)) + 1);
}  // end calc_ranks


std::vector<int> calc_unique(std::vector<int> vectorv) {
  
  // Define unordered_set and copy the input vector into it.
  // The set contains only unique elements.
  std::unordered_set<int> se_t(vectorv.begin(), vectorv.end());
  // Define output vector and copy the set into it.
  std::vector<int> output(se_t.begin(), se_t.end());
  
  return output;
  
}  // end calc_unique


int main() {
  
	int input;
  std::string stri_ng;
  std::vector<int> inputs;
  // std::vector<int> inputs = {7, 5, 16, 8, 16, 8};
  
  // std::cout << "Enter string: ";
  // std::getline (std::cin, stri_ng);
  // std::cout << "String is " << stri_ng << std::endl;
  
  // std::cout << "Enter integer: ";
  // std::getline (std::cin, stri_ng);
  // 
  // if (stri_ng.empty()) {
  //   std::cout << "Input was empty" << std::endl;
  // } else {
  //   input = std::stoi(stri_ng);
  //   std::cout << "Integer is " << input << std::endl;
  // }  // end if
  
  stri_ng = "0";
  while (!stri_ng.empty()) {
    std::cout << "Enter integer: ";
    std::getline (std::cin, stri_ng);
    if (stri_ng.empty()) {
      std::cout << "Input was empty" << std::endl;
    } else {
      input = std::stoi(stri_ng);
      inputs.push_back(input);
    }  // end if
  }  // end while
  
  // std::stringstream ss (s);
  // ss >> i;
  // std::cout << "Integer is " << i << std::endl;
  
  
	std::vector<int> uniqu_e = calc_unique(inputs);

	std::cout << "The sum of the inputs = " << std::accumulate(inputs.begin(), inputs.end(), 0) << std::endl;
	// std::cout << "The sum of the inputs = " << sum_it(inputs) << std::endl;
	
	std::cout << "These are the unique elements:" << std::endl;
	for (auto &ele_ment: uniqu_e) {
		std::cout << ele_ment << " ";
	}  // end for
	
	// return sum_it(inputs);
	
}  // end main


