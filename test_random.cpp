////////////////////////////
// Test C++ random number code using Standard Template Library (STL)
////////////////////////////

// Compile this C++ file in Mac terminal using g++:
//  cd /Users/jerzy/Develop/Rcpp/
//  g++ -std=c++11 test_random.cpp random.cpp -o test
// Run it in Mac terminal: ./test
// 
// Compile this C++ file using MinGW:
// C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\test_random.cpp -o test

using namespace std;
#include "random.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <time.h>



int main() {
  
	// int input;
	size_t ndata;
  // std::string stri_ng;
  
  std::cout << "Enter number of random numbers: ";
  std::cin >> ndata;
  
  // Test for time().
  // std::cout << "Time in seconds now:" << time(0) << std::endl;
  // Or
  // time_t seedv;
  // time(&seedv);
  // std::cout << "Time in seconds now:" << seedv << std::endl;
  
  // Test for random_device.
  // std::random_device de_vice;
  // std::cout << "The random seed:" << de_vice() << std::endl;
  
  // Produce random seed sequence using system clock.
  // std::seed_seq seq_uence{time(0)};
  // std::vector<std::uint32_t> seed_s(ndata);
  // std::vector<std::uint32_t> seed_s = calc_seeds(3);
  // seq_uence.generate(seed_s.begin(), seed_s.end());
  // std::vector<std::uint32_t> seed_s = calc_seeds(ndata);
  // std::cout << "Random seeds:" << std::endl;
  // for (auto elemv : seed_s) {
    // std::cout << elemv << std::endl;
  // }  // end for

  // Don't know how to pass the seed_s into an object of class std::default_random_engine
  // std::default_random_engine random_generator{seed_s};  // this doesn't work
  // std::uniform_real_distribution<double> random_dist(0, 1);
  // std::default_random_engine random_generator{seed_s};
  // std::vector<double> datav(ndata);
  // std::generate(datav.begin(), datav.end(), [&]() { return random_dist(random_generator); });
  
  // Test for calc_seed()
  // int seedv = calc_seed();
  // std::cout << "Enter seed for random number generator: ";
  // std::cin >> seedv;
  // Test for random number generator
  // std::mt19937 random_generator(seedv);
  // std::default_random_engine random_generator{seedv};
  // std::uniform_real_distribution<double> random_dist(0, 1);
  std::vector<double> datav(ndata);
  // std::generate(datav.begin(), datav.end(), [&]() { return random_dist(random_generator); });
  // std::cout << "Uniform (0, 1): " << random_dist(random_generator);
  
  // Test for calc_random().
  // std::vector<int> datav = {7, 5, 16, 8, 16, 8};
  datav = calc_random_double(ndata, calc_seed());
  std::cout << "Random numbers:" << std::endl;
  for (auto elemv: datav) {
    std::cout << elemv << " ";
  }  // end for
  std::cout << std::endl;
  
}  // end main


