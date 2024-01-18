////////////////////////////
// Test C++ random number code using Standard Template Library (STL)
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/random.cpp")

using namespace std;
#include "random.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <time.h>


// The function calc_seed() calculates random seeds from the system clock,
// because random_device doesn't work under Windows.
// Generating true random numbers is more expensive than generating pseudorandom 
// numbers, so usually a single truly random number is generated and then passed 
// as a seed to a pseudorandom generator. 
// http://www.pcg-random.org/posts/simple-portable-cpp-seed-entropy.html
// https://stackoverflow.com/questions/34490599/c11-how-to-set-seed-using-random
// https://stackoverflow.com/questions/18880654/why-do-i-get-the-same-sequence-for-every-run-with-stdrandom-device-with-mingw
// https://stackoverflow.com/questions/22522829/how-to-properly-initialize-a-c11-stdseed-seq
// [[Rcpp::export]]
std::uint32_t calc_seed(int ndata) {
  
  // Produce random seeds using random_device - but random_device doesn't work under Windows.
  // std::random_device devicer;
  // std::seed_seq sequencer{devicer()};
  // 
  // Or produce a single random seed
  // Get time in seconds.
  time_t seedv;
  time(&seedv);
  
  // Return last 2 seconds.
  // return (seedv % 100);
  return (seedv);
  
}  // end calc_seed



// The function calc_random_int() calculates a vector of random double numbers,
// between the values minv and maxv.
// It uses the STL algorithm std::generate().
// The STL algorithm std::generate() is similar to the R functional apply(),
// but it can only call unary functions that don't accept any arguments.
// The difference between std::transform() and std::generate() is that 
// std::generate() only calls functions which return a value, while 
// std::transform() can also call functions which perform an operation.
// https://www.cprogramming.com/tutorial/statickeyword.html
// In R, run (example):
// calc_random_int(11, 1705352332, 0, 10)
// [[Rcpp::export]]
std::vector<int> calc_random_int(size_t ndata, unsigned int seedv, int minv, int maxv) {
  
  // Define the type of random distribution (uniform, Gaussian, etc.)
  std::uniform_int_distribution<int> random_dist(minv, maxv);
  // std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
  
  // Define the random number calculation engine (generator)
  // The random generator is static to ensure it's initialized only once.
  static std::default_random_engine random_generator{seedv};
  // Or
  // static std::mt19937 random_generator(seedv);
  
  // Define the vector of random numbers
  std::vector<int> datav(ndata);
  // Populate the vector of random numbers
  std::generate(datav.begin(), datav.end(),
                // Lambda function generates a single integer random number - it doesn't take any arguments
                [&]() { return random_dist(random_generator); }
                // [&]() { return ((random_dist(random_generator) % (maxv - minv)) + minv); }
  );
  
  return datav;
  
}  // end calc_random_int



// The function calc_random_int2() calculates a vector of random integer numbers,
// between the values minv and maxv.
// It uses the STL algorithm std::generate().
// The STL algorithm std::generate() is similar to the R functional apply().
// https://www.cprogramming.com/tutorial/statickeyword.html
// In R, run (example):
// calc_random_int2(11, 1705352332, 0, 10)
// [[Rcpp::export]]
std::vector<int> calc_random_int2(size_t ndata, int seedv, int minv, int maxv) {
  
  // Define the vector of random numbers
  std::vector<int> datav(ndata);
  // Set the seed
  srand(seedv);
  // Populate the vector of random numbers
  // std::generate(datav.begin(), datav.end(), rand);
  std::generate(datav.begin(), datav.end(),
                // Lambda function generates a single integer random number.
                // The brackets [] are used to pass in variables from the outer scope of the lambda function.
                // The "&" passes the outer scope variables by reference.
                // The brackets () are used to pass in variables from the std::generate() STL algorithm.
                [&maxv, &minv]() { return ((rand() % (maxv - minv)) + minv); });
  
  return datav; 
  
}  // end calc_random_int2



// The function calc_random_double() calculates a vector of random double numbers,
// between the values minv and maxv.
// It uses the STL algorithm std::generate().
// https://www.cprogramming.com/tutorial/statickeyword.html
// In R, run (example):
// calc_random_double(11, 1705352332, 0, 10)
// [[Rcpp::export]]
std::vector<double> calc_random_double(size_t ndata, unsigned int seedv, double minv, double maxv) {
  
  // Define the type of random distribution (uniform, Gaussian, etc.)
  // The random generator is static to ensure it's initialized only once.
  std::uniform_real_distribution<double> random_dist(minv, maxv);
  
  // Define the random number calculation engine (generator)
  // The random generator is static to ensure it's initialized only once.
  static std::default_random_engine random_generator{seedv};
  // Or
  // static std::mt19937 random_generator(seedv);
  
  // Define the vector of random numbers
  std::vector<double> datav(ndata);
  // Populate the vector of random numbers
  std::generate(datav.begin(), datav.end(),
                // Lambda function generates a single integer random number - it doesn't take any arguments
                // The brackets [] are used to pass in variables from the outer scope of the lambda function.
                // The "&" passes the outer scope variables by reference.
                // The brackets () are used to pass in variables from the std::generate() STL algorithm.
                [&]() { return random_dist(random_generator); });
  
  return datav;
  
}  // end calc_random_double

