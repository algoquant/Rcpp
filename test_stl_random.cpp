////////////////////////////
// Test C++ random number code using Standard Template Library (STL)
////////////////////////////

// Compile this C++ file using MinGW:
// C:\Rtools\mingw_64\bin\g++ -std=c++11 C:\Develop\R\Rcpp\test_stl_random.cpp -o test

using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <time.h>


// The function count_er() counts the number of times it was called.
// It ceates a static integer variable cou_nt, and advances it every 
// time it is called.
// cou_nt is static so it remains alive outside the scope of count_er(),
// between calls to count_er().
// The function count_er() is defined as a static function as 
// an illustration - it doesn't have to be static.
// A static function in C is only visible to those functions in the same source file. 
// Making a function static limits its scope to functions from the same source file. 
// https://www.cprogramming.com/tutorial/statickeyword.html
// [[Rcpp::export]]
static void count_er(int see_d=1) {
  
  // cou_nt is static so it's initialized only once the first time count_er() is called.
  static int cou_nt = see_d;
  
  std::cout << "The function count_er() was called " << cou_nt << " times." << std::endl;
  
  cou_nt++;
  
}  // end count_er



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
std::uint32_t calc_seed(int num_el = 1) {
  
  // Produce random seeds using random_device - but random_device doesn't work under Windows.
  // std::random_device de_vice;
  // std::seed_seq seq_uence{de_vice()};
  // 
  // Or produce a single random seed
  // Get time in seconds.
  time_t see_d;
  time(&see_d);
  
  // Return last 2 seconds.
  // return (see_d % 100);
  return (see_d);
  
}  // end calc_seed



// The function calc_random_int() calculates a vector of random double numbers,
// between the values mi_n and ma_x.
// It uses the STL algorithm std::generate().
// The STL algorithm std::generate() is similar to the R functional apply(),
// but it can only call unary functions that don't accept any arguments.
// The difference between std::transform() and std::generate() is that 
// std::generate() only calls functions which return a value, while 
// std::transform() can also call functions which perform an operation.
// https://www.cprogramming.com/tutorial/statickeyword.html
// [[Rcpp::export]]
std::vector<int> calc_random_int(size_t num_el, int see_d=1, int mi_n=0, int ma_x=100) {

  // Define the type of random distribution (uniform, Gaussian, etc.)
  std::uniform_int_distribution<int> random_dist(mi_n, ma_x);
  // std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
  
  // Define the random number calculation engine (generator)
  // The random generator is static to ensure it's initialized only once.
  static std::default_random_engine random_generator{see_d};
  // Or
  // static std::mt19937 random_generator(see_d);
  
  // Define the vector of random numbers
  std::vector<int> da_ta(num_el);
  // Populate the vector of random numbers
  std::generate(da_ta.begin(), da_ta.end(),
                // Lambda function generates a single integer random number - it doesn't take any arguments
                [&]() { return random_dist(random_generator); }
                // [&]() { return ((random_dist(random_generator) % (ma_x - mi_n)) + mi_n); }
  );

  return da_ta;

}  // end calc_random_int



// The function calc_random_int2() calculates a vector of random integer numbers,
// between the values mi_n and ma_x.
// It uses the STL algorithm std::generate().
// The STL algorithm std::generate() is similar to the R functional apply().
// https://www.cprogramming.com/tutorial/statickeyword.html
// [[Rcpp::export]]
std::vector<int> calc_random_int2(size_t num_el, int see_d=1, int mi_n=0, int ma_x=100) {
  
  // Define the vector of random numbers
  std::vector<int> da_ta(num_el);
  // Set the seed
  srand(see_d);
  // Populate the vector of random numbers
  // std::generate(da_ta.begin(), da_ta.end(), rand);
  std::generate(da_ta.begin(), da_ta.end(),
                // Lambda function generates a single integer random number.
                // The brackets [] are used to pass in variables from the outer scope of the lambda function.
                // The "&" passes the outer scope variables by reference.
                // The brackets () are used to pass in variables from the std::generate() STL algorithm.
                [&ma_x, &mi_n]() { return ((rand() % (ma_x - mi_n)) + mi_n); });
  
  return da_ta; 
  
}  // end calc_random_int2



// The function calc_random_double() calculates a vector of random double numbers,
// between the values mi_n and ma_x.
// It uses the STL algorithm std::generate().
// https://www.cprogramming.com/tutorial/statickeyword.html
// [[Rcpp::export]]
std::vector<double> calc_random_double(size_t num_el, int see_d=1, double mi_n=0.0, double ma_x=1.0) {
  
  // Define the type of random distribution (uniform, Gaussian, etc.)
  // The random generator is static to ensure it's initialized only once.
  std::uniform_real_distribution<double> random_dist(mi_n, ma_x);
  
  // Define the random number calculation engine (generator)
  // The random generator is static to ensure it's initialized only once.
  static std::default_random_engine random_generator{see_d};
  // Or
  // static std::mt19937 random_generator(see_d);
  
  // Define the vector of random numbers
  std::vector<double> da_ta(num_el);
  // Populate the vector of random numbers
  std::generate(da_ta.begin(), da_ta.end(),
                // Lambda function generates a single integer random number - it doesn't take any arguments
                // The brackets [] are used to pass in variables from the outer scope of the lambda function.
                // The "&" passes the outer scope variables by reference.
                // The brackets () are used to pass in variables from the std::generate() STL algorithm.
                [&]() { return random_dist(random_generator); });
  
  return da_ta;
  
}  // end calc_random_double






int main() {
  
	// int in_put;
	size_t num_el;
  // std::string stri_ng;
  
  std::cout << "Enter number of random numbers: ";
  std::cin >> num_el;
  
  // Test for count_er().
  // std::cout << "The function count_er() was called as follows:" << std::endl;
  // for (int ele_ment = 0; ele_ment < num_el; ele_ment++) {
  //   // std::cout << "count_er was called as follows:" << std::endl;
  //   count_er();
  // }  // end for
  
  // Test for time().
  // std::cout << "Time in seconds now:" << time(0) << std::endl;
  // Or
  // time_t see_d;
  // time(&see_d);
  // std::cout << "Time in seconds now:" << see_d << std::endl;
  
  // Test for random_device.
  // std::random_device de_vice;
  // std::cout << "The random seed:" << de_vice() << std::endl;
  
  // Produce random seed sequence using system clock.
  // std::seed_seq seq_uence{time(0)};
  // std::vector<std::uint32_t> seed_s(num_el);
  // std::vector<std::uint32_t> seed_s = calc_seeds(3);
  // seq_uence.generate(seed_s.begin(), seed_s.end());
  // std::vector<std::uint32_t> seed_s = calc_seeds(num_el);
  // std::cout << "Random seeds:" << std::endl;
  // for (auto ele_ment : seed_s) {
    // std::cout << ele_ment << std::endl;
  // }  // end for

  // Don't know how to pass the seed_s into an object of class std::default_random_engine
  // std::default_random_engine random_generator{seed_s};  // this doesn't work
  // std::uniform_real_distribution<double> random_dist(0, 1);
  // std::default_random_engine random_generator{seed_s};
  // std::vector<double> da_ta(num_el);
  // std::generate(da_ta.begin(), da_ta.end(), [&]() { return random_dist(random_generator); });
  
  // Test for calc_seed()
  // int see_d = calc_seed();
  // std::cout << "Enter seed for random number generator: ";
  // std::cin >> see_d;
  // Test for random number generator
  // std::mt19937 random_generator(see_d);
  // std::default_random_engine random_generator{see_d};
  // std::uniform_real_distribution<double> random_dist(0, 1);
  std::vector<double> da_ta(num_el);
  // std::generate(da_ta.begin(), da_ta.end(), [&]() { return random_dist(random_generator); });
  // std::cout << "Uniform (0, 1): " << random_dist(random_generator);
  
  // Test for calc_random().
  // std::vector<int> da_ta = {7, 5, 16, 8, 16, 8};
  da_ta = calc_random_double(num_el, calc_seed());
  std::cout << "Random numbers:" << std::endl;
  for (auto ele_ment: da_ta) {
    std::cout << ele_ment << " ";
  }  // end for
  
}  // end main


