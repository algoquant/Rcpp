////////////////////////////
// Test functions using Armadillo and the Standard Template Library (STL)
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_fun.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;
// Use STL
// using namespace std;
// #include <functional>
// #include <string>
// #include <map>
// #include <iostream>
// [[Rcpp::plugins(cpp11)]]


////////////////////////////
// Misc functions

// [[Rcpp::export]]
std::vector<double> decode_it(Rcpp::List encodel) {
  
  // Extract vector of encoded data
  std::vector<double> codev = encodel["data"];
  // Extract vector of data counts (repeats)
  std::vector<int> countv = encodel["counts"];
  // Define the output decoded vector
  std::vector<double> decodev;
  decodev.reserve(std::accumulate(countv.begin(), countv.end(), 0));
  
  // Perform loop over the codev and countv vectors
  for (std::size_t it = 0; it < codev.size(); it++) {
    for (int j = 0; j < countv[it]; j++) {
      decodev.push_back(codev[it]);
    }  // end for
  }  // end for
  
  return decodev;
  
}  // end decode_it

arma::uvec calc_endpoints(arma::uword length,  // The length of the sequence
                          arma::uword step = 1,  // The number of periods between neighboring end points
                          arma::uword stub = 0,  // The first non-zero end point
                          bool stubs = true) {  // Include stub intervals?
  
  // Number of initial end points
  arma::uword numpts = length / step + 3;
  // Define the end points
  arma::uvec endp;
  endp.zeros(numpts);
  // Define the last end point
  arma::uword lastp = length - 1;
  
  // Calculate the initial end points - including extra end points at the end
  if (stub == 0) {
    for (arma::uword it = 0; it < numpts; ++it) {
      endp[it] = it*step;
    }  // end for
  } else if ((stub > 0) & (stubs)) {
    for (arma::uword it = 1; it < numpts; ++it) {
      endp[it] = stub + (it-1)*step;
    }  // end for
  } else {
    for (arma::uword it = 0; it < numpts; ++it) {
      endp[it] = stub + it*step;
    }  // end for
  }  // end if
  // std::cout << "endp = " << arma::conv_to<rowvec>::from(endp) << std::endl;
  
  // arma::uvec endp = arma::regspace<uvec>(stub, step, lastp + step);
  // Find the index of the largest element of endp which is less than lastp
  arma::uword endpp = 0;
  for (arma::uword it = 0; endp[it] < lastp; ++it) {
    endpp++;
  }  // end for
  // std::cout << "endpp = " << endpp << std::endl;
  
  // Trim the initial end points at the end - remove extra end points at the end
  // Subset endp to include the smallest element of endp which is equal to or greater than lastp
  endp = endp.subvec(0, endpp);
  
  // Set the stub intervals at the end
  if (stubs) {
    // Include stub intervals
    // Set the last end point to lastp - last element of endp
    endp[endpp] = lastp;
  } else {
    // Do not include the last end point - no stub interval at the end
    // Exclude the last element greater than lastp
    if (endp[endpp] > lastp) {
      endp = endp.subvec(0, endpp-1);
    }  // end if
  }  // end if
  
  return endp;
  
}  // end calc_endpoints


arma::mat diffit(const arma::mat& tseries, arma::sword lagg = 1, bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_rows-1);
  arma::uword ncols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    // Matrix difference without padding
    arma::mat diffmat = (tseries.rows(lagg, nrows) - tseries.rows(0, nrows - lagg));
    if (pad_zeros) {
      // Pad diffmat with zeros at the front
      return arma::join_cols(arma::zeros<mat>(lagg, ncols), diffmat);
    } else {
      // Don't pad the output
      return diffmat;
    }  // end if pad_zeros
  } else {
    // Negative lag
    // Matrix difference without padding
    arma::mat diffmat = (tseries.rows(0, nrows + lagg) - tseries.rows(-lagg, nrows));
    if (pad_zeros) {
      // Pad diffmat with zeros at the back
      return arma::join_cols(diffmat, arma::zeros<mat>(-lagg, ncols));
    } else {
      // Don't pad the output
      return diffmat;
    }  // end if pad_zeros
  }  // end if lagg
  
}  // end diffit



arma::mat calc_var_ag(const arma::mat& tseries, 
                      arma::uword step = 1) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return arma::var(diffit(tseries, 1, false));
  else {
    // Allocate aggregations, end points, and variance.
    arma::uword nrows = tseries.n_rows;
    arma::mat aggs;
    arma::uvec endp;
    // The number of rows is equal to step so that loop works for stub=0
    arma::mat vars(step, tseries.n_cols);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub < step; stub++) {
      endp = calc_endpoints(nrows, step, stub, false);
      // endp = arma::regspace<uvec>(stub, step, nrows + step);
      // endp = endp.elem(find(endp < nrows));
      aggs = tseries.rows(endp);
      vars.row(stub) = arma::var(diffit(aggs, 1, false));
    }  // end for
    return mean(vars);
  }  // end if
  
}  // end calc_var_ag


// [[Rcpp::export]]
arma::mat calc_hurst(const arma::mat& tseries, 
                     const arma::vec& aggv) {

  // If only single agg value then calculate the Hurst exponents from a single data point
  if (aggv.n_rows == 1) {
    return 0.5*arma::log(calc_var_ag(tseries, aggv[0])/calc_var_ag(tseries, 1))/log(aggv[0]);
  }  // end if
  
  // Allocate the objects
  arma::uword nrows = aggv.n_rows;
  arma::mat volv(nrows, tseries.n_cols, fill::zeros);
  
  // Calculate the log volatility at the agg points
  for (arma::uword it=0; it < nrows; it++) {
    volv.row(it) = 0.5*arma::log(calc_var_ag(tseries, aggv[it]));
  }  // end for
  
  // Calculate the log of the agg points
  arma::mat agglog = arma::log(aggv);
  
  // Calculate the Hurst exponents from the regression slopes
  arma::mat varagg = arma::var(agglog);
  return (arma::cov(volv, agglog).t())/varagg[0];
  
}  // end calc_hurst



////////////////////////////
// STL functionals
// STL has functionals like std::transform()
// Functionals are functions which accept functions as an argument. 

// Some other STL functionals
// std::operator()
// std::accumulate()


// Helper function extracts the second elements from a pair.
int get_val(std::pair<std::string, int> const &pair) {return pair.second;}

// The function map_out() copies elements from a map to a vector. 
// It uses the STL functional std::transform()
// https://thispointer.com/how-to-copy-all-values-from-a-map-to-a-vector-in-c/
// [[Rcpp::export]]
std::vector<int> map_out(std::string stringv = "blah") {
  
  // Define map
  std::map <std::string, int> mapv;
  mapv["first"] = 1;
  mapv["second"] = 2;
  mapv["third"] = 3;
  
  std::vector<int> tseries;
  tseries.reserve(mapv.size());
  
  // std::cout << "Value: " << stringv << std::endl;
  // Copy all values from a map to vector using transform() and a function pointer
  std::transform(mapv.begin(), mapv.end(), std::back_inserter(tseries), &get_val);
  
  return tseries;
  
}  // end map_out



// The function print_string() prints the elements of a vector of strings 
// using the STL algorithm std::for_each().
// [[Rcpp::export]]
void print_string(std::vector<std::string> tseries) {
  
  std::for_each(tseries.begin(), tseries.end(),
                // Lambda function prints a single string
                [](std::string stringv) { std::cout << stringv << ", "; });
  std::cout << std::endl;
  
}  // end print_string

// square_it() is a non-exported function which squares a double.
// It can be used by other functions.
double square_it(double x) { return x*x; }

// The STL functional std::transform() can be used to apply a function over a vector.
// The function square_vec() squares the elements of a numeric vector 
// by calling a lambda function using the functional std::transform().
// Lambda functions are anonymous functions which can be passed to functionals.
// [[Rcpp::export]]
std::vector<double> square_vec(const std::vector<double> tseries) {
  std::vector<double> outv(tseries.size());
  
  std::transform(tseries.begin(), tseries.end(), outv.begin(), 
                 // Pass function square_it() to functional std::transform()
                 // square_it
                 // Or pass a lambda function
                 [](double x) { return x*x; }
  ); // end std::transform()
  
  return outv;
}  // end square_vec



// Define a non-exported function which squares two doubles.
// It can be used by other functions.
// [[Rcpp::export]]
double square_two(double x, double y) { return (x*x + y*y); }

// The STL functional std::transform() can be used to apply a function over a vector.
// The function square_two_vec() squares and sums the elements of two numeric vectors 
// by calling a lambda function using the functional std::transform().
// Lambda functions are anonymous functions which can be passed to functionals.
// [[Rcpp::export]]
std::vector<double> square_two_vec(const std::vector<double> x, const std::vector<double> y) {
  std::vector<double> z(x.size());
  std::transform(x.begin(), x.end(), y.begin(), z.begin(),  
                 // Pass function square_two() to functional std::transform()
                 // square_two
                 // Or pass a lambda function
                 [](double x, double y) { return (x*x + y*y); }
  ); // end std::transform()
  return z;
}  // end square_two_vec


// Do these examples:
// https://stackoverflow.com/questions/24017617/using-functions-with-multiple-arguments-with-lapply-in-rcpp
// https://stackoverflow.com/questions/55715001/how-to-pass-lambda-function-in-a-c-functor



// The function match_it() reproduces the R function findInterval().
// The function match_it() matches its inputs with a vector of break points.
// tseries is a vector of inputs to be matched with the break points.
// breakv is a vector of break points.
// The matches are the indices of the break points closest to the tseries.
// [[Rcpp::export]]
std::vector<int> match_it(std::vector<double> tseries, std::vector<double> breakv) {
  
  // Allocate vector of matchv: the break points that match the tseries.
  std::vector<int> matchv(tseries.size());
  // Allocate iterators (pointers) to the tseries and matchv
  std::vector<int>::iterator matchit;
  std::vector<double>::iterator inpit, breakit;
  
  // Loop over the vector tseries and calculate the matchv
  for (inpit = tseries.begin(), matchit = matchv.begin(); inpit != tseries.end(); ++inpit, ++matchit) {
    // Find closest break point to the input
    breakit = std::upper_bound(breakv.begin(), breakv.end(), *inpit);
    // Calculate the index of the closest break point
    *matchit = std::distance(breakv.begin(), breakit);
  }  // end for
  
  return matchv;
  
}  // end match_it



// The function read_back() returns the reverse of its 
// input vector using a reverse_iterator.
// [[Rcpp::export]]
std::vector<int> read_back(std::vector<int> tseries) {
  
  // Define vectors
  std::vector<int> revec;
  
  // Initialize first value
  std::vector<int>::reverse_iterator revit;
  
  for (revit = tseries.rbegin(); revit != tseries.rend(); ++revit) {
    revec.push_back(*revit);
  }
  
  return revec;
  
}  // end read_back



// Create a functional that accepts a function as an argument.
// Pass the function as argument to functional.
// None of the functionals below work in Rcpp.

// void print_it( void(*func_arg)() ) {
// Or
// void print_it(std::function<void()> func_arg) {
//   func_arg();
// }  // end print_it

// void call_it(void func_arg())
// { func_arg(); }

// void call_it(int a,int b, void func_arg(int, int))
// { func_arg(a, b); }

// template <typename Callable>
// void call_it(Callable func_arg) {func_arg();}

void print_it(int a,int b)
{
  std::cout << "first arg: " << a << std::endl;
  std::cout << "second arg: " << b << std::endl;
}



// Define functions that print something
// [[Rcpp::export]]
void first() {
  std::cout << "first" << std::endl;
}  // end first

// [[Rcpp::export]]
void second() {
  std::cout << "second" << std::endl;
}  // end second

// [[Rcpp::export]]
void third() {
  std::cout << "third" << std::endl;
}  // end third




// Convert a string to a function name and run it
// use map::find for lookups

// typedef void (*func_ptr)(void);
// 
// std::map<std::string, func_ptrfunctions;
// 
// void fun1()
// {
//   std::std::cout << "fun1\n";
// }
// 
// void fun2()
// {
//   std::std::cout << "fun2\n";
// }
// 
// int main()
// {
//   functions["fun1"] = &fun1;
//   functions["fun2"] = &fun2;
//   
//   std::string name;
//   std::cin >name;
//   
//   functions[name](); // invoke
// }

