////////////////////////////////////////////////////////////
// Utility functions using Armadillo and the Standard Template Library (STL)
////////////////////////////////////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/util_fun.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
// Use STL
// using namespace std;
// #include <functional>
// #include <string>
// #include <map>
// #include <iostream>
// [[Rcpp::plugins(cpp11)]]


////////////////////////////////////////////////////////////
// Functions miscellaneous
////////////////////////////////////////////////////////////

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


void print_it(int a,int b)
{
  std::cout << "first arg: " << a << std::endl;
  std::cout << "second arg: " << b << std::endl;
}




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




////////////////////////////////////////////////////////////
// Iterators
////////////////////////////////////////////////////////////


// The function remove_dup() removes consecutive duplicate elements 
// from the input vector of strings.
// It doesn't remove all duplicate elements.  
// It doesn't remove duplicate elements which don't neighbor each other.
// It uses the STL algorithm std::unique().
// [[Rcpp::export]]
std::vector<std::string> remove_dup(std::vector<std::string> stringv) {
  
  int ndata = stringv.size();
  // Define vector iterator
  std::vector<std::string>::iterator stringit;
  // Define vector of output strings
  std::vector<std::string> outv = stringv;
  outv.reserve(ndata);
  
  // Remove consecutive duplicate elements
  stringit = std::unique(outv.begin(), outv.end());
  // Resize the output vector
  outv.resize(std::distance(outv.begin(), stringit));
  
  return outv;
  
}  // end remove_dup



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



// The function select_it() selects elements of a vector using an STL iterator.
// [[Rcpp::export]]
void select_it(std::vector<double> tseries, int shiftv) {
  
  // Define STL iterator
  std::vector<double>::iterator itv = tseries.begin();
  // Select element from the front of the vector
  itv = next(itv, shiftv-1);
  
  // Print iterator value
  std::cout << "The " << shiftv << "th ";
  std::cout << "vector element from the front is = "; 
  std::cout << *itv << " "; 
  std::cout << std::endl;

  // Select element from the end of the vector
  itv = tseries.end();
  itv = prev(itv, shiftv);

  // Print iterator value
  std::cout << "The " << shiftv << "th ";
  std::cout << "vector element from the end is = "; 
  std::cout << *itv << " "; 
  std::cout << std::endl;
  
  // return *itv;
}  // end select_it



////////////////////////////////////////////////////////////
// STL set structure
// https://thispointer.com/stdset-tutorial-part-1-set-usage-details-with-default-sorting-criteria/
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// The function not_dup() returns a Boolean vector with TRUE indicating 
// that an element is not a duplicate of an element with a smaller subscript.
// It uses the STL set structure.
// The method .insert() adds an element to a set.
// The .insert().second is a Boolean equal to TRUE if the element is new, 
// i.e it's not already in the set.
// https://thispointer.com/different-ways-to-insert-elements-in-an-unordered_set-in-c11/
// [[Rcpp::export]]
std::vector<bool> not_dup(std::vector<int> tseries) {
  
  // int ndata = tseries.size();
  // Define iterator for input vector
  // std::vector<int>::iterator inpit;
  // Define Boolean output vector
  std::vector<bool> outv(tseries.size());
  // Define Boolean iterator for output vector
  std::vector<bool>::iterator notit = outv.begin();
  // Define unordered_set with unique elements
  std::unordered_set<int> setv;
  
  // Copy the elements of the input vector into a set
  for (auto inpit: tseries) {
    // for (auto inpit = tseries.begin(); inpit != tseries.end(); ++inpit) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    *notit++ = setv.insert(inpit).second;
  }  // end for
  
  return outv;
  
}  // end not_dup


// The function calc_unique() returns a vector with the unique
// elements of the numeric input vector.
// It uses the STL set structure, without any explicit loop.
// It extracts the unique elements in two steps.
// First it calls a constructor for an unordered_set and copies 
// the input vector into a set.
// The STL set structure contains only unique elements, so any 
// duplicate elements are ignored.
// Then it calls a constructor for the output vector and copies 
// the set structure elements into it.
// The STL set can't be exported to Rcpp but can be used internally.
// The function calc_unique() is slightly faster than the function 
// calc_unique_loop() which uses an explicit loop.
// https://www.techiedelight.com/convert-set-vector-cpp/
// https://codeyarns.com/2010/07/16/c-stl-inserting-vector-into-set/
// [[Rcpp::export]]
std::vector<double> calc_unique(std::vector<double> tseries) {
  
  // Define unordered_set and copy the input vector into it.
  // The set contains only unique elements.
  std::unordered_set<double> setv(tseries.begin(), tseries.end());
  // Define output vector and copy the set into it.
  std::vector<double> outv(setv.begin(), setv.end());
  
  return outv;
  
}  // end calc_unique



// The function calc_unique_loop() returns a vector with the unique
// elements of the numeric input vector.
// It uses the STL set structure, and an explicit loop.
// The method .insert() adds an element to a set.
// The .insert().second is a Boolean equal to TRUE if the element is new, 
// i.e it's not already in the set.
// The function calc_unique_loop() is about as fast as the R function unique().
// https://www.techiedelight.com/remove-duplicates-vector-cpp/
// [[Rcpp::export]]
std::vector<double> calc_unique_loop(std::vector<double> tseries) {
  
  // Define iterator for input vector
  // std::vector<double>::iterator inpit;
  // Define output vector
  std::vector<double> outv;
  // Define unordered_set with unique elements
  std::unordered_set<double> setv;
  
  // Copy the elements of the input vector into a set
  for (auto inpit: tseries) {
    // for (auto inpit = tseries.begin(); inpit != tseries.end(); ++inpit) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    if (setv.insert(inpit).second)
      outv.push_back(inpit);
  }  // end for
  
  return outv;

}  // end calc_unique_loop


// The function calc_unique_int() returns a vector with the unique
// elements of the integer input vector.
// It uses the STL set structure.
// [[Rcpp::export]]
std::vector<int> calc_unique_int(std::vector<int> tseries) {
  
  // Define unordered_set and copy the input vector into it.
  // The set contains only unique elements.
  std::unordered_set<int> setv(tseries.begin(), tseries.end());
  // Define output vector and copy the set into it.
  std::vector<int> outv(setv.begin(), setv.end());
  
  return outv;
  
}  // end calc_unique_int



// Uses an explicit loop.
// The function calc_unique_intt() returns a vector with the unique
// elements of the integer input vector.
// It uses the STL set structure, and an explicit loop.
// The method .insert() adds an element to a set.
// The .insert().second is a Boolean equal to TRUE if the element is new, 
// i.e if it's not already in the set.
// The function calc_unique() is about as fast as the R function unique().
// [[Rcpp::export]]
std::vector<int> calc_unique_intt(std::vector<int> tseries) {
  
  // Define iterator for input vector
  // std::vector<int>::iterator inpit;
  // Define output vector
  std::vector<int> outv;
  // Define unordered_set with unique elements
  std::unordered_set<int> setv;
  
  // Copy the elements of the input vector into a set
  for (auto inpit: tseries) {
    // for (auto inpit = tseries.begin(); inpit != tseries.end(); ++inpit) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    if (setv.insert(inpit).second)
      outv.push_back(inpit);
  }  // end for
  
  return outv;
  
}  // end calc_unique_intt



// The function calc_unique_str() returns a vector with the unique
// elements of the input vector of strings.
// It uses the STL set structure.
// It extracts the unique elements in two steps.
// First it copies the elements of the input vector into a set.
// The STL set structure contains only unique elements, so any 
// duplicate elements are ignored.
// Then it copies the unique set elements to the output vector.
// The STL set can't be exported to Rcpp but can be used internally.
// The function calc_unique_str() is slower than the R function unique().
// https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
// [[Rcpp::export]]
std::vector<std::string> calc_unique_str(std::vector<std::string> stringv) {

  int ndata = stringv.size();
  // Define vector iterator
  // std::vector<std::string>::iterator itv;
  // Define vector of output strings
  std::vector<std::string> outv;
  outv.reserve(ndata);
  // Define unordered_set of strings
  std::unordered_set<std::string> setv;
  // Define unordered_set iterator
  // std::unordered_set<std::string>::iterator uniquev;
  
  // Copy the elements of the input vector into a set
  for (auto stringit: stringv) {
    // for (auto stringit = stringv.begin(); stringit != stringv.end(); ++stringit) {
    setv.insert(stringit);
  }  // end for
  
  // Iterate through all the elements in a set and display the value.
  // std::cout << "These are the unique strings in the set:" << std::endl;
  // Copy the elements of the set to the output vector
  for (auto uniquev: setv) {
    // for (auto uniquev=setv.begin(); uniquev!=setv.end(); ++uniquev) {
    //   std::cout << *uniquev << std::endl;
    outv.push_back(uniquev);
  }  // end for
  
  return outv;

}  // end calc_unique_str




////////////////////////////////////////////////////////////
// STL map structure
// https://www.geeksforgeeks.org/map-associative-containers-the-c-standard-template-library-stl/
////////////////////////////////////////////////////////////

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


// The function map_it() maps strings to integers.
// It uses the STL map structure.
// [[Rcpp::export]]
int map_it(std::string stringv) {
  
  // Define map
  std::map <std::string, int> string_to_int;
  string_to_int["first"] = 1;
  string_to_int["second"] = 2;
  string_to_int["third"] = 3;

  std::cout << "Value: " << string_to_int[stringv] << std::endl;
  
  // Map the stringv argument to integer
  return string_to_int[stringv];
  // return string_to_int.find(stringv);
  
}  // end map_it


// Define function to map string to function.
// Doesn't work.
// void map_fun(std::string stringv) {
//   
//   // Define map
//   map <std::string, int> string_to_fun;
//   string_to_fun["first"] = first;
//   string_to_fun["second"] = second;
//   string_to_fun["third"] = third;
//   
//   string_to_fun[stringv];
//   
// }  // end map_fun


// The function calc_table() calculates a contingency table 
// of the integer input vector.
// It uses the STL map structure.
// [[Rcpp::export]]
std::vector<int> calc_table(std::vector<int> tseries) {

  // Define iterator for input vector
  // std::vector<int>::iterator inpit;
  // Define map iterator
  // std::map <int, int>::iterator itv;
  // Define map
  std::map <int, int> mapv;
  
  // Copy the elements of the input vector into a map
  for (auto inpit: tseries) {
    // for (auto inpit = tseries.begin(); inpit != tseries.end(); ++inpit) {
    mapv[inpit]++;
  }  // end for
  
  // Define contingency table
  std::vector<int> tablev;
  tablev.reserve(mapv.size());
  
  // Below are four different methods of copying the map to the output vector
  
  // Explicit for loop:
  // https://thispointer.com/how-to-copy-all-values-from-a-map-to-a-vector-in-c/
  for (auto mapit: mapv) {
    tablev.push_back(mapit.second);
  }  // end for
  
  // Old style loop:
  // for (mapit=mapv.begin(); mapit!=mapv.end(); ++mapit) {
  //   tablev.push_back(mapit->second);
  // }  // end for
  
  // Copy the map to the output vector using std::transform() STL algorithm.
  // The STL algorithm std::transform() is similar to the R functional apply().
  // std::transform(mapv.begin(), mapv.end(), std::back_inserter(tablev),
                 // Lambda function extracts the second field from the map row.
                 // The std::transform() STL algorithm passes in rows (pairs) of the mapv variable
                 // through the brackets ().
  //                [] (std::pair<int, int> const &rowv) {return rowv.second;}
  // );

  // Copy the map to the output vector using std::for_each() STL algorithm.
  // The STL algorithm std::for_each() is similar to a for loop.
  // std::for_each(mapv.begin(), mapv.end(),
                // Lambda function extracts the second field from the map row.
                // The std::for_each() STL algorithm passes in rows (pairs) of the mapv variable 
                // through the brackets ().
                // The brackets [] are used to pass in variables from the outer scope of the lambda function.
                // The "&" passes the outer scope variables by reference.
  //               [&tablev](std::pair<int, int>  rowv) {
  //                 tablev.push_back(rowv.second);
  //               }
  // );
  
  return tablev;

}  // end calc_table



////////////////////////////////////////////////////////////
// STL functionals
////////////////////////////////////////////////////////////

// Functionals are functions which accept functions as an argument. 
// STL algorithms are similar to R functionals.
// An example is the STL algorithm std::transform() - similar to the R functional apply().

// Some other STL functionals are:
// std::operator()
// std::accumulate()


////////////////////////////////////////////////////////////
//' Calculate the sum of the elements of a single-column \emph{time series},
 //' \emph{matrix}, or a \emph{vector} using \code{RcppArmadillo}.
 //' 
 //' @param \code{tseries} A single-column \emph{time series}, \emph{matrix}, or
 //'   a \emph{vector}.
 //'
 //' @return A single numeric value.
 //'
 //' @details
 //'   The function \code{sum_it()} calculates the sum of the elements of a
 //'   single-column \emph{time series}, \emph{matrix}, or a \emph{vector}.
 //'   
 //'   The function \code{sum_it()} calls the \code{STL} \code{C++} function
 //'   \code{std::accumulate()}.
 //' 
 //'   The function \code{sum_it()} is slower than the \code{R} function
 //'   \code{sum()}.
 //'   
 //' @examples
 //' \dontrun{
 //' # Create a vector of data
 //' datav <- rnorm(1e3)
 //' # Calculate the ranks of the elements using R code and RcppArmadillo
 //' all.equal(sum(datav), drop(HighFreq::sum_it(datav)))
 //' # Compare the speed of R code with RcppArmadillo
 //' library(microbenchmark)
 //' summary(microbenchmark(
 //'   Rcode=sum(datav),
 //'   Rcpp=HighFreq::sum_it(datav),
 //'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
 //' }
 //' 
 // [[Rcpp::export]]
 double sum_it(arma::vec tseries) {
   
   // Sum up elements of a vector using STL accumulate
   return std::accumulate(tseries.begin(), tseries.end(), 0.0);
   
   // Old code below
   // Sum up elements of a vector using STL iterator
   // int total = 0;
   
   // for (auto itv: tseries) {
   // Old-style loop
   // for (auto itv = tseries.begin(); itv != tseries.end(); ++itv) {
   // total += itv;
   // }  // end for
 }  // end sum_it



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




////////////////////////////////////////////////////////////
// STL algorithms
////////////////////////////////////////////////////////////

// The function sort_num() sorts the elements of a numeric vector.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<double> sort_num(std::vector<double> tseries) {
  
  std::sort(tseries.begin(), tseries.end());
  return tseries;
  
}  // end sort_num

// The function sort_string() sorts the elements of a vector of strings.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<std::string> sort_string(std::vector<std::string> tseries) {
  
  std::sort(tseries.begin(), tseries.end());
  return tseries;
  
}  // end sort_string


// Define a comparison function.
// [[Rcpp::export]]
bool is_greater(const std::vector<double> tseries, int i1, int i2) {
  return (tseries[i1] < tseries[i2]);
}  // end is_greater


// Define a comparison functor as a class.
class is_greater_functor {
  
public:
  // Overloaded operator - actual function
  bool operator()(std::vector<double> tseries, int i1, int i2) { return (tseries[i1] < tseries[i2]); }
  
};  // end is_greater_functor


// Define a comparison functor as a class.
class is_greater_functor2 {
  
  std::vector<double> tseries;
  
public:
  // Constructor
  is_greater_functor2(std::vector<double> tseries) : tseries(tseries) {}
  
  // Overloaded operator - actual function
  bool operator()(int i1, int i2) { return (tseries[i1] < tseries[i2]); }
  
};  // end is_greater_functor2



// The function sort_index() calculates the sort index of a numeric vector.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<int> sort_index(const std::vector<double> tseries) {
  
  int ndata = tseries.size();
  // Define index of integers along tseries
  std::vector<int> indeks(ndata);
  // Fill the vector indeks with a sequence of consecutive integers.
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(indeks.begin(), indeks.end(), 0);
  // Sort the index according to the order of tseries
  // is_greater_functor is_greater2;
  sort(indeks.begin(), indeks.end(), 
       // Lambda function defines sorting order.
       // The brackets [] are used to pass in variables from the outer scope of the lambda function.
       // The "&" passes the outer scope variables by reference.
       [&tseries](int i1, int i2) {return tseries[i1] < tseries[i2];}
       // Or call functor is_greater_functor2() - elegant but very slow!!!
       // is_greater_functor2(tseries)
       // Or call function is_greater() using std::bind() - very slow!!!
       // std::bind(is_greater, tseries, std::placeholders::1, std::placeholders::2)
       // Call function is_greater() using lambda function - also very slow!!!
       // [&tseries](int i1, int i2) {return is_greater(tseries, i1, i2);}
       // Or call function is_greater() using std::bind() - very slow!!!
       // std::bind(is_greater2, tseries, std::placeholders::1, std::placeholders::2)
  );
  
  return indeks;
  
}  // end sort_index


// The function calc_unique_sort() returns a vector with the unique
// elements of the input vector of strings.
// It uses the STL algorithms std::sort() and std::unique().
// It extracts the unique elements in two steps.
// First it sorts the elements of the input vector, and then 
// it removes consecutive duplicate elements.
// The function calc_unique_sort() is slower than the R function unique().
// [[Rcpp::export]]
std::vector<std::string> calc_unique_sort(std::vector<std::string> stringv) {
  
  int ndata = stringv.size();
  // Define vector iterator
  std::vector<std::string>::iterator stringit;
  // Define vector of output strings
  std::vector<std::string> outv = stringv;
  outv.reserve(ndata);
  // Sort the elements of the input vector
  std::sort(outv.begin(), outv.end());
  // Remove consecutive duplicate elements
  stringit = std::unique(outv.begin(), outv.end());
  // Resize the output vector
  outv.resize(std::distance(outv.begin(), stringit));
  
  return outv;
  
}  // end calc_unique_sort



