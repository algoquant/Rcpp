////////////////////////////
// Functions to test C++ syntax using the Standard Template Library (STL)
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_stl.cpp")

#include <Rcpp.h>
using namespace Rcpp;
// Use STL
using namespace std;
// #include <functional>
// #include <string>
// #include <map>
// #include <iostream>
// [[Rcpp::plugins(cpp11)]]


////////////////////////////
// STL iterators (pointers)
// https://www.geeksforgeeks.org/iterators-c-stl/


// The function sum_it() sums up the elements of a vector using an STL iterator.
// In practcie it's better to simply use std::accumulate().
// [[Rcpp::export]]
int sum_it(std::vector<int> vec_tor) {
  
  // Define STL iterator
  // std::vector<int>::iterator it_er;
  
  // Sum up elements of a vector using STL iterator
  int total = 0;
  for (auto it_er: vec_tor) {
    // Old-style loop
    // for (auto it_er = vec_tor.begin(); it_er != vec_tor.end(); ++it_er) {
    total += it_er;
  }  // end for
  
  return total;
}  // end sum_it



// The function select_it() selects elements of a vector using an STL iterator.
// [[Rcpp::export]]
void select_it(std::vector<double> vec_tor, int shif_t) {
  
  // Define STL iterator
  std::vector<double>::iterator it_er = vec_tor.begin();
  // Select element from the front of the vector
  it_er = next(it_er, shif_t-1);
  
  // Print iterator value
  std::cout << "The " << shif_t << "th ";
  std::cout << "vector element from the front is = "; 
  std::cout << *it_er << " "; 
  std::cout << std::endl;

  // Select element from the end of the vector
  it_er = vec_tor.end();
  it_er = prev(it_er, shif_t);

  // Print iterator value
  std::cout << "The " << shif_t << "th ";
  std::cout << "vector element from the end is = "; 
  std::cout << *it_er << " "; 
  std::cout << std::endl;
  
  // return *it_er;
}  // end select_it




// The function encode_it() implements a run length encoder
// using a reverse_iterator.
// The run length encoding of a vector consists of a vector
// of consecutive input values and their counts (the number 
// of times the same value repeats in succession).
// Run-length encoding (RLE) is a data compression algorithm
// which encodes the data in two vectors: the consecutive data 
// values and their counts.  If a data value occurs several 
// times in succession then it is recorded only once and its 
// corresponding count is equal to the number of times it 
// occurs.
// Run-length encoding is different from a contingency table.
// [[Rcpp::export]]
List encode_it(std::vector<double> in_puts) {
  
  // Define vector of input data
  std::vector<double> da_ta;
  // Define vector of data counts (repeats)
  std::vector<int> count_s;
  
  // Define iterators
  std::vector<double>::iterator in_put = in_puts.begin();
  std::vector<int>::reverse_iterator cou_nt = count_s.rbegin();
  
  // Initialise the data
  double la_st = *in_put;
  da_ta.push_back(la_st);
  count_s.push_back(1);
  
  // Perform loop over in_puts using the iterator
  for (in_put = in_puts.begin() + 1; in_put != in_puts.end(); ++in_put) {
    if (la_st == *in_put) {
      (*cou_nt)++;
    } else {
      da_ta.push_back(*in_put);
      count_s.push_back(1);
      cou_nt = count_s.rbegin();
      la_st = *in_put;
    }  // end if
  }  // end for
  
  return List::create(
    _["data"] = da_ta,
    _["counts"] = count_s
  );
}  // end encode_it




// wippp
// The function decode_it() decodes a vector from its run length encoding.
// [[Rcpp::export]]
std::vector<double> decode_it(List in_puts) {
  
  // Define vector of input data
  std::vector<double> da_ta;
  // Define vector of data counts (repeats)
  // std::vector<int> count_s;
  
  // Define iterators
  // std::vector<double>::iterator in_put = in_puts.begin();
  // std::vector<int>::reverse_iterator cou_nt = count_s.rbegin();
  // 
  // // Initialise the data
  // double la_st = *in_put;
  // da_ta.push_back(la_st);
  // count_s.push_back(1);
  // 
  // // Perform loop over in_puts using the iterator
  // for (in_put = in_puts.begin() + 1; in_put != in_puts.end(); ++in_put) {
  //   if (la_st == *in_put) {
  //     (*cou_nt)++;
  //   } else {
  //     da_ta.push_back(*in_put);
  //     count_s.push_back(1);
  //     cou_nt = count_s.rbegin();
  //     la_st = *in_put;
  //   }  // end if
  // }  // end for
  
  return da_ta;
}  // end decode_it




////////////////////////////
// STL set structure
// https://thispointer.com/stdset-tutorial-part-1-set-usage-details-with-default-sorting-criteria/


// The function not_dup() returns a Boolean vector with TRUE indicating 
// that an element is a duplicate of an element with a smaller subscript.
// It uses the STL set structure.
// The method .insert() adds an element to a set.
// The .insert().second is a Boolean equal to TRUE if the element is new, 
// i.e it's not already in the set.
// https://thispointer.com/different-ways-to-insert-elements-in-an-unordered_set-in-c11/
// [[Rcpp::export]]
std::vector<bool> not_dup(std::vector<int> vec_tor) {
  
  // int num_el = vec_tor.size();
  // Define iterator for input vector
  // std::vector<int>::iterator ele_ment;
  // Define Boolean output vector
  std::vector<bool> out_put(vec_tor.size());
  // Define Boolean iterator for output vector
  std::vector<bool>::iterator is_dup = out_put.begin();
  // Define unordered_set with unique elements
  std::unordered_set<int> se_t;
  
  // Copy the elements of the input vector into a set
  for (auto ele_ment: vec_tor) {
    // for (auto ele_ment = vec_tor.begin(); ele_ment != vec_tor.end(); ++ele_ment) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    *is_dup++ = se_t.insert(ele_ment).second;
  }  // end for
  
  return out_put;
  
}  // end not_dup


// wippp
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
std::vector<double> calc_unique(std::vector<double> vec_tor) {
  
  // Define unordered_set and copy the input vector into it.
  // The set contains only unique elements.
  std::unordered_set<double> se_t(vec_tor.begin(), vec_tor.end());
  // Define output vector and copy the set into it.
  std::vector<double> out_put(se_t.begin(), se_t.end());
  
  return out_put;
  
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
std::vector<double> calc_unique_loop(std::vector<double> vec_tor) {
  
  // Define iterator for input vector
  // std::vector<double>::iterator ele_ment;
  // Define output vector
  std::vector<double> out_put;
  // Define unordered_set with unique elements
  std::unordered_set<double> se_t;
  
  // Copy the elements of the input vector into a set
  for (auto ele_ment: vec_tor) {
    // for (auto ele_ment = vec_tor.begin(); ele_ment != vec_tor.end(); ++ele_ment) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    if (se_t.insert(ele_ment).second)
      out_put.push_back(ele_ment);
  }  // end for
  
  return out_put;

}  // end calc_unique_loop



// The function calc_unique_int() returns a vector with the unique
// elements of the integer input vector.
// It uses the STL set structure.
// The method .insert() adds an element to a set.
// The .insert().second is a Boolean equal to TRUE if the element is new, 
// i.e it's not already in the set.
// The function calc_unique() is about as fast as the R function unique().
// [[Rcpp::export]]
std::vector<int> calc_unique_int(std::vector<int> vec_tor) {
  
  // Define iterator for input vector
  // std::vector<int>::iterator ele_ment;
  // Define output vector
  std::vector<int> out_put;
  // Define unordered_set with unique elements
  std::unordered_set<int> se_t;
  
  // Copy the elements of the input vector into a set
  for (auto ele_ment: vec_tor) {
    // for (auto ele_ment = vec_tor.begin(); ele_ment != vec_tor.end(); ++ele_ment) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    if (se_t.insert(ele_ment).second)
      out_put.push_back(ele_ment);
  }  // end for
  
  return out_put;
  
}  // end calc_unique_int



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
std::vector<std::string> calc_unique_str(std::vector<std::string> string_s) {

  int num_el = string_s.size();
  // Define vector iterator
  // std::vector<std::string>::iterator stri_ng;
  // Define vector of output strings
  std::vector<std::string> out_put;
  out_put.reserve(num_el);
  // Define unordered_set of strings
  std::unordered_set<std::string> se_t;
  // Define unordered_set iterator
  // std::unordered_set<std::string>::iterator uniqu_e;
  
  // Copy the elements of the input vector into a set
  for (auto stri_ng: string_s) {
    // for (auto stri_ng = string_s.begin(); stri_ng != string_s.end(); ++stri_ng) {
    se_t.insert(stri_ng);
  }  // end for
  
  // Iterate through all the elements in a set and display the value.
  // std::cout << "These are the unique strings in the set:" << std::endl;
  // Copy the elements of the set to the output vector
  for (auto uniqu_e: se_t) {
    // for (auto uniqu_e=se_t.begin(); uniqu_e!=se_t.end(); ++uniqu_e) {
    //   std::cout << *uniqu_e << std::endl;
    out_put.push_back(uniqu_e);
  }  // end for
  
  return out_put;

}  // end calc_unique_str




////////////////////////////
// STL map structure
// https://www.geeksforgeeks.org/map-associative-containers-the-c-standard-template-library-stl/


// The function map_it() maps strings to integers.
// It uses the STL map structure.
// [[Rcpp::export]]
int map_it(std::string stri_ng) {
  
  // Define map
  std::map <std::string, int> string_to_int;
  string_to_int["first"] = 1;
  string_to_int["second"] = 2;
  string_to_int["third"] = 3;

  std::cout << "Value: " << string_to_int[stri_ng] << std::endl;
  
  // Map the stri_ng argument to integer
  return string_to_int[stri_ng];
  // return string_to_int.find(stri_ng);
  
}  // end map_it


// Define function to map string to function.
// Doesn't work.
// void map_fun(std::string stri_ng) {
//   
//   // Define map
//   map <std::string, int> string_to_fun;
//   string_to_fun["first"] = first;
//   string_to_fun["second"] = second;
//   string_to_fun["third"] = third;
//   
//   string_to_fun[stri_ng];
//   
// }  // end map_fun


// The function calc_table() calculates a contingency table 
// of the integer input vector.
// It uses the STL map structure.
// [[Rcpp::export]]
std::vector<int> calc_table(std::vector<int> vec_tor) {

  // Define iterator for input vector
  // std::vector<int>::iterator ele_ment;
  // Define map iterator
  // std::map <int, int>::iterator cou_nt;
  // Define map
  std::map <int, int> ma_p;
  
  // Copy the elements of the input vector into a map
  for (auto ele_ment: vec_tor) {
    // for (auto ele_ment = vec_tor.begin(); ele_ment != vec_tor.end(); ++ele_ment) {
    ma_p[ele_ment]++;
  }  // end for
  
  // Define contingency table
  std::vector<int> ta_ble;
  ta_ble.reserve(ma_p.size());
  
  // Below are four different methods of copying the map to the output vector
  
  // Explicit for loop:
  // https://thispointer.com/how-to-copy-all-values-from-a-map-to-a-vector-in-c/
  for (auto cou_nt: ma_p) {
    ta_ble.push_back(cou_nt.second);
  }  // end for
  
  // Old style loop:
  // for (cou_nt=ma_p.begin(); cou_nt!=ma_p.end(); ++cou_nt) {
  //   ta_ble.push_back(cou_nt->second);
  // }  // end for
  
  // Copy the map to the output vector using std::transform() STL algorithm.
  // The STL algorithm std::transform() is similar to the R functional apply().
  // std::transform(ma_p.begin(), ma_p.end(), std::back_inserter(ta_ble),
                 // Lambda function extracts the second field from the map row.
                 // The std::transform() STL algorithm passes in rows (pairs) of the ma_p variable
                 // through the brackets ().
  //                [] (std::pair<int, int> const &r_ow) {return r_ow.second;}
  // );

  // Copy the map to the output vector using std::for_each() STL algorithm.
  // The STL algorithm std::for_each() is similar to a for loop.
  // std::for_each(ma_p.begin(), ma_p.end(),
                // Lambda function extracts the second field from the map row.
                // The std::for_each() STL algorithm passes in rows (pairs) of the ma_p variable 
                // through the brackets ().
                // The brackets [] are used to pass in variables from the outer scope of the lambda function.
                // The "&" passes the outer scope variables by reference.
  //               [&ta_ble](std::pair<int, int>  r_ow) {
  //                 ta_ble.push_back(r_ow.second);
  //               }
  // );
  
  return ta_ble;

}  // end calc_table




////////////////////////////
// STL algorithms


// The function sort_num() sorts the elements of a numeric vector.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<double> sort_num(std::vector<double> in_puts) {
  
  std::sort(in_puts.begin(), in_puts.end());
  return in_puts;
  
}  // end sort_num

// The function sort_string() sorts the elements of a vector of strings.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<std::string> sort_string(std::vector<std::string> in_puts) {
  
  std::sort(in_puts.begin(), in_puts.end());
  return in_puts;
  
}  // end sort_string


// Define a comparison function.
// [[Rcpp::export]]
bool is_greater(const std::vector<double> in_puts, int i1, int i2) {
  return (in_puts[i1] < in_puts[i2]);
}  // end is_greater


// Define a comparison functor as a class.
class is_greater_functor {
  
public:
  // Overloaded operator - actual function
  bool operator()(std::vector<double> in_puts, int i1, int i2) { return (in_puts[i1] < in_puts[i2]); }
  
};  // end is_greater_functor


// Define a comparison functor as a class.
class is_greater_functor2 {
  
  std::vector<double> in_puts;
  
public:
  // Constructor
  is_greater_functor2(std::vector<double> in_puts) : in_puts(in_puts) {}
  
  // Overloaded operator - actual function
  bool operator()(int i1, int i2) { return (in_puts[i1] < in_puts[i2]); }
  
};  // end is_greater_functor2



// The function sort_index() calculates the sort index of a numeric vector.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<int> sort_index(const std::vector<double> in_puts) {
  
  int num_el = in_puts.size();
  // Define index of integers along in_puts
  std::vector<int> in_dex(num_el);
  // Fill the vector in_dex with a sequence of consecutive integers.
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(in_dex.begin(), in_dex.end(), 0);
  // Sort the index according to the order of in_puts
  // is_greater_functor is_greater2;
  sort(in_dex.begin(), in_dex.end(), 
       // Lambda function defines sorting order.
       // The brackets [] are used to pass in variables from the outer scope of the lambda function.
       // The "&" passes the outer scope variables by reference.
       [&in_puts](int i1, int i2) {return in_puts[i1] < in_puts[i2];}
       // Or call functor is_greater_functor2() - elegant but very slow!!!
       // is_greater_functor2(in_puts)
       // Or call function is_greater() using std::bind() - very slow!!!
       // std::bind(is_greater, in_puts, std::placeholders::_1, std::placeholders::_2)
       // Call function is_greater() using lambda function - also very slow!!!
       // [&in_puts](int i1, int i2) {return is_greater(in_puts, i1, i2);}
       // Or call function is_greater() using std::bind() - very slow!!!
       // std::bind(is_greater2, in_puts, std::placeholders::_1, std::placeholders::_2)
  );
  
  return in_dex;
  
}  // end sort_index



// The function calc_ranks() calculates the ranks of the numeric vector elements.
// It uses the STL algorithm std::sort().
// The ranks of the elements are equal to the reverse permutation index.
// The reverse permutation index is calculated in two steps, applying the STL algorithm std::sort() twice.
// First, the permutation index is calculated, and then the ranks are calculated by sorting the sequence 
// of consecutive integers according to the order of the permutation index.
// [[Rcpp::export]]
std::vector<int> calc_ranks(const std::vector<double> &in_puts) {
  
  size_t num_el = in_puts.size();
  // size_t num_el = sizeof(in_puts);
  
  // Define index of integers along in_puts
  std::vector<int> in_dex(num_el);
  // Define the ranks of the vector elements
  std::vector<int> rank_s(num_el);
  // Fill the vectors with a sequence of consecutive integers.
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(in_dex.begin(), in_dex.end(), 0);
  std::iota(rank_s.begin(), rank_s.end(), 0);
  
  // Calculate the permutation index by sorting the sequence according to the order of in_puts
  std::sort(in_dex.begin(), in_dex.end(), 
       // Lambda function defines sorting order.
       // The brackets [] are used to pass in variables from the outer scope of the lambda function.
       // The "&" passes the outer scope variables by reference.
       [&in_puts](int i1, int i2) {return in_puts[i1] < in_puts[i2];});
  
  // Calculate the ranks by sorting the sequence according to the order of the permutation index
  std::sort(rank_s.begin(), rank_s.end(), 
       // Lambda function defines sorting order.
       // The brackets [] are used to pass in variables from the outer scope of the lambda function.
       // The "&" passes the outer scope variables by reference.
       [&in_dex](int i1, int i2) {return in_dex[i1] < in_dex[i2];});
  
  return rank_s;
  
}  // end calc_ranks



// The function remove_dup() removes consecutive duplicate elements 
// from the input vector of strings.
// It doesn't remove all duplicate elements.  
// It doesn't remove duplicate elements which don't neighbor each other.
// It uses the STL algorithm std::unique().
// [[Rcpp::export]]
std::vector<std::string> remove_dup(std::vector<std::string> string_s) {
  
  int num_el = string_s.size();
  // Define vector iterator
  std::vector<std::string>::iterator stri_ng;
  // Define vector of output strings
  std::vector<std::string> out_put = string_s;
  out_put.reserve(num_el);
  
  // Remove consecutive duplicate elements
  stri_ng = std::unique(out_put.begin(), out_put.end());
  // Resize the output vector
  out_put.resize(std::distance(out_put.begin(), stri_ng));

  return out_put;
  
}  // end remove_dup



// The function calc_unique_sort() returns a vector with the unique
// elements of the input vector of strings.
// It uses the STL algorithms std::sort() and std::unique().
// It extracts the unique elements in two steps.
// First it sorts the elements of the input vector, and then 
// it removes consecutive duplicate elements.
// The function calc_unique_sort() is slower than the R function unique().
// [[Rcpp::export]]
std::vector<std::string> calc_unique_sort(std::vector<std::string> string_s) {
  
  int num_el = string_s.size();
  // Define vector iterator
  std::vector<std::string>::iterator stri_ng;
  // Define vector of output strings
  std::vector<std::string> out_put = string_s;
  out_put.reserve(num_el);
  // Sort the elements of the input vector
  std::sort(out_put.begin(), out_put.end());
  // Remove consecutive duplicate elements
  stri_ng = std::unique(out_put.begin(), out_put.end());
  // Resize the output vector
  out_put.resize(std::distance(out_put.begin(), stri_ng));
  
  return out_put;
  
}  // end calc_unique_sort




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
std::vector<int> map_out(std::string stri_ng = "blah") {
  
  // Define map
  std::map <std::string, int> ma_p;
  ma_p["first"] = 1;
  ma_p["second"] = 2;
  ma_p["third"] = 3;
  
  std::vector<int> vec_tor;
  vec_tor.reserve(ma_p.size());
  
  // std::cout << "Value: " << stri_ng << std::endl;
  // Copy all values from a map to vector using transform() and a function pointer
  std::transform(ma_p.begin(), ma_p.end(), std::back_inserter(vec_tor), &get_val);
  
  return vec_tor;
  
}  // end map_out



// The function print_string() prints the elements of a vector of strings 
// using the STL algorithm std::for_each().
// [[Rcpp::export]]
void print_string(std::vector<std::string> in_puts) {
  
  std::for_each(in_puts.begin(), in_puts.end(),
                // Lambda function prints a single string
                [](std::string stri_ng) { std::cout << stri_ng << ", "; });
  std::cout << std::endl;

}  // end print_string



// Example of sorting a defined class structure.
// https://thispointer.com/stl-algorithm-stdsort-tutorial-and-example/

// Define a Person class.
// The Person class can't be exported to Rcpp but can be used internally.
class Person {
public:
  // Fields
  std::string m_name;
  int m_id;
  // Constructor
  Person(std::string name, int id) : m_name(name), m_id(id) {}
  
  // Comparison operator to compare person's IDs.
  bool operator < (const Person &per_son) {
    if (m_id < per_son.m_id)
      return true;
    else
      return false;
  }  // end comparison operator
};  // end Person class


// Person comparison function.
struct compare_names {
  bool operator()(const Person &first, const Person &sec) 
  {
    if (first.m_name < sec.m_name)
      return true;
    else
      return false;
  }
};  // end compare_names


// The function sort_ids() sorts persons according to their IDs.
// The Person class can't be exported to Rcpp but can be used internally.
// [[Rcpp::export]]
void sort_ids(std::vector<int> id_s) {
  // Define a vector of persons
  std::vector<Person> person_s = { Person("aaa", id_s[0]), 
                                   Person("kkk", id_s[1]),
                                   Person("ddd", id_s[2]), 
                                   Person("abc", id_s[3]) };
  
  // Sort the vector of persons according to their IDs.
  std::sort(person_s.begin(), person_s.end());
  
  std::cout << "Sorted Persons List based on ID\n";
  std::for_each(person_s.begin(), person_s.end(), 
                [](Person &per_son) { std::cout << per_son.m_id << " :: " << per_son.m_name << std::endl; });

}  // end sort_ids


// Define a vector of persons
// std::vector<int> id_s = c(3, 2, 4, 1);

// The Person class can't be exported to Rcpp but can be used internally.
std::vector<Person> create_persons(std::vector<int> id_s) {

  std::vector<Person> person_s = { Person("aaa", id_s[0]),
                                   Person("kkk", id_s[1]),
                                   Person("ddd", id_s[2]),
                                   Person("abc", id_s[3]) };
  return person_s;

}  // end create_persons


// The function sort_names() sorts persons according to their names.
// The Person class can't be exported to Rcpp but can be used internally.
// [[Rcpp::export]]
void sort_names(std::vector<int> id_s) {
    
  std::vector<Person> person_s = { Person("aaa", id_s[0]),
                                   Person("kkk", id_s[1]),
                                   Person("ddd", id_s[2]),
                                   Person("abc", id_s[3]) };
  
  // Sort the vector of persons according to their names.
  std::sort(person_s.begin(), person_s.end(), 
            // Compare using compare_names()
            // compare_names()
            // Or pass a lambda function
            // Lambda function defines sorting order.
            [](Person &first, Person &sec) {
              if (first.m_name < sec.m_name)
                return true; 
              else 
                return false;}
  );
  // std::sort(person_s.begin(), person_s.end());
  
  std::cout << "Sorted Persons List based on name\n";
  std::for_each(person_s.begin(), person_s.end(), 
                [](Person &per_son) { std::cout << per_son.m_id << " :: " << per_son.m_name << std::endl; });
  
}  // end sort_names



// square_it() is a non-exported function which squares a double.
// It can be used by other functions.
double square_it(double x) { return x*x; }

// The STL functional std::transform() can be used to apply a function over a vector.
// The function square_vec() squares the elements of a numeric vector 
// by calling a lambda function using the functional std::transform().
// Lambda functions are anonymous functions which can be passed to functionals.
// [[Rcpp::export]]
std::vector<double> square_vec(const std::vector<double> vec_tor) {
  std::vector<double> out_put(vec_tor.size());
  
  std::transform(vec_tor.begin(), vec_tor.end(), out_put.begin(), 
                 // Pass function square_it() to functional std::transform()
                 // square_it
                 // Or pass a lambda function
                 [](double x) { return x*x; }
  ); // end std::transform()
  
  return out_put;
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

// [[Rcpp::export]]
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


// The function match_it() reproduces the R function findInterval().
// The function match_it() matches its inputs with a vector of break points.
// in_puts is a vector of inputs to be matched with the break points.
// break_s is a vector of break points.
// The matches are the indices of the break points closest to the inputs.
// [[Rcpp::export]]
std::vector<int> match_it(std::vector<double> in_puts, std::vector<double> break_s) {
  
  // Allocate vector of match_es: the break points that match the in_puts.
  std::vector<int> match_es(in_puts.size());
  // Allocate iterators (pointers) to the in_puts and match_es
  std::vector<int>::iterator ma_tch;
  std::vector<double>::iterator in_put, brea_k;

  // Loop over the vectors of in_puts and calculate the match_es
  for (in_put = in_puts.begin(), ma_tch = match_es.begin(); in_put != in_puts.end(); ++in_put, ++ma_tch) {
    // Find closest break point to the in_put
    brea_k = std::upper_bound(break_s.begin(), break_s.end(), *in_put);
    // Calculate the index of the closest break point
    *ma_tch = std::distance(break_s.begin(), brea_k);
  }  // end for
  
  return match_es;
  
}  // end match_it



// The function read_back() returns the reverse of its 
// input vector using a reverse_iterator.
// [[Rcpp::export]]
std::vector<int> read_back(std::vector<int> in_puts) {
  
  // Define vectors
  std::vector<int> re_verse;

  // Initialise first value
  std::vector<int>::reverse_iterator in_put;
  
  for(in_put = in_puts.rbegin(); in_put != in_puts.rend(); ++in_put) {
    re_verse.push_back(*in_put);
  }
  
  return re_verse;
  
}  // end read_back




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


