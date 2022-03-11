////////////////////////////
// Functions to test C++ syntax using the Standard Template Library (STL)
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_stl.cpp")

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
// It may be simpler to use std::accumulate().
// [[Rcpp::export]]
int sum_it(std::vector<int> vectorv) {
  
  // Define STL iterator
  // std::vector<int>::iterator it_er;
  
  // Sum up elements of a vector using STL iterator
  int total = 0;
  for (auto it_er: vectorv) {
    // Old-style loop
    // for (auto it_er = vectorv.begin(); it_er != vectorv.end(); ++it_er) {
    total += it_er;
  }  // end for
  
  return total;
}  // end sum_it



// The function select_it() selects elements of a vector using an STL iterator.
// [[Rcpp::export]]
void select_it(std::vector<double> vectorv, int shiftv) {
  
  // Define STL iterator
  std::vector<double>::iterator it_er = vectorv.begin();
  // Select element from the front of the vector
  it_er = next(it_er, shiftv-1);
  
  // Print iterator value
  std::cout << "The " << shiftv << "th ";
  std::cout << "vector element from the front is = "; 
  std::cout << *it_er << " "; 
  std::cout << std::endl;

  // Select element from the end of the vector
  it_er = vectorv.end();
  it_er = prev(it_er, shiftv);

  // Print iterator value
  std::cout << "The " << shiftv << "th ";
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
List encode_it(std::vector<double> inputs) {
  
  // Define vector of input data
  std::vector<double> datav;
  // Define vector of data counts (repeats)
  std::vector<int> count_s;
  
  // Define iterators
  std::vector<double>::iterator input = inputs.begin();
  std::vector<int>::reverse_iterator cou_nt = count_s.rbegin();
  
  // Initialise the data
  double lastv = *input;
  datav.push_back(lastv);
  count_s.push_back(1);
  
  // Perform loop over inputs using the iterator
  for (input = inputs.begin() + 1; input != inputs.end(); ++input) {
    if (lastv == *input) {
      (*cou_nt)++;
    } else {
      datav.push_back(*input);
      count_s.push_back(1);
      cou_nt = count_s.rbegin();
      lastv = *input;
    }  // end if
  }  // end for
  
  return List::create(
    _["data"] = datav,
    _["counts"] = count_s
  );
}  // end encode_it




// wippp
// The function decode_it() decodes a vector from its run length encoding.
// [[Rcpp::export]]
std::vector<double> decode_it(List inputs) {
  
  // Define vector of input data
  std::vector<double> datav;
  // Define vector of data counts (repeats)
  // std::vector<int> count_s;
  
  // Define iterators
  // std::vector<double>::iterator input = inputs.begin();
  // std::vector<int>::reverse_iterator cou_nt = count_s.rbegin();
  // 
  // // Initialise the data
  // double lastv = *input;
  // datav.push_back(lastv);
  // count_s.push_back(1);
  // 
  // // Perform loop over inputs using the iterator
  // for (input = inputs.begin() + 1; input != inputs.end(); ++input) {
  //   if (lastv == *input) {
  //     (*cou_nt)++;
  //   } else {
  //     datav.push_back(*input);
  //     count_s.push_back(1);
  //     cou_nt = count_s.rbegin();
  //     lastv = *input;
  //   }  // end if
  // }  // end for
  
  return datav;
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
std::vector<bool> not_dup(std::vector<int> vectorv) {
  
  // int num_el = vectorv.size();
  // Define iterator for input vector
  // std::vector<int>::iterator ele_ment;
  // Define Boolean output vector
  std::vector<bool> output(vectorv.size());
  // Define Boolean iterator for output vector
  std::vector<bool>::iterator is_dup = output.begin();
  // Define unordered_set with unique elements
  std::unordered_set<int> se_t;
  
  // Copy the elements of the input vector into a set
  for (auto ele_ment: vectorv) {
    // for (auto ele_ment = vectorv.begin(); ele_ment != vectorv.end(); ++ele_ment) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    *is_dup++ = se_t.insert(ele_ment).second;
  }  // end for
  
  return output;
  
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
std::vector<double> calc_unique(std::vector<double> vectorv) {
  
  // Define unordered_set and copy the input vector into it.
  // The set contains only unique elements.
  std::unordered_set<double> se_t(vectorv.begin(), vectorv.end());
  // Define output vector and copy the set into it.
  std::vector<double> output(se_t.begin(), se_t.end());
  
  return output;
  
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
std::vector<double> calc_unique_loop(std::vector<double> vectorv) {
  
  // Define iterator for input vector
  // std::vector<double>::iterator ele_ment;
  // Define output vector
  std::vector<double> output;
  // Define unordered_set with unique elements
  std::unordered_set<double> se_t;
  
  // Copy the elements of the input vector into a set
  for (auto ele_ment: vectorv) {
    // for (auto ele_ment = vectorv.begin(); ele_ment != vectorv.end(); ++ele_ment) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    if (se_t.insert(ele_ment).second)
      output.push_back(ele_ment);
  }  // end for
  
  return output;

}  // end calc_unique_loop



// The function calc_unique_int() returns a vector with the unique
// elements of the integer input vector.
// It uses the STL set structure.
// The method .insert() adds an element to a set.
// The .insert().second is a Boolean equal to TRUE if the element is new, 
// i.e it's not already in the set.
// The function calc_unique() is about as fast as the R function unique().
// [[Rcpp::export]]
std::vector<int> calc_unique_int(std::vector<int> vectorv) {
  
  // Define iterator for input vector
  // std::vector<int>::iterator ele_ment;
  // Define output vector
  std::vector<int> output;
  // Define unordered_set with unique elements
  std::unordered_set<int> se_t;
  
  // Copy the elements of the input vector into a set
  for (auto ele_ment: vectorv) {
    // for (auto ele_ment = vectorv.begin(); ele_ment != vectorv.end(); ++ele_ment) {
    // The .insert().second is a Boolean equal to TRUE if the element 
    // is new, i.e it's not already in the set.
    if (se_t.insert(ele_ment).second)
      output.push_back(ele_ment);
  }  // end for
  
  return output;
  
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
  std::vector<std::string> output;
  output.reserve(num_el);
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
    output.push_back(uniqu_e);
  }  // end for
  
  return output;

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
std::vector<int> calc_table(std::vector<int> vectorv) {

  // Define iterator for input vector
  // std::vector<int>::iterator ele_ment;
  // Define map iterator
  // std::map <int, int>::iterator cou_nt;
  // Define map
  std::map <int, int> ma_p;
  
  // Copy the elements of the input vector into a map
  for (auto ele_ment: vectorv) {
    // for (auto ele_ment = vectorv.begin(); ele_ment != vectorv.end(); ++ele_ment) {
    ma_p[ele_ment]++;
  }  // end for
  
  // Define contingency table
  std::vector<int> tablev;
  tablev.reserve(ma_p.size());
  
  // Below are four different methods of copying the map to the output vector
  
  // Explicit for loop:
  // https://thispointer.com/how-to-copy-all-values-from-a-map-to-a-vector-in-c/
  for (auto cou_nt: ma_p) {
    tablev.push_back(cou_nt.second);
  }  // end for
  
  // Old style loop:
  // for (cou_nt=ma_p.begin(); cou_nt!=ma_p.end(); ++cou_nt) {
  //   tablev.push_back(cou_nt->second);
  // }  // end for
  
  // Copy the map to the output vector using std::transform() STL algorithm.
  // The STL algorithm std::transform() is similar to the R functional apply().
  // std::transform(ma_p.begin(), ma_p.end(), std::back_inserter(tablev),
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
  //               [&tablev](std::pair<int, int>  r_ow) {
  //                 tablev.push_back(r_ow.second);
  //               }
  // );
  
  return tablev;

}  // end calc_table




////////////////////////////
// STL algorithms


// The function sort_num() sorts the elements of a numeric vector.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<double> sort_num(std::vector<double> inputs) {
  
  std::sort(inputs.begin(), inputs.end());
  return inputs;
  
}  // end sort_num

// The function sort_string() sorts the elements of a vector of strings.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<std::string> sort_string(std::vector<std::string> inputs) {
  
  std::sort(inputs.begin(), inputs.end());
  return inputs;
  
}  // end sort_string


// Define a comparison function.
// [[Rcpp::export]]
bool is_greater(const std::vector<double> inputs, int i1, int i2) {
  return (inputs[i1] < inputs[i2]);
}  // end is_greater


// Define a comparison functor as a class.
class is_greater_functor {
  
public:
  // Overloaded operator - actual function
  bool operator()(std::vector<double> inputs, int i1, int i2) { return (inputs[i1] < inputs[i2]); }
  
};  // end is_greater_functor


// Define a comparison functor as a class.
class is_greater_functor2 {
  
  std::vector<double> inputs;
  
public:
  // Constructor
  is_greater_functor2(std::vector<double> inputs) : inputs(inputs) {}
  
  // Overloaded operator - actual function
  bool operator()(int i1, int i2) { return (inputs[i1] < inputs[i2]); }
  
};  // end is_greater_functor2



// The function sort_index() calculates the sort index of a numeric vector.
// It uses the STL algorithm std::sort().
// [[Rcpp::export]]
std::vector<int> sort_index(const std::vector<double> inputs) {
  
  int num_el = inputs.size();
  // Define index of integers along inputs
  std::vector<int> indeks(num_el);
  // Fill the vector indeks with a sequence of consecutive integers.
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(indeks.begin(), indeks.end(), 0);
  // Sort the index according to the order of inputs
  // is_greater_functor is_greater2;
  sort(indeks.begin(), indeks.end(), 
       // Lambda function defines sorting order.
       // The brackets [] are used to pass in variables from the outer scope of the lambda function.
       // The "&" passes the outer scope variables by reference.
       [&inputs](int i1, int i2) {return inputs[i1] < inputs[i2];}
       // Or call functor is_greater_functor2() - elegant but very slow!!!
       // is_greater_functor2(inputs)
       // Or call function is_greater() using std::bind() - very slow!!!
       // std::bind(is_greater, inputs, std::placeholders::1, std::placeholders::2)
       // Call function is_greater() using lambda function - also very slow!!!
       // [&inputs](int i1, int i2) {return is_greater(inputs, i1, i2);}
       // Or call function is_greater() using std::bind() - very slow!!!
       // std::bind(is_greater2, inputs, std::placeholders::1, std::placeholders::2)
  );
  
  return indeks;
  
}  // end sort_index



// The function calc_ranks() calculates the ranks of the numeric vector elements.
// It uses the STL algorithm std::sort().
// The ranks of the elements are equal to the reverse permutation index.
// The reverse permutation index is calculated in two steps, applying the STL algorithm std::sort() twice.
// First, the permutation index is calculated, and then the ranks are calculated by sorting the sequence 
// of consecutive integers according to the order of the permutation index.
// [[Rcpp::export]]
std::vector<int> calc_ranks(const std::vector<double> &inputs) {
  
  size_t num_el = inputs.size();
  // size_t num_el = sizeof(inputs);
  
  // Define index of integers along inputs
  std::vector<int> indeks(num_el);
  // Define the ranks of the vector elements
  std::vector<int> ranks(num_el);
  // Fill the vectors with a sequence of consecutive integers.
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(indeks.begin(), indeks.end(), 0);
  std::iota(ranks.begin(), ranks.end(), 0);
  
  // Calculate the permutation index by sorting the sequence according to the order of inputs
  std::sort(indeks.begin(), indeks.end(), 
       // Lambda function defines sorting order.
       // The brackets [] are used to pass in variables from the outer scope of the lambda function.
       // The "&" passes the outer scope variables by reference.
       [&inputs](int i1, int i2) {return inputs[i1] < inputs[i2];});
  
  // Calculate the ranks by sorting the sequence according to the order of the permutation index
  std::sort(ranks.begin(), ranks.end(), 
       // Lambda function defines sorting order.
       // The brackets [] are used to pass in variables from the outer scope of the lambda function.
       // The "&" passes the outer scope variables by reference.
       [&indeks](int i1, int i2) {return indeks[i1] < indeks[i2];});
  
  return ranks;
  
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
  std::vector<std::string> output = string_s;
  output.reserve(num_el);
  
  // Remove consecutive duplicate elements
  stri_ng = std::unique(output.begin(), output.end());
  // Resize the output vector
  output.resize(std::distance(output.begin(), stri_ng));

  return output;
  
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
  std::vector<std::string> output = string_s;
  output.reserve(num_el);
  // Sort the elements of the input vector
  std::sort(output.begin(), output.end());
  // Remove consecutive duplicate elements
  stri_ng = std::unique(output.begin(), output.end());
  // Resize the output vector
  output.resize(std::distance(output.begin(), stri_ng));
  
  return output;
  
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
  
  std::vector<int> vectorv;
  vectorv.reserve(ma_p.size());
  
  // std::cout << "Value: " << stri_ng << std::endl;
  // Copy all values from a map to vector using transform() and a function pointer
  std::transform(ma_p.begin(), ma_p.end(), std::back_inserter(vectorv), &get_val);
  
  return vectorv;
  
}  // end map_out



// The function print_string() prints the elements of a vector of strings 
// using the STL algorithm std::for_each().
// [[Rcpp::export]]
void print_string(std::vector<std::string> inputs) {
  
  std::for_each(inputs.begin(), inputs.end(),
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
std::vector<double> square_vec(const std::vector<double> vectorv) {
  std::vector<double> output(vectorv.size());
  
  std::transform(vectorv.begin(), vectorv.end(), output.begin(), 
                 // Pass function square_it() to functional std::transform()
                 // square_it
                 // Or pass a lambda function
                 [](double x) { return x*x; }
  ); // end std::transform()
  
  return output;
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
// inputs is a vector of inputs to be matched with the break points.
// break_s is a vector of break points.
// The matches are the indices of the break points closest to the inputs.
// [[Rcpp::export]]
std::vector<int> match_it(std::vector<double> inputs, std::vector<double> break_s) {
  
  // Allocate vector of match_es: the break points that match the inputs.
  std::vector<int> match_es(inputs.size());
  // Allocate iterators (pointers) to the inputs and match_es
  std::vector<int>::iterator ma_tch;
  std::vector<double>::iterator input, brea_k;

  // Loop over the vectors of inputs and calculate the match_es
  for (input = inputs.begin(), ma_tch = match_es.begin(); input != inputs.end(); ++input, ++ma_tch) {
    // Find closest break point to the input
    brea_k = std::upper_bound(break_s.begin(), break_s.end(), *input);
    // Calculate the index of the closest break point
    *ma_tch = std::distance(break_s.begin(), brea_k);
  }  // end for
  
  return match_es;
  
}  // end match_it



// The function read_back() returns the reverse of its 
// input vector using a reverse_iterator.
// [[Rcpp::export]]
std::vector<int> read_back(std::vector<int> inputs) {
  
  // Define vectors
  std::vector<int> re_verse;

  // Initialise first value
  std::vector<int>::reverse_iterator input;
  
  for(input = inputs.rbegin(); input != inputs.rend(); ++input) {
    re_verse.push_back(*input);
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


