////////////////////////////
// Utility functions using the Standard Template Library (STL)
////////////////////////////


double add2(double x1, double x2);
double mult2(double x1, double x2);
std::vector<double> mult_vec_lambda(double factorv, std::vector<double> vectorv);
std::vector<double> run_mean(double lambda, std::vector<double> vectorv);
std::vector<double> run_int(double lambda, std::vector<int> vectorv);


////////////////////////////////////////////////////////////
// Functions templates
////////////////////////////////////////////////////////////

// Define template function - accepts template arguments and returns template
// Calculate the maximum of two numbers
template <class numericc>
numericc calc_max(numericc a, numericc b) {
  numericc result;
  result = (a>b)? a : b;
  return (result);
}

// Define template functional - accepts a function template as argument
// The functional run_fun3() accepts a function of two variables of any type, executes it, and returns a std::vector.
template <typename FuncType, class VarType1, class VarType2>
std::vector<double> run_fun3(FuncType func, VarType1 factorv, VarType2 vectorv) {
  return (func)(factorv, vectorv);
}  // end run_fun3
