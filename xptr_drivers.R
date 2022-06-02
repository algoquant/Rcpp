########################
### Driver scripts for RcppXPtrUtils functions
### Package RcppXPtrUtils creates function pointers
########################

library(Rcpp)
library(RcppXPtrUtils)

# Compile Rcpp functions
Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_xptr.cpp")

# Create function pointer
func_cpp <- RcppXPtrUtils::cppXPtr(
  "double foo(double lambda, double x) { return lambda*x; }"
)

# Run function that accepts the function pointer
run_cpp(func_cpp, 5, 1:4)

run_cpp(func_cpp, 0.5, 1:4)



