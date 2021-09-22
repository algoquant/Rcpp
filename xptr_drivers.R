########################
### Driver scripts for RcppXPtrUtils functions
### Package RcppXPtrUtils creates function pointers
########################

library(Rcpp)
library(RcppXPtrUtils)

# Compile Rcpp functions
Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_xptr.cpp")

# Create function pointer
func_cpp <- RcppXPtrUtils::cppXPtr(
  "double foo(double fac_tor, double x) { return fac_tor*x; }"
)

# Run function that accepts the function pointer
run_cpp(func_cpp, 5, 1:4)

run_cpp(func_cpp, 0.5, 1:4)



