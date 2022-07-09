########################
### Scripts for running Rcpp code for Rcpp::XPtr function pointers.
########################

# Compile Rcpp functions
Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_xptr.cpp")

# Run a function of a vector based on a string
drop(runfunvec("fun1", 11:5))
drop(runfunvec("fun2", 11:5))
# Produces error because no such function
drop(runfunvec("fun222", 11:5))


# Run a function of a matrix based on a string
datav <- matrix(21:5)
drop(runfunmat("run_mean", datav, 0.9))
drop(runfunmat("run_max", datav, 0.9))


# Pass an R function into Rcpp using an Rcpp Function object
execute_r(rnorm, 3)
execute_r(seq, 2)


# Code for RcppXPtrUtils functions.
# Package RcppXPtrUtils creates function pointers
# https://gallery.rcpp.org/articles/passing-cpp-function-pointers-rcppxptrutils/

# This works:
# Create function pointer from C++ code
rexpfun <- RcppXPtrUtils::cppXPtr(
  "SEXP foo(int n, double rate) { return rexp(n, rate); }"
)
# Check if function pointer is valid
RcppXPtrUtils::checkXPtr(rexpfun, "SEXP", c("int", "double"))
# Run function that accepts the function pointer
execute_cpp(rexpfun, 5)

# But this doesn't work: dexp() instead of rexp()
# Create function pointer from C++ code
doublefun <- RcppXPtrUtils::cppXPtr(
  "SEXP foo(double x, double rate) { return dexp(x, rate); }"
)

# Can't run execute_cpp() for anything but rexp()
# Or any modified version of execute_cpp() to account for different arguments


# This works:
# Create function pointer from C++ code
addfun <- RcppXPtrUtils::cppXPtr(
  "double foo(double x1, double x2) { return (x1+x2); }"
)
# Check if function pointer is valid
RcppXPtrUtils::checkXPtr(addfun, "double", c("double", "double"))

# Create function pointer from C++ code
multfun <- RcppXPtrUtils::cppXPtr(
  "double foo(double lambda, double x) { return lambda*x; }"
)

# Run function that accepts the function pointer
run_cpp(addfun, 5, 1:4)
run_cpp(multfun, 5, 1:4)
# This doesn't work:
run_cpp(mult2, 5, 1:4)



