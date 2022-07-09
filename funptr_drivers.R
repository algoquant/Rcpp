########################
### Scripts for running Rcpp code for Rcpp::XPtr function pointers.
########################

# Compile Rcpp functions
Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_funptr.cpp")

# Run a function of two double numbers based on a string
# Using a functional without using a function pointer
runfun2("add2", 11.5, 2)
runfun2("mult2", 11.5, 2)
# Produces error because no such function
runfun2("fun222", 11.5, 2)

# Run a function of two double numbers based on a string
# Using a function pointer instead of a functional
runfun("add2", 11.5, 2)
runfun("mult2", 11.5, 2)


# Run a function of a vector based on a string
# Using a function pointer instead of a functional
drop(runfunvec("fun1", 11:5))
drop(runfunvec("fun2", 11:5))
# Produces error because no such function
drop(runfunvec("fun222", 11:5))

