########################
### C++ STL driver scripts for testing functions
########################

library(microbenchmark)
# Compile Rcpp functions
Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_stl.cpp")


for (i in 1:5) count_er()

da_ta <- sample(1001)
all.equal(sort_index(da_ta)+1, order(da_ta))
summary(microbenchmark(
  rcode=order(da_ta),
  stl=sort_index(da_ta),
  times=10))[, c(1, 4, 5)]


da_ta <- round(runif(1e2), 2)
all.equal(drop(HighFreq::calc_ranks(da_ta)), calc_ranks(da_ta)+1)
summary(microbenchmark(
  rcode=rank(da_ta),
  arma=HighFreq::calc_ranks(da_ta),
  stl=calc_ranks(da_ta),
  times=10))[, c(1, 4, 5)]


da_ta <- letters[sample(7)]
print_string(sort_string(da_ta))

da_ta <- sample(3, 21, replace=TRUE)
encode_it(da_ta)
calc_table(da_ta)
all.equal(calc_table(da_ta), as.numeric(table(da_ta)))
summary(microbenchmark(
  rcode=table(da_ta),
  stl=calc_table(da_ta),
  times=10))[, c(1, 4, 5)]

da_ta <- runif(30, 1, 100)
break_s <- runif(10, 1, 100)
break_s <- sort(break_s)
all.equal(match_it(da_ta, break_s), findInterval(da_ta, break_s))
summary(microbenchmark(
  rcode=findInterval(da_ta, break_s),
  stl=match_it(da_ta, break_s),
  times=10))[, c(1, 4, 5)]


da_ta <- sample(11, 1e3, replace=TRUE)
all.equal(da_ta[not_dup(da_ta)], unique(da_ta))
summary(microbenchmark(
  rcode=unique(da_ta),
  stl=not_dup(da_ta),
  times=10))[, c(1, 4, 5)]


all.equal(sort(calc_unique(da_ta)), sort(unique(da_ta)))
summary(microbenchmark(
  rcode=unique(da_ta),
  stl_loop=calc_unique_loop(da_ta),
  stl=calc_unique(da_ta),
  times=10))[, c(1, 4, 5)]


sort_ids(c(3, 2, 4, 1))
sort_names(c(3, 2, 4, 1))

string_s <- sample(letters, 1e3, replace=TRUE)

all.equal(sort(calc_unique_str(string_s)), sort(unique(string_s)))
summary(microbenchmark(
  rcode=unique(string_s),
  stl=calc_unique_str(string_s),
  times=10))[, c(1, 4, 5)]

all.equal(calc_unique_sort(string_s), sort(unique(string_s)))
summary(microbenchmark(
  rcode=unique(string_s),
  stl=calc_unique_sort(string_s),
  times=10))[, c(1, 4, 5)]


