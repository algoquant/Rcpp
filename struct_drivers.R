########################
### Driver scripts for testing C++ Structures
########################

library(microbenchmark)
# Compile Rcpp functions
# Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_struct.cpp")

datav <- sample(1001)
all.equal(sort_index(datav)+1, order(datav))
summary(microbenchmark(
  rcode=order(datav),
  stl=sort_index(datav),
  times=10))[, c(1, 4, 5)]


datav <- round(runif(1e2), 2)
all.equal(drop(HighFreq::calc_ranks(datav)), calc_ranks(datav)+1)
summary(microbenchmark(
  rcode=rank(datav),
  arma=HighFreq::calc_ranks(datav),
  stl=calc_ranks(datav),
  times=10))[, c(1, 4, 5)]


datav <- letters[sample(7)]
print_string(sort_string(datav))

datav <- sample(3, 21, replace=TRUE)
encode_it(datav)
calc_table(datav)
all.equal(calc_table(datav), as.numeric(table(datav)))
summary(microbenchmark(
  rcode=table(datav),
  stl=calc_table(datav),
  times=10))[, c(1, 4, 5)]

datav <- runif(30, 1, 100)
breakv <- runif(10, 1, 100)
breakv <- sort(breakv)
all.equal(match_it(datav, breakv), findInterval(datav, breakv))
summary(microbenchmark(
  rcode=findInterval(datav, breakv),
  stl=match_it(datav, breakv),
  times=10))[, c(1, 4, 5)]


datav <- sample(11, 1e3, replace=TRUE)
all.equal(datav[not_dup(datav)], unique(datav))
summary(microbenchmark(
  rcode=unique(datav),
  stl=not_dup(datav),
  times=10))[, c(1, 4, 5)]


all.equal(sort(calc_unique(datav)), sort(unique(datav)))
summary(microbenchmark(
  rcode=unique(datav),
  stl_loop=calc_unique_loop(datav),
  stl=calc_unique(datav),
  times=10))[, c(1, 4, 5)]


sort_ids(c(3, 2, 4, 1))
sort_names(c(3, 2, 4, 1))

stringv <- sample(letters, 1e3, replace=TRUE)

all.equal(sort(calc_unique_str(stringv)), sort(unique(stringv)))
summary(microbenchmark(
  rcode=unique(stringv),
  stl=calc_unique_str(stringv),
  times=10))[, c(1, 4, 5)]

all.equal(calc_unique_sort(stringv), sort(unique(stringv)))
summary(microbenchmark(
  rcode=unique(stringv),
  stl=calc_unique_sort(stringv),
  times=10))[, c(1, 4, 5)]


