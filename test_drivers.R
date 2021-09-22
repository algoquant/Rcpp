########################
### RcppArmadillo driver scripts
########################


#########
### Scripts for ad-hoc tests

library(HighFreq)

# Compile Rcpp functions
Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/test_temp.cpp")



#########
### Scripts for calling RcppArmadillo functions


#########
### Scripts for testing HighFreq::diff_vec() and HighFreq::diff_it()

re_turns <- matrix(rnorm(1e3), nc=1)
lagg <- 3
re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "IEF")])
re_turns <- na.omit(rutils::etf_env$re_turns$VTI)

foo <- lag_it(re_turns, lagg=3, pad_zeros=FALSE)
bar <- lag_vec(re_turns, lagg=3, pad_zeros=FALSE)
all.equal(drop(foo), drop(bar), check.attributes=FALSE)


look_back <- 22
stu_b <- 21
calc_var(re_turns)
calc_var(re_turns, ste_p=look_back)

# Calculate rolling sums at each point
c_sum <- roll_sum(re_turns, look_back=look_back)
r_sum <- rutils::roll_sum(re_turns, look_back=look_back)
all.equal(c_sum, coredata(r_sum), check.attributes=FALSE)
r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
lag_sum <- rbind(matrix(numeric(2*look_back), nc=2), r_sum[1:(NROW(r_sum) - look_back), ])
r_sum <- (r_sum - lag_sum)
all.equal(c_sum, r_sum, check.attributes=FALSE)

# Calculate rolling sums at end points
c_sum <- roll_sum(re_turns, look_back=look_back, stu_b=stu_b)
end_p <- (stu_b + look_back*(0:(NROW(re_turns) %/% look_back)))
end_p <- end_p[end_p < NROW(re_turns)]
r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
r_sum <- r_sum[end_p+1, ]
# r_sum <- diff_it(r_sum)
lag_sum <- rbind(numeric(2), r_sum[1:(NROW(r_sum) - 1), ])
r_sum <- (r_sum - lag_sum)
all.equal(c_sum, r_sum, check.attributes=FALSE)

# Calculate rolling sums at end points - pass in end_points
c_sum <- roll_sum(re_turns, end_points=end_p)
all.equal(c_sum, r_sum, check.attributes=FALSE)


# Create exponentially decaying weights
weight_s <- exp(-0.2*(1:11))
weight_s <- matrix(weight_s/sum(weight_s), nc=1)
# Calculate rolling weighted sum
# weight_ed <- HighFreq::roll_conv(re_turns, weight_s)
c_sum <- roll_sum(re_turns, weight_s=weight_s)
# Calculate rolling weighted sum using filter()
filter_ed <- filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
# Compare both methods
all.equal(c_sum[-(1:11), ], filter_ed[-(1:11), ], check.attributes=FALSE)

# Calculate rolling weighted sums at end points
c_sum <- roll_sum(re_turns, end_points=end_p, weight_s=weight_s)
all.equal(c_sum, filter_ed[end_p+1, ], check.attributes=FALSE)


library(microbenchmark)
summary(microbenchmark(
  rcpp=roll_sum(re_turns, look_back=look_back),
  rcode=HighFreq::roll_sum(re_turns, look_back=look_back),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary
summary(microbenchmark(
  rcpp=roll_sum(re_turns, look_back=22, stu_b=21),
  rcode=HighFreq::roll_sum(re_turns, look_back=look_back),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# Loop over stu_b
foo <- lapply(0:21, roll_rets, re_turns=re_turns, look_back=22)
bar <- sapply(foo, sd)
hist(bar)


# Compare HighFreq::diff_vec() with rutils::diff_it()
all.equal(drop(HighFreq::diff_vec(re_turns, lagg=lagg, padd=TRUE)),
          rutils::diff_it(re_turns, lagg=lagg))


# Compare HighFreq::diff_it() with rutils::diff_it()
re_turns <- matrix(re_turns, nc=5)
all.equal(HighFreq::diff_it(re_turns, lagg=lagg, padd=TRUE),
          rutils::diff_it(re_turns, lagg=lagg))



#########
### Scripts for testing HighFreq::agg_ohlc() and HighFreq::roll_ohlc()

oh_lc <- coredata(rutils::etf_env$VTI[, c(1, 5)])
all.equal(drop(HighFreq::agg_ohlc(oh_lc)),
          c(oh_lc[1, 1], max(oh_lc[, 1]), min(oh_lc[, 1]), oh_lc[NROW(oh_lc), 1], sum(oh_lc[, 2])), 
          check.attributes=FALSE)

oh_lc <- coredata(rutils::etf_env$VTI[, 1:5])
all.equal(drop(HighFreq::agg_ohlc(oh_lc)),
          c(oh_lc[1, 1], max(oh_lc[, 2]), min(oh_lc[, 3]), oh_lc[NROW(oh_lc), 4], sum(oh_lc[, 5])), 
          check.attributes=FALSE)

oh_lc <- rutils::etf_env$VTI[, 1:5]
end_points <- rutils::calc_endpoints(oh_lc, inter_val=25)
foo <- HighFreq::roll_ohlc(oh_lc, end_points-1)
# bar <- rutils::to_period(oh_lc[, 1:5], end_points=end_points)
bar <- .Call("toPeriod", oh_lc, as.integer(end_points), TRUE, NCOL(oh_lc), FALSE, FALSE, colnames(oh_lc), PACKAGE="xts")
all.equal(foo, coredata(bar), check.attributes=FALSE)

library(microbenchmark)
summary(microbenchmark(
  rcpp=HighFreq::roll_ohlc(oh_lc, end_points-1),
  rcode=rutils::to_period(oh_lc, end_points=end_points),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary



#########
### Calculate the order index, to resort back to unsorted order
un_sorted <- round(runif(10), 2)
sor_ted <- sort(un_sorted)
# Calculate the index to sort
sor_t <- (1 + drop(order_index(un_sorted)))
all.equal(un_sorted[sor_t], sor_ted)
# Calculate the index to sort
un_sort <- (1 + drop(order_index(order_index(un_sorted))))
all.equal(sor_ted[un_sort], un_sorted)

# Calculate the ranks of the elements
rank_s <- seq_along(un_sorted)[sor_t][sor_t]
all.equal(sor_ted[rank_s], un_sorted)
rank_s <- sor_t[sor_t]




#########
### Scripts for calling RcppArmadillo functions for manipulating vectors and matrices

# Compile Rcpp functions
Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/armadillo_functions.cpp")


## sum_na() sum_if() conditional sums Rcpp functions

# Create synthetic data
vec_tor <- 1:100
vec_tor[sample(1:100, 5)] <- NA

sum(is.na(vec_tor))
sum_na(vec_tor)
sum_na_stl(vec_tor)

# Benchmark Rcpp sum_na() function
library(microbenchmark)
summary(microbenchmark(
  sum_na=sum_na(vec_tor),
  sum_na_stl=sum_na_stl(vec_tor),
  sum_is_na=sum(is.na(vec_tor)),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# sum(is.na()) is 5 times faster than Rcpp
#         expr    mean median
# 1     sum_na 3778.50   3910
# 2 sum_na_stl 3592.92   3422
# 3  sum_is_na  728.87    490

sum_if(vec_tor, 5)
sum_if_cpp(vec_tor, 5)
sum_if_stl(vec_tor, 5)
sum(vec_tor < 5)

summary(microbenchmark(
  sum_if_cpp=sum_if_cpp(vec_tor, 5),
  sum_if=sum_if(vec_tor, 5),
  sum_if_stl=sum_if_stl(vec_tor, 5),
  r_code=sum(vec_tor < 5),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# sum(vec_tor < 5) is over 2 times faster than Rcpp
#         expr    mean median
# 1 sum_if_cpp 2424.76   2444
# 2     sum_if 2419.90   2444
# 3 sum_if_stl 2185.26   1956
# 4     r_code 1056.44    978


## which() Rcpp functions

# Create synthetic data
vec_tor <- round(runif(16), 2)
mat_rix <- matrix(round(runif(16), 2), nc=4)
bool_ean <- sample(c(TRUE, rep(FALSE, 9)), size=1e3, replace=TRUE)

# whi_ch3(bool_ean)
all.equal(whi_ch3(bool_ean), whi_ch4(bool_ean))

# Benchmark Rcpp which functions
library(microbenchmark)
summary(microbenchmark(
  whi_ch32=whi_ch32(bool_ean),
  whi_ch33=whi_ch33(bool_ean),
  whi_ch34=whi_ch34(bool_ean),
  whi_ch=whi_ch(bool_ean),
  whi_ch2=whi_ch2(bool_ean),
  whi_ch4=whi_ch4(bool_ean),
  whi_ch3=whi_ch3(bool_ean),
  which=which(bool_ean),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: which() is fastest followed by whi_ch3():
#      expr     mean median
# 1 whi_ch5 59.05335 57.181
# 2  whi_ch  7.16539  6.843
# 3 whi_ch2  8.05976  7.332
# 4 whi_ch4  4.17935  3.911
# 5 whi_ch3  3.32402  2.934
# 6   which  2.28303  2.444


## select elements and assign values

# sub-matrix Rcpp functions
sub_mat(mat_rix=mat_rix, row_num=c(1, 3), col_num=1:2)
sub_mat(mat_rix=mat_rix, row_num=1:2, col_num=1:2)
sub_mat_cast(mat_rix=mat_rix, row_num=1:2, col_num=1:2)

library(microbenchmark)
# microbenchmark shows: sub_mat() is slightly faster
summary(microbenchmark(
  sub_mat=sub_mat(mat_rix=mat_rix, row_num=1:2, col_num=1:2),
  sub_mat_cast=sub_mat_cast(mat_rix=mat_rix, row_num=1:2, col_num=1:2),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary


select_sub_mat(mat_rix=mat_rix, 0.4, 0)
find_sub_mat(mat_rix=mat_rix, 0.4, 0)

summary(microbenchmark(
  select_sub_mat=select_sub_mat(mat_rix=mat_rix, 0.4, 0),
  find_sub_mat=find_sub_mat(mat_rix=mat_rix, 0.4, 0),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: both about the same
#             expr    mean median
# 1 select_sub_mat 3.21168  2.933
# 2   find_sub_mat 2.79620  2.445


# function to assign values to selected vector elements
sub_assign(vec_tor=vec_tor, in_dex=c(2, 4, 6), da_ta=c(3, 5, 7))
# function to find selected vector elements and to assign values
find_assign_vec(vec_tor=vec_tor, fi_nd=0.5, da_ta=0.1)
# function to find selected vector elements and to assign values
find_assign_vec_point(vec_tor=vec_tor, fi_nd=0.5, da_ta=0.1)
# Rcpp function to assign values to selected matrix elements
find_assign_mat(mat_rix=mat_rix, fi_nd=0.5, da_ta=1)
# Rcpp function to assign values to selected matrix elements
find_extract_mat(mat_rix=mat_rix, fi_nd=0.8)

library(microbenchmark)
summary(microbenchmark(
  in_place=find_assign_vec_point(vec_tor=vec_tor, fi_nd=0.5, da_ta=1),
  find_assign=find_assign_vec(vec_tor=vec_tor, fi_nd=0.5, da_ta=1),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# find_assign_vec_point() is slightly faster than find_assign_vec()
#          expr     mean  median
# 1    in_place 540.4299 543.465
# 2 find_assign 940.3258 815.197


# column compare Rcpp functions
compare_col(mat_rix=mat_rix, 0.5, 1)
compare_col_arma(mat_rix=mat_rix, 0.5)
compare_col_arma(mat_rix=mat_rix, 0.5, 1)
compare_col_armaa(mat_rix=mat_rix, 0.5, 1)

library(microbenchmark)
summary(microbenchmark(
  compare_col=compare_col(mat_rix=mat_rix, 0.5, 1),
  compare_col_arma=compare_col_arma(mat_rix=mat_rix, 0.5, 1),
  compare_col_armaa=compare_col_armaa(mat_rix=mat_rix, 0.5, 1),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: compare_col() is fastest
#                expr    mean median
# 1       compare_col 1525.52   1466
# 2  compare_col_arma 2366.22   2444
# 3 compare_col_armaa 2024.09   1956

# which column Rcpp function
which_col(mat_rix=mat_rix, 0.5, 2)
which_col(mat_rix=mat_rix, 0.5)



## Calculate the rolling sum over a vector

Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/roll_sum.cpp")

vec_tor <- rnorm(1e6)
all.equal(HighFreq::roll_sum(vec_tor, look_back=11)[-(1:10)],
          RcppRoll::roll_sum(vec_tor, n=11))
all.equal(HighFreq::roll_sum(vec_tor, look_back=11),
          drop(roll_sum_arma(vec_tor, look_back=11)))
all.equal(rutils::roll_sum(matrix(vec_tor, nc=1), look_back=11),
          roll_sum_arma(vec_tor, look_back=11))

library(microbenchmark)
summary(microbenchmark(
  rcpp=HighFreq::roll_sum(vec_tor, look_back=11),
  RcppRoll=RcppRoll::roll_sum(vec_tor, n=11),
  arma=roll_sum_arma(vec_tor, look_back=11),
  r_code=rutils::roll_sum(matrix(vec_tor, nc=1), look_back=11),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# HighFreq::roll_sum() in pure Rcpp is faster than RcppRoll::roll_sum()
# which is faster than roll_sum_arma().
# HighFreq::roll_sum() in pure Rcpp is 
# over 6 times faster than rutils::roll_sum() in vectorized R
#       expr      mean    median
# 1     rcpp  9.165521  6.592201
# 2 RcppRoll  9.919361  9.940151
# 3     arma 18.360971 16.809851
# 4   r_code 66.861341 54.078750


## Calculate the rolling weighted sum over a vector

Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/roll_wsum.cpp")


vec_tor <- as.numeric(rutils::env_etf$VTI[, 6])
wei_ghts <- exp(-0.2*1:11)
wei_ghts <- wei_ghts/sum(wei_ghts)
# compare R with Rcpp
weight_ed <- HighFreq::roll_wsum(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
filter_ed <- filter(x=vec_tor, filter=wei_ghts, method="convolution", sides=1)
all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))

# compare Rcpp with RcppArmadillo: agrees for filter wei_ghts <- c(1, rep(1e-5, 10))
#  but different for exponentially decaying weights
filter_ed <- roll_wsum_armaa(vec_tor, rev(wei_ghts))
all.equal(as.numeric(filter_ed), as.numeric(weight_ed))
round(as.numeric(tail(weight_ed, 22)), 2)
round(as.numeric(tail(filter_ed, 22)), 2)

library(microbenchmark)
summary(microbenchmark(
  rcpp=HighFreq::roll_wsum(vec_tor=vec_tor, wei_ghts=rev(wei_ghts)),
  r_code=filter(x=vec_tor, filter=wei_ghts, method="convolution", sides=1, circular=TRUE),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary


# microbenchmark shows: 
# roll_wsum() in pure Rcpp is a little slower than roll_wsum_arma()
# and about 5 times faster than filter() in vectorized R
#     expr     mean  median
# 1   rcpp  91.4411  76.975
# 2 r_code 513.0144 459.646



## Calculate the lag of a vector

Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/roll_var.cpp")

vec_tor <- rnorm(1e6)

all.equal(drop(lagvec_arma(vec_tor, lagg=11)),
          lagvec_rcpp(vec_tor, lagg=11))
all.equal(drop(lagvec_arma(vec_tor, lagg=11)),
          rutils::lag_it(vec_tor, lagg=11))

library(microbenchmark)
summary(microbenchmark(
  rcpp=lagvec_rcpp(vec_tor, lagg=11),
  arma=lagvec_arma(vec_tor, lagg=11),
  lag_it=rutils::lag_it(vec_tor, lagg=11),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# lagvec_rcpp() in pure Rcpp is about twice as fast as lagvec_arma(), 
# which is slightly faster than the vectorized function lag_it() in native R
#     expr     mean   median
# 1   rcpp  5.22094  2.98555
# 2   arma  9.62885  7.05545
# 3 lag_it 12.56096 10.60035


## Calculate the variance of a vector


all.equal(variance_rcpp(vec_tor),
          variance_arma(vec_tor))

library(microbenchmark)
summary(microbenchmark(
  arma=variance_arma(vec_tor),
  rcpp=variance_rcpp(vec_tor),
  r_code=var(vec_tor),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# variance_arma() is slightly faster than variance_rcpp() 
# in pure Rcpp, and is more than twice as fast as the 
# vectorized function var() in native R
#     expr     mean   median
# 1   arma 1.767471 1.747152
# 2   rcpp 2.262741 2.273301
# 3 r_code 4.944931 4.780251


## Calculate the rolling variance over a vector

Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/roll_var.cpp")

re_turns <- rnorm(1e3)

all.equal(roll_var_rcpp(re_turns, look_back=21),
          drop(roll_var_arma(re_turns, look_back=21)))
all.equal(drop(roll_var_arma(re_turns, look_back=21))[-(1:20)],
          RcppRoll::roll_var(re_turns, n=21))

library(microbenchmark)
summary(microbenchmark(
  arma=roll_var_arma(re_turns, look_back=21),
  rcpp=roll_var_rcpp(re_turns, look_back=21),
  RcppRoll=RcppRoll::roll_var(re_turns, n=21),
  # r_code=var(re_turns),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# The RcppArmadillo function roll_var_arma() is twice as fast 
# as the pure Rcpp function roll_var_rcpp() which is faster 
# than RcppRoll::roll_sum().
# 
#       expr      mean    median
# 1     arma  70.57456  63.20685
# 2     rcpp 200.11045 182.29265
# 3 RcppRoll 273.72994 261.03385


## Calculate the rolling variance over a column matrix

re_turns <- na.omit(rutils::etf_env$re_turns$VTI)

all.equal(drop(HighFreq::roll_var(re_turns, ste_p=1, look_back=11))[-(1:10)],
          drop(RcppRoll::roll_var(re_turns, n=11)), 
          check.attributes=FALSE)
# end_p <- rutils::calc_endpoints(re_turns, inter_val=5)
end_p <- drop(HighFreq::calc_endpoints(NROW(re_turns), ste_p=5)) + 1
n_rows <- NROW(end_p)
# start_p <- c(0, 0, 0, end_p[1:(n_rows-3)])
start_p <- drop(HighFreq::calc_startpoints(end_p, look_back=3))
foo <- sapply(1:n_rows, function(ep) {
  var(re_turns[start_p[ep]:end_p[ep]])
})  # end sapply
bar <- HighFreq::roll_var(re_turns, ste_p=5, look_back=3)
all.equal(foo, drop(bar), check.attributes=FALSE)

## Calculate the variance of OHLC

oh_lc <- rutils::etf_env$VTI
all.equal(variance_ohlc(oh_lc),
          HighFreq::calc_variance_ohlc_r(oh_lc))

## Calculate the rolling variance over OHLC
oh_lc <- rutils::etf_env$VTI[1:31, ]
in_dex <- c(1, diff(xts::.index(oh_lc)))
var_rolling <- HighFreq::roll_var_ohlc(oh_lc,
                                       ste_p=1,
                                       look_back=2,
                                       calc_method="yang_zhang",
                                       scal_e=TRUE,
                                       in_dex=in_dex)

var_rollingn <- roll_var_ohlc(oh_lc, 
                             ste_p=1,
                             look_back=2,
                             calc_method="yang_zhang",
                             scal_e=TRUE,
                             in_dex=in_dex)

var_rollingp <- roll_var_ohlcp(oh_lc, 
                              look_back=2,
                              calc_method="yang_zhang",
                              scal_e=TRUE,
                              in_dex=in_dex)

n_rows <- NROW(oh_lc)
lag_close = HighFreq::lag_it(oh_lc[, 4])
# var_rollingr <- numeric(n_rows)
end_p <- drop(HighFreq::calc_endpoints(n_rows, 1)) + 1
start_p <- drop(HighFreq::calc_startpoints(end_p, 2))
n_pts <- NROW(end_p)
var_rollingr <- sapply(2:n_pts, function(it) {
  ran_ge <- start_p[it]:end_p[it]
  sub_ohlc = oh_lc[ran_ge, ]
  sub_close = lag_close[ran_ge]
  sub_index = in_dex[ran_ge]
  HighFreq::calc_var_ohlc(sub_ohlc, lag_close=sub_close, in_dex=sub_index, scal_e=TRUE)
})  # end sapply
var_rollingr <- c(0, var_rollingr)
foo <- cbind(var_rolling, var_rollingp, var_rollingn[, 1], var_rollingr)
tail(foo)


# Multistep
var_rolling <- HighFreq::roll_var_ohlc(oh_lc,
                                       ste_p=3,
                                       look_back=2,
                                       calc_method="yang_zhang",
                                       scal_e=TRUE,
                                       in_dex=in_dex)
end_p <- drop(HighFreq::calc_endpoints(n_rows, 3)) + 1
start_p <- drop(HighFreq::calc_startpoints(end_p, 2))
n_pts <- NROW(end_p)
var_rollingr <- sapply(2:n_pts, function(it) {
  ran_ge <- start_p[it]:end_p[it]
  sub_ohlc = oh_lc[ran_ge, ]
  sub_close = lag_close[ran_ge]
  sub_index = in_dex[ran_ge]
  HighFreq::calc_var_ohlc(sub_ohlc, lag_close=sub_close, in_dex=sub_index, scal_e=TRUE)
})  # end sapply
var_rollingr <- c(0, var_rollingr)
all.equal(drop(var_rolling), var_rollingr)
foo <- cbind(var_rolling, var_rollingr)
tail(foo)





## microbenchmark
library(microbenchmark)
summary(microbenchmark(
  arma=variance_arma(oh_lc),
  rcpp=variance_rcpp(oh_lc),
  r_code=var(oh_lc),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# variance_arma() is slightly faster than variance_rcpp() 
# in pure Rcpp, and is more than twice as fast as the 
# vectorized function var() in native R
#     expr     mean   median
# 1   arma 1.767471 1.747152
# 2   rcpp 2.262741 2.273301
# 3 r_code 4.944931 4.780251





#########
### Scripts for calling RcppArmadillo functions for matrix algebra

## de-mean the columns of a matrix

summary(microbenchmark(
  demean_mat=demean_mat(mat_rix),
  demean_arma=demean_arma(mat_rix),
  apply=(apply(mat_rix, 2, mean)),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# demean_mat() is over 5 times faster than demean_arma()
# and over 20 times faster than apply()
#           expr      mean    median
# 1   demean_mat  1.206325  1.188584
# 2  demean_arma  9.909479  5.964911
# 3        apply 44.555462  25.05482


## bind the columns of two matrices

mat_rix1 <- matrix(runif(1e6), nc=1e3)
mat_rix2 <- matrix(runif(1e6), nc=1e3)
# cbind(mat_rix1, mat_rix2)
all.equal(cbind_rcpp(mat_rix1, mat_rix2), cbind_arma(mat_rix1, mat_rix2))

summary(microbenchmark(
  cbind_arma=cbind_arma(mat_rix1, mat_rix2),
  cbind_rcpp=cbind_rcpp(mat_rix1, mat_rix2),
  cbind=cbind(mat_rix1, mat_rix2),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# cbind_rcpp() is as fast as cbind(), and more that 2 times faster than
# cbind_arma().
#         expr      mean    median
# 1 cbind_arma 12.893332 12.414150
# 2 cbind_rcpp  5.943275  4.813715
# 3      cbind  5.829133  4.906573


## Calculate the inner (dot) product of two vectors.

vec1 <- runif(1e5)
vec2 <- runif(1e5)

vec_prod(vec1, vec2)
vec1 %*% vec2

summary(microbenchmark(
  vec_prod=vec_prod(vec1, vec2),
  r_code=(vec1 %*% vec2),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# vec_prod() is several times faster than %*%, especially for longer vectors.
#       expr     mean   median
# 1   vec_prod 110.7067 110.4530
# 2   r_code 585.5127 591.3575


## Calculate the product of a matrix times a vector.

mat_rix <- matrix(runif(1e7), nc=1e5)
vec1 <- runif(1e5)
all.equal(mat_rix %*% vec1, mat_vec_prod(vec_tor=vec1, mat_rix=mat_rix))

summary(microbenchmark(
  mat_vec_prod=mat_vec_prod(vec_tor=vec1, mat_rix=mat_rix),
  r_code=(mat_rix %*% vec1),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# mat_vec_prod() is 3 times faster than %*%, for matrix with 100,000 columns.
#        expr      mean    median
# 1 mat_vec_prod  7.299448  7.180375
# 2    r_code 21.133891 21.048730


vec2 <- runif(1e2)
all.equal(drop(vec2 %*% (mat_rix %*% vec1)), mat_2vec_prod(vec2, mat_rix, vec1))

summary(microbenchmark(
  mat_2vec_prod=mat_2vec_prod(vec2, mat_rix, vec1),
  r_code=(vec2 %*% (mat_rix %*% vec1)),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# mat_2vec_prod() is 3 times faster than %*%, for matrix with 100,000 columns.
#            expr      mean    median
# 1   mat_2vec_prod  7.138696  7.071877
# 2        r_code 20.826379 20.678520


## Calculate product of matrix and vectors
# multiply the matrix elements *by* the vector elements

# Multiply matrix columns
mat_rix <- matrix(round(runif(25e4), 2), nc=5e2)
vec_tor <- round(runif(5e2), 2)
prod_uct <- vec_tor*mat_rix
mult_vec_mat(vec_tor, mat_rix)
all.equal(prod_uct, mat_rix)
summary(microbenchmark(
    r_cpp=mult_vec_mat(vec_tor, mat_rix),
    r_code=vec_tor*mat_rix,
    times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# Multiply matrix rows
mat_rix <- matrix(round(runif(25e4), 2), nc=5e2)
vec_tor <- round(runif(5e2), 2)
prod_uct <- t(vec_tor*t(mat_rix))
mult_vec_mat(vec_tor, mat_rix, by_col=FALSE)
all.equal(prod_uct, mat_rix)
library(microbenchmark)
summary(microbenchmark(
    r_cpp=mult_vec_mat(vec_tor, mat_rix, by_col=FALSE),
    r_code=t(vec_tor*t(mat_rix)),
    times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# multiply the matrix elements *by* the elements of two vectors

mat_rix <- matrix(runif(1e7), nc=1e5)
vec1 <- runif(NCOL(mat_rix))
vec2 <- runif(NROW(mat_rix))
prod_uct <- t(t(vec2*mat_rix)*vec1)
mult_vec2_mat_copy(vec2, mat_rix, vec1)
all.equal(mat_rix, prod_uct)

summary(microbenchmark(
  mult_vec2_mat_copy=mult_vec2_mat_copy(vec2, mat_rix, vec1),
  mult_vec2_mat_rcpp=mult_vec2_mat_rcpp(vec2, mat_rix, vec1),
  mult_vec2_mat_rcpp2=mult_vec2_mat_rcpp2(vec2, mat_rix, vec1),
  r_code=(t(t(vec2*mat_rix)*vec1)),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# mult_vec2_mat_copy() is over 2 times faster than %*% and t(), for matrix
# with 100,000 columns.
#                expr      mean    median
# 1 mult_vec2_mat_copy  73.65367  73.50842
# 2  mult_vec2_mat_rcpp 101.39165 100.44875
# 3 mult_vec2_mat_rcpp2 612.48159 612.98899
# 4            r_code 182.74140 174.80584



## matrix inversion correlation PCA

# Create random matrix
mat_rix <- matrix(rnorm(500), nc=5)
# Calculate correlation matrix
matrix_cor <- cor(mat_rix)
matrix_cor <- get_cor(mat_rix)

library(microbenchmark)
summary(microbenchmark(
  cor_mat=cor(mat_rix),
  cor_arma=get_cor(mat_rix),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# get_cor() is over 4 times faster than cor()
#       expr     mean  median
# 1  cor_mat 36.06386 32.7460
# 2 cor_arma 21.12839  7.8205


# Calculate eigen matrix
eigen_r <- eigen(cor(mat_rix))
ei_gen <- calc_eigen(mat_rix)
all.equal(eigen_r$values, sort(drop(ei_gen$values), decreasing=TRUE))
all.equal(abs(eigen_r$vectors), 
          abs(ei_gen$vectors)[, match(round(abs(eigen_r$vectors[1, ]), 3), round(abs(ei_gen$vectors[1, ]), 3))])

summary(microbenchmark(
  eigen_r=eigen(cor(mat_rix)),
  eigen_arma=calc_eigen(mat_rix),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# calc_eigen() is about 18 times faster than eigen() plus cor()
#         expr      mean   median
# 1    eigen_r 435.81817 429.3475
# 2 eigen_arma  39.67553  23.9480


# de-mean (center) and scale the mat_rix columns
mat_rix <- t(t(mat_rix) - colMeans(mat_rix))
mat_rix <- t(t(mat_rix) / sqrt(colSums(mat_rix^2)/(NROW(mat_rix)-1)))
# Calculate PCA
pc_a <- get_pca(mat_rix)

summary(microbenchmark(
  eigen_r=eigen(cor(mat_rix)),
  eigen_arma=calc_eigen(mat_rix),
  pc_a=get_pca(mat_rix),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# calc_eigen() is about 6 times faster than get_pca()
#         expr      mean   median
# 1    eigen_r 436.10160 429.5920
# 2 eigen_arma  24.15352  23.4595
# 3       pc_a 148.45158 145.6415



# Create random positive semi-definite matrix
mat_rix <- matrix(runif(25), nc=5)
mat_rix <- t(mat_rix) %*% mat_rix
# perform matrix inversion
matrix_inv <- solve(mat_rix)
matrix_inv <- invspd_arma(mat_rix)
matrix_inv <- invspd_rcpp(mat_rix)

library(microbenchmark)
summary(microbenchmark(
  inv_mat=inv_mat(mat_rix),
  invspd_arma=invspd_arma(mat_rix),
  invspd_rcpp=invspd_rcpp(mat_rix),
  solve=solve(mat_rix),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# inv_mat() is over 10 times faster than solve()
# invspd_arma() is over 7 times faster than solve()
#                expr     mean median
# 1           inv_mat  3.42669  2.933
# 2       invspd_arma  4.68759  3.911
# 3       invspd_rcpp  4.74625  3.911
# 4             solve 32.00254 31.280


## Multivariate linear regression

# define design matrix with explanatory variables
len_gth <- 100
n_var <- 5
de_sign <- matrix(rnorm(n_var*len_gth), nc=n_var)
# response equals linear form plus error terms
noise <- rnorm(len_gth, sd=0.5)
weight_s <- rnorm(n_var)
res_ponse <- -3 + de_sign %*% weight_s + noise
# perform multivariate regression using lm()
reg_model <- lm(res_ponse ~ de_sign)
coef(reg_model)
sum_mary <- summary(reg_model)
# perform multivariate regression using HighFreq::calc_lm()
reg_model_arma <- HighFreq::calc_lm(res_ponse=res_ponse, de_sign=de_sign)
reg_model_arma$coefficients
all.equal(reg_model_arma$coefficients[, "coeff"], unname(coef(reg_model)))
all.equal(unname(reg_model_arma$coefficients), unname(sum_mary$coefficients))
all.equal(drop(reg_model_arma$residuals), unname(reg_model$residuals))
all.equal(unname(reg_model_arma$stats), c(sum_mary$r.squared, unname(sum_mary$fstatistic[1])))

# library(MASS)
# multivariate regression using MASS::ginv() in lm_r()
lm_r <- function(res_ponse, de_sign) {
  # solve for regression betas
  de_sign <- cbind(rep(1, NROW(de_sign)), de_sign)
  beta_s <- MASS::ginv(de_sign) %*% res_ponse
  fit_ted <- drop(de_sign %*% beta_s)
  # Calculate residuals
  resid_uals <- drop(res_ponse - fit_ted)
  # variance of residuals
  deg_free <- len_gth-NCOL(de_sign)
  resid_var <- sum(resid_uals^2)/deg_free
  # explanatory matrix squared
  explain_squared <- crossprod(de_sign)
  # Calculate covariance matrix of betas
  beta_covar <- resid_var*MASS::ginv(explain_squared)
  beta_sd <- sqrt(diag(beta_covar))
  # Calculate t-values of betas
  beta_tvals <- drop(beta_s)/beta_sd
  # Calculate two-sided p-values of betas
  beta_pvals <- 2*pt(-abs(beta_tvals), df=deg_free)
  cbind(beta_s, beta_sd, beta_tvals, beta_pvals)
}  # end lm_r
lm_r(res_ponse, de_sign)

library(microbenchmark)
summary(microbenchmark(
  lm_arma=HighFreq::calc_lm(res_ponse, de_sign),
  lm_r=lm_r(res_ponse, de_sign),
  lm=lm(res_ponse ~ de_sign),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# HighFreq::calc_lm() is over 10 times faster than lm() and over 3 times faster than
# lm_r().
#      expr       mean    median
# 1 lm_arma   99.31485   98.4795
# 2    lm_r  328.17102  324.5155
# 3      lm 1070.44432 1036.8345


## Calculate Z-scores from rolling multivariate regression using RcppArmadillo

# Calculate Z-scores from rolling time series regression using RcppArmadillo
look_back <- 11
clo_se <- Cl(rutils::env_etf$VTI)
date_s <- xts::.index(clo_se)

z_scores <- HighFreq::roll_zscores(res_ponse=as.numeric(clo_se), 
                                   de_sign=matrix(as.numeric(date_s), nc=1), 
                                   look_back=look_back)
# Plot dygraph
hist(z_scores)
dygraphs::dygraph(xts(z_scores, index(clo_se)), 
                  main="Z-scores of VTI time series regressions")

# Calculate Z-scores from rolling multivariate regression using RcppArmadillo

z_scores <- HighFreq::roll_zscores(res_ponse=res_ponse, de_sign=de_sign, look_back=look_back)

# Calculate z-scores in R from rolling multivariate regression using lm()
z_scores_r <- sapply(1:NROW(de_sign), function(ro_w) {
  if (ro_w==1) return(0)
  start_point <- max(1, ro_w-look_back+1)
  sub_response <- res_ponse[start_point:ro_w]
  sub_design <- de_sign[start_point:ro_w, ]
  reg_model <- lm(sub_response ~ sub_design)
  resid_uals <- reg_model$residuals
  resid_uals[NROW(resid_uals)]/sd(resid_uals)
})  # end sapply
all.equal(unname(z_scores[-(1:look_back)]), unname(z_scores_r[-(1:look_back)]))


library(microbenchmark)
summary(microbenchmark(
  lm_arma=HighFreq::roll_zscores(res_ponse=res_ponse, de_sign=de_sign, look_back=look_back),
  lm_r=sapply(1:NROW(de_sign), function(ro_w) {
    if (ro_w==1) return(0)
    start_point <- max(1, ro_w-look_back+1)
    sub_response <- res_ponse[start_point:ro_w]
    sub_design <- de_sign[start_point:ro_w, ]
    reg_model <- lm(sub_response ~ sub_design)
    resid_uals <- reg_model$residuals
    resid_uals[NROW(resid_uals)]/sd(resid_uals)
  }),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# roll_zscores() is about 18 times faster than sapply().
#      expr       mean    median
# 1 lm_arma   5.668367  5.612995
# 2    lm_r 101.479011 98.942891



## split-apply-combine procedure

# Create synthetic data
vec_tor <- sample(1:5, 1e3, replace=TRUE)
fac_tor <- sample(1:5, 1e3, replace=TRUE)
mat_rix <- matrix(runif(2e3), nc=2)
mat_rix <- cbind(vec_tor, mat_rix)

# The function tapply_arma() performs aggregations over a vector using a factor.
# It produces the same result as the R code: 
#   tapply(X=vec_tor, INDEX=fac_tor, FUN=NROW)

tapply_arma(vec_tor, fac_tor)
tapply(X=vec_tor, INDEX=fac_tor, FUN=NROW)
all.equal(drop(tapply_arma(vec_tor, fac_tor)), as.numeric(tapply(X=vec_tor, INDEX=fac_tor, FUN=NROW)))

summary(microbenchmark(
  tapply_arma=tapply_arma(vec_tor, fac_tor),
  tapply=tapply(X=vec_tor, INDEX=fac_tor, FUN=NROW),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# tapply_arma() is almost 4 times faster than tapply(), for vector with 1,000
# elements, but it loses its advantage for longer vectors.
#          expr      mean   median
# 1 tapply_arma  45.19329  44.4750
# 2      tapply 178.35619 173.0095



# The function apply_agg() performs aggregations over a matrix using its
# first column as a factor.
# It produces the same result as the R code: 
#     sapply(X=unique(mat_rix[, 1]), FUN=function(mat_rix[, -1]))

tapply(X=mat_rix[, 2], INDEX=mat_rix[, 1], FUN=mean)

all.equal(sort(apply_agg(mat_rix)), 
          sort(sapply(X=unique(mat_rix[, 1]), FUN=function(x) {
            foo <- mat_rix[which(mat_rix[, 1] == x), -1,  drop=FALSE]
            sum(apply(foo, 1, prod))
            # sum(foo)
          })))

summary(microbenchmark(
  apply_agg=apply_agg(mat_rix),
  sapply=sapply(X=unique(mat_rix[, 1]), FUN=function(x) {
    foo <- mat_rix[which(mat_rix[, 1] == x), -1,  drop=FALSE]
    sum(apply(foo, 1, prod))
    # sum(foo)
  }),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# apply_agg() is over 40 times faster than sapply(), for matrix with 1,000 
# rows.
#        expr       mean   median
# 1 apply_agg   66.95125   54.494
# 2    sapply 2433.73195 2351.748

