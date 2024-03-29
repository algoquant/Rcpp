########################
### RcppArmadillo driver scripts
########################


#########
### Scripts for ad-hoc tests

# Compile Rcpp functions
Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_arma.cpp")



#########
### Scripts for calling RcppArmadillo functions

# Calculate the order index, to resort back to unsorted order
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
Rcpp::sourceCpp(file="Users/jerzy/Develop/Rcpp/armadillo_functions.cpp")


## sum_na() sum_if() conditional sums Rcpp functions

# Create synthetic data
vectorv <- 1:100
vectorv[sample(1:100, 5)] <- NA

sum(is.na(vectorv))
sum_na(vectorv)
sum_na_stl(vectorv)

# Benchmark Rcpp sum_na() function
library(microbenchmark)
summary(microbenchmark(
  sum_na=sum_na(vectorv),
  sum_na_stl=sum_na_stl(vectorv),
  sum_is_na=sum(is.na(vectorv)),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# sum(is.na()) is 5 times faster than Rcpp
#         expr    mean median
# 1     sum_na 3778.50   3910
# 2 sum_na_stl 3592.92   3422
# 3  sum_is_na  728.87    490

sum_if(vectorv, 5)
sum_if_cpp(vectorv, 5)
sum_if_stl(vectorv, 5)
sum(vectorv < 5)

summary(microbenchmark(
  sum_if_cpp=sum_if_cpp(vectorv, 5),
  sum_if=sum_if(vectorv, 5),
  sum_if_stl=sum_if_stl(vectorv, 5),
  rcode=sum(vectorv < 5),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# sum(vectorv < 5) is over 2 times faster than Rcpp
#         expr    mean median
# 1 sum_if_cpp 2424.76   2444
# 2     sum_if 2419.90   2444
# 3 sum_if_stl 2185.26   1956
# 4     rcode 1056.44    978


## which() Rcpp functions

# Create synthetic data
vectorv <- round(runif(16), 2)
matrixv <- matrix(round(runif(16), 2), nc=4)
booleanv <- sample(c(TRUE, rep(FALSE, 9)), size=1e3, replace=TRUE)

# whi_ch3(booleanv)
all.equal(whi_ch3(booleanv), whi_ch4(booleanv))

# Benchmark Rcpp which functions
library(microbenchmark)
summary(microbenchmark(
  whi_ch32=whi_ch32(booleanv),
  whi_ch33=whi_ch33(booleanv),
  whi_ch34=whi_ch34(booleanv),
  whi_ch=whi_ch(booleanv),
  whi_ch2=whi_ch2(booleanv),
  whi_ch4=whi_ch4(booleanv),
  whi_ch3=whi_ch3(booleanv),
  which=which(booleanv),
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
sub_mat(matrixv=matrixv, row_num=c(1, 3), col_num=1:2)
sub_mat(matrixv=matrixv, row_num=1:2, col_num=1:2)
sub_mat_cast(matrixv=matrixv, row_num=1:2, col_num=1:2)

library(microbenchmark)
# microbenchmark shows: sub_mat() is slightly faster
summary(microbenchmark(
  sub_mat=sub_mat(matrixv=matrixv, row_num=1:2, col_num=1:2),
  sub_mat_cast=sub_mat_cast(matrixv=matrixv, row_num=1:2, col_num=1:2),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary


select_sub_mat(matrixv=matrixv, 0.4, 0)
find_sub_mat(matrixv=matrixv, 0.4, 0)

summary(microbenchmark(
  select_sub_mat=select_sub_mat(matrixv=matrixv, 0.4, 0),
  find_sub_mat=find_sub_mat(matrixv=matrixv, 0.4, 0),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: both about the same
#             expr    mean median
# 1 select_sub_mat 3.21168  2.933
# 2   find_sub_mat 2.79620  2.445


# function to assign values to selected vector elements
sub_assign(vectorv=vectorv, in_dex=c(2, 4, 6), da_ta=c(3, 5, 7))
# function to find selected vector elements and to assign values
find_assign_vec(vectorv=vectorv, fi_nd=0.5, da_ta=0.1)
# function to find selected vector elements and to assign values
find_assign_vec_point(vectorv=vectorv, fi_nd=0.5, da_ta=0.1)
# Rcpp function to assign values to selected matrix elements
find_assign_mat(matrixv=matrixv, fi_nd=0.5, da_ta=1)
# Rcpp function to assign values to selected matrix elements
find_extract_mat(matrixv=matrixv, fi_nd=0.8)

library(microbenchmark)
summary(microbenchmark(
  in_place=find_assign_vec_point(vectorv=vectorv, fi_nd=0.5, da_ta=1),
  find_assign=find_assign_vec(vectorv=vectorv, fi_nd=0.5, da_ta=1),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# find_assign_vec_point() is slightly faster than find_assign_vec()
#          expr     mean  median
# 1    in_place 540.4299 543.465
# 2 find_assign 940.3258 815.197


# column compare Rcpp functions
compare_col(matrixv=matrixv, 0.5, 1)
compare_col_arma(matrixv=matrixv, 0.5)
compare_col_arma(matrixv=matrixv, 0.5, 1)
compare_col_armaa(matrixv=matrixv, 0.5, 1)

library(microbenchmark)
summary(microbenchmark(
  compare_col=compare_col(matrixv=matrixv, 0.5, 1),
  compare_col_arma=compare_col_arma(matrixv=matrixv, 0.5, 1),
  compare_col_armaa=compare_col_armaa(matrixv=matrixv, 0.5, 1),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: compare_col() is fastest
#                expr    mean median
# 1       compare_col 1525.52   1466
# 2  compare_col_arma 2366.22   2444
# 3 compare_col_armaa 2024.09   1956

# which column Rcpp function
which_col(matrixv=matrixv, 0.5, 2)
which_col(matrixv=matrixv, 0.5)



## Calculate the rolling sum over a vector

Rcpp::sourceCpp(file="Users/jerzy/Develop/Rcpp/roll_sum.cpp")

vectorv <- rnorm(1e6)
all.equal(HighFreq::roll_sum(vectorv, look_back=11)[-(1:10)],
          RcppRoll::roll_sum(vectorv, n=11))
all.equal(HighFreq::roll_sum(vectorv, look_back=11),
          drop(roll_sum_arma(vectorv, look_back=11)))
all.equal(rutils::roll_sum(matrix(vectorv, nc=1), look_back=11),
          roll_sum_arma(vectorv, look_back=11))

library(microbenchmark)
summary(microbenchmark(
  rcpp=HighFreq::roll_sum(vectorv, look_back=11),
  RcppRoll=RcppRoll::roll_sum(vectorv, n=11),
  arma=roll_sum_arma(vectorv, look_back=11),
  rcode=rutils::roll_sum(matrix(vectorv, nc=1), look_back=11),
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
# 4   rcode 66.861341 54.078750


## Calculate the rolling weighted sum over a vector

Rcpp::sourceCpp(file="Users/jerzy/Develop/Rcpp/roll_wsum.cpp")


vectorv <- as.numeric(rutils::env_etf$VTI[, 6])
wei_ghts <- exp(-0.2*1:11)
wei_ghts <- wei_ghts/sum(wei_ghts)
# compare R with Rcpp
weight_ed <- HighFreq::roll_wsum(vectorv=vectorv, wei_ghts=rev(wei_ghts))
filter_ed <- filter(x=vectorv, filter=wei_ghts, method="convolution", sides=1)
all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))

# compare Rcpp with RcppArmadillo: agrees for filter wei_ghts <- c(1, rep(1e-5, 10))
#  but different for exponentially decaying weights
filter_ed <- roll_wsum_armaa(vectorv, rev(wei_ghts))
all.equal(as.numeric(filter_ed), as.numeric(weight_ed))
round(as.numeric(tail(weight_ed, 22)), 2)
round(as.numeric(tail(filter_ed, 22)), 2)

library(microbenchmark)
summary(microbenchmark(
  rcpp=HighFreq::roll_wsum(vectorv=vectorv, wei_ghts=rev(wei_ghts)),
  rcode=filter(x=vectorv, filter=wei_ghts, method="convolution", sides=1, circular=TRUE),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary


# microbenchmark shows: 
# roll_wsum() in pure Rcpp is a little slower than roll_wsum_arma()
# and about 5 times faster than filter() in vectorized R
#     expr     mean  median
# 1   rcpp  91.4411  76.975
# 2 rcode 513.0144 459.646



## Calculate the lag of a vector

Rcpp::sourceCpp(file="Users/jerzy/Develop/Rcpp/roll_var.cpp")

vectorv <- rnorm(1e6)

all.equal(drop(lagvec_arma(vectorv, lagg=11)),
          lagvec_rcpp(vectorv, lagg=11))
all.equal(drop(lagvec_arma(vectorv, lagg=11)),
          rutils::lag_it(vectorv, lagg=11))

library(microbenchmark)
summary(microbenchmark(
  rcpp=lagvec_rcpp(vectorv, lagg=11),
  arma=lagvec_arma(vectorv, lagg=11),
  lag_it=rutils::lag_it(vectorv, lagg=11),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# lagvec_rcpp() in pure Rcpp is about twice as fast as lagvec_arma(), 
# which is slightly faster than the vectorized function lag_it() in native R
#     expr     mean   median
# 1   rcpp  5.22094  2.98555
# 2   arma  9.62885  7.05545
# 3 lag_it 12.56096 10.60035


## Calculate the variance of a vector


all.equal(variance_rcpp(vectorv),
          variance_arma(vectorv))

library(microbenchmark)
summary(microbenchmark(
  arma=variance_arma(vectorv),
  rcpp=variance_rcpp(vectorv),
  rcode=var(vectorv),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# variance_arma() is slightly faster than variance_rcpp() 
# in pure Rcpp, and is more than twice as fast as the 
# vectorized function var() in native R
#     expr     mean   median
# 1   arma 1.767471 1.747152
# 2   rcpp 2.262741 2.273301
# 3 rcode 4.944931 4.780251


## Calculate the rolling variance over a vector

Rcpp::sourceCpp(file="Users/jerzy/Develop/Rcpp/roll_var.cpp")

vectorv <- rnorm(1e6)
all.equal(roll_var_rcpp(vectorv, look_back=21),
          drop(roll_var_arma(vectorv, look_back=21)))
all.equal(drop(roll_var_arma(vectorv, look_back=21))[-(1:20)],
          RcppRoll::roll_var(vectorv, n=21))

library(microbenchmark)
summary(microbenchmark(
  arma=roll_var_arma(vectorv, look_back=21),
  rcpp=roll_var_rcpp(vectorv, look_back=21),
  RcppRoll=RcppRoll::roll_var(vectorv, n=21),
  # rcode=var(vectorv),
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



## Calculate the variance of OHLC

oh_lc <- rutils::etf_env$VTI
all.equal(variance_ohlc(oh_lc),
          HighFreq::calc_variance_ohlc_r(oh_lc))

library(microbenchmark)
summary(microbenchmark(
  arma=variance_arma(oh_lc),
  rcpp=variance_rcpp(oh_lc),
  rcode=var(oh_lc),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# variance_arma() is slightly faster than variance_rcpp() 
# in pure Rcpp, and is more than twice as fast as the 
# vectorized function var() in native R
#     expr     mean   median
# 1   arma 1.767471 1.747152
# 2   rcpp 2.262741 2.273301
# 3 rcode 4.944931 4.780251



#########
### Scripts for calling RcppArmadillo functions for matrix algebra

## de-mean the columns of a matrix

summary(microbenchmark(
  demean_mat=demean_mat(matrixv),
  demean_arma=demean_arma(matrixv),
  apply=(apply(matrixv, 2, mean)),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# demean_mat() is over 5 times faster than demean_arma()
# and over 20 times faster than apply()
#           expr      mean    median
# 1   demean_mat  1.206325  1.188584
# 2  demean_arma  9.909479  5.964911
# 3        apply 44.555462  25.05482


## bind the columns of two matrices

matrixv1 <- matrix(runif(1e6), nc=1e3)
matrixv2 <- matrix(runif(1e6), nc=1e3)
# cbind(matrixv1, matrixv2)
all.equal(cbind_rcpp(matrixv1, matrixv2), cbind_arma(matrixv1, matrixv2))

summary(microbenchmark(
  cbind_arma=cbind_arma(matrixv1, matrixv2),
  cbind_rcpp=cbind_rcpp(matrixv1, matrixv2),
  cbind=cbind(matrixv1, matrixv2),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# cbind_rcpp() is as fast as cbind(), and more that 2 times faster than
# cbind_arma().
#         expr      mean    median
# 1 cbind_arma 12.893332 12.414150
# 2 cbind_rcpp  5.943275  4.813715
# 3      cbind  5.829133  4.906573



## Calculate the outer product of two vectors.

summary(microbenchmark(
  Rcpp=outer_vec(bar, bar),
  Rcode=outer(bar, bar),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary


## Calculate the inner (dot) product of two vectors.

vec1 <- runif(1e5)
vec2 <- runif(1e5)

vec_prod(vec1, vec2)
vec1 %*% vec2

summary(microbenchmark(
  vec_prod=vec_prod(vec1, vec2),
  rcode=(vec1 %*% vec2),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# vec_prod() is several times faster than %*%, especially for longer vectors.
#       expr     mean   median
# 1   vec_prod 110.7067 110.4530
# 2   rcode 585.5127 591.3575


## Calculate the product of a matrix times a vector.

matrixv <- matrix(runif(1e7), nc=1e5)
vec1 <- runif(1e5)
all.equal(matrixv %*% vec1, mat_vec_prod(vectorv=vec1, matrixv=matrixv))

summary(microbenchmark(
  mat_vec_prod=mat_vec_prod(vectorv=vec1, matrixv=matrixv),
  rcode=(matrixv %*% vec1),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# mat_vec_prod() is 3 times faster than %*%, for matrix with 100,000 columns.
#        expr      mean    median
# 1 mat_vec_prod  7.299448  7.180375
# 2    rcode 21.133891 21.048730


vec2 <- runif(1e2)
all.equal(drop(vec2 %*% (matrixv %*% vec1)), mat_2vec_prod(vec2, matrixv, vec1))

summary(microbenchmark(
  mat_2vec_prod=mat_2vec_prod(vec2, matrixv, vec1),
  rcode=(vec2 %*% (matrixv %*% vec1)),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# mat_2vec_prod() is 3 times faster than %*%, for matrix with 100,000 columns.
#            expr      mean    median
# 1   mat_2vec_prod  7.138696  7.071877
# 2        rcode 20.826379 20.678520


## Calculate product of matrix and vectors
# multiply the matrix elements *by* the vector elements

# Multiply matrix columns
matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
vectorv <- round(runif(5e2), 2)
prod_uct <- vectorv*matrixv
mult_vec_mat(vectorv, matrixv)
all.equal(prod_uct, matrixv)
summary(microbenchmark(
    r_cpp=mult_vec_mat(vectorv, matrixv),
    rcode=vectorv*matrixv,
    times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# Multiply matrix rows
matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
vectorv <- round(runif(5e2), 2)
prod_uct <- t(vectorv*t(matrixv))
mult_vec_mat(vectorv, matrixv, by_col=FALSE)
all.equal(prod_uct, matrixv)
library(microbenchmark)
summary(microbenchmark(
    r_cpp=mult_vec_mat(vectorv, matrixv, by_col=FALSE),
    rcode=t(vectorv*t(matrixv)),
    times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# multiply the matrix elements *by* the elements of two vectors

matrixv <- matrix(runif(1e7), nc=1e5)
vec1 <- runif(NCOL(matrixv))
vec2 <- runif(NROW(matrixv))
prod_uct <- t(t(vec2*matrixv)*vec1)
mult_vec2_mat_copy(vec2, matrixv, vec1)
all.equal(matrixv, prod_uct)

summary(microbenchmark(
  mult_vec2_mat_copy=mult_vec2_mat_copy(vec2, matrixv, vec1),
  mult_vec2_mat_rcpp=mult_vec2_mat_rcpp(vec2, matrixv, vec1),
  mult_vec2_mat_rcpp2=mult_vec2_mat_rcpp2(vec2, matrixv, vec1),
  rcode=(t(t(vec2*matrixv)*vec1)),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# mult_vec2_mat_copy() is over 2 times faster than %*% and t(), for matrix
# with 100,000 columns.
#                expr      mean    median
# 1 mult_vec2_mat_copy  73.65367  73.50842
# 2  mult_vec2_mat_rcpp 101.39165 100.44875
# 3 mult_vec2_mat_rcpp2 612.48159 612.98899
# 4            rcode 182.74140 174.80584



## matrix inversion correlation PCA

# create random matrix
matrixv <- matrix(rnorm(500), nc=5)
# calculate correlation matrix
matrix_cor <- cor(matrixv)
matrix_cor <- get_cor(matrixv)

library(microbenchmark)
summary(microbenchmark(
  cor_mat=cor(matrixv),
  cor_arma=get_cor(matrixv),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# get_cor() is over 4 times faster than cor()
#       expr     mean  median
# 1  cor_mat 36.06386 32.7460
# 2 cor_arma 21.12839  7.8205


# calculate eigen matrix
eigen_r <- eigen(cor(matrixv))
ei_gen <- calc_eigen(matrixv)
all.equal(eigen_r$values, sort(drop(ei_gen$values), decreasing=TRUE))
all.equal(abs(eigen_r$vectors), 
          abs(ei_gen$vectors)[, match(round(abs(eigen_r$vectors[1, ]), 3), round(abs(ei_gen$vectors[1, ]), 3))])

summary(microbenchmark(
  eigen_r=eigen(cor(matrixv)),
  eigen_arma=calc_eigen(matrixv),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# calc_eigen() is about 18 times faster than eigen() plus cor()
#         expr      mean   median
# 1    eigen_r 435.81817 429.3475
# 2 eigen_arma  39.67553  23.9480


# de-mean (center) and scale the matrixv columns
matrixv <- t(t(matrixv) - colMeans(matrixv))
matrixv <- t(t(matrixv) / sqrt(colSums(matrixv^2)/(NROW(matrixv)-1)))
# calculate PCA
pc_a <- get_pca(matrixv)

summary(microbenchmark(
  eigen_r=eigen(cor(matrixv)),
  eigen_arma=calc_eigen(matrixv),
  pc_a=get_pca(matrixv),
  times=100))[, c(1, 4, 5)]  # end microbenchmark summary

# microbenchmark shows: 
# calc_eigen() is about 6 times faster than get_pca()
#         expr      mean   median
# 1    eigen_r 436.10160 429.5920
# 2 eigen_arma  24.15352  23.4595
# 3       pc_a 148.45158 145.6415



# create random positive semi-definite matrix
matrixv <- matrix(runif(25), nc=5)
matrixv <- t(matrixv) %*% matrixv
# perform matrix inversion
matrix_inv <- solve(matrixv)
matrix_inv <- invspd_arma(matrixv)
matrix_inv <- invspd_rcpp(matrixv)

library(microbenchmark)
summary(microbenchmark(
  inv_mat=inv_mat(matrixv),
  invspd_arma=invspd_arma(matrixv),
  invspd_rcpp=invspd_rcpp(matrixv),
  solve=solve(matrixv),
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
  # calculate residuals
  resid_uals <- drop(res_ponse - fit_ted)
  # variance of residuals
  deg_free <- len_gth-NCOL(de_sign)
  resid_var <- sum(resid_uals^2)/deg_free
  # explanatory matrix squared
  explain_squared <- crossprod(de_sign)
  # calculate covariance matrix of betas
  beta_covar <- resid_var*MASS::ginv(explain_squared)
  beta_sd <- sqrt(diag(beta_covar))
  # calculate t-values of betas
  beta_tvals <- drop(beta_s)/beta_sd
  # calculate two-sided p-values of betas
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

# calculate Z-scores from rolling time series regression using RcppArmadillo
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

# calculate Z-scores from rolling multivariate regression using RcppArmadillo

z_scores <- HighFreq::roll_zscores(res_ponse=res_ponse, de_sign=de_sign, look_back=look_back)

# calculate z-scores in R from rolling multivariate regression using lm()
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
vectorv <- sample(1:5, 1e3, replace=TRUE)
fac_tor <- sample(1:5, 1e3, replace=TRUE)
matrixv <- matrix(runif(2e3), nc=2)
matrixv <- cbind(vectorv, matrixv)

# The function tapply_arma() performs aggregations over a vector using a factor.
# It produces the same result as the R code: 
#   tapply(X=vectorv, INDEX=fac_tor, FUN=NROW)

tapply_arma(vectorv, fac_tor)
tapply(X=vectorv, INDEX=fac_tor, FUN=NROW)
all.equal(drop(tapply_arma(vectorv, fac_tor)), as.numeric(tapply(X=vectorv, INDEX=fac_tor, FUN=NROW)))

summary(microbenchmark(
  tapply_arma=tapply_arma(vectorv, fac_tor),
  tapply=tapply(X=vectorv, INDEX=fac_tor, FUN=NROW),
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
#     sapply(X=unique(matrixv[, 1]), FUN=function(matrixv[, -1]))

tapply(X=matrixv[, 2], INDEX=matrixv[, 1], FUN=mean)

all.equal(sort(apply_agg(matrixv)), 
          sort(sapply(X=unique(matrixv[, 1]), FUN=function(x) {
            foo <- matrixv[which(matrixv[, 1] == x), -1,  drop=FALSE]
            sum(apply(foo, 1, prod))
            # sum(foo)
          })))

summary(microbenchmark(
  apply_agg=apply_agg(matrixv),
  sapply=sapply(X=unique(matrixv[, 1]), FUN=function(x) {
    foo <- matrixv[which(matrixv[, 1] == x), -1,  drop=FALSE]
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

