////////////////////////////
// Utility functions using the Armadillo library
////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// Use STL
// using namespace std;
// #include <functional>
// #include <string>
// #include <map>
// #include <iostream>
// [[Rcpp::plugins(cpp11)]]


arma::vec lag_vec(const arma::vec& tseries, arma::sword lagg = 1, bool pad_zeros = true);

arma::mat lagit(const arma::mat& tseries, arma::sword lagg = 1, bool pad_zeros = true);

arma::mat diffit(const arma::mat& tseries, arma::sword lagg = 1, bool pad_zeros = true);

arma::vec diff_vec(const arma::vec& tseries, arma::uword lagg = 1, bool pad_zeros = true);

arma::uvec calc_endpoints(arma::uword length, arma::uword step = 1, arma::uword stub = 0);

arma::uvec calc_startpoints(arma::uvec endp, arma::uword look_back);

arma::uvec calc_ranks(arma::vec tseries);

arma::uvec calc_ranks_stl(arma::vec tseries);

double sum_it(arma::vec tseries);

arma::uvec roll_count(const arma::uvec& tseries);

Rcpp::List encode_it(arma::vec tseries);

std::vector<double> decode_it(Rcpp::List encodel);

