////////////////////////////
// Function roll_skew() for calculating the rolling skewness estimator 
// (third moment) using RcppArmadillo.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="C:/Develop/R/Rcpp/roll_skew.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


////////////////////////////////////////////////////////////
// Switch statement in calc_skew() uses C++ enum type.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
enum skew_type {Pearson, Quantile, Nonparametric};
// Map string to C++ enum type for switch statement.
skew_type calc_skew_type(const std::string& typ_e) {
  if (typ_e == "Pearson" || typ_e == "pearson" || typ_e == "p") 
    return skew_type::Pearson;
  else if (typ_e == "Quantile" || typ_e == "quantile" || typ_e == "q")
    return skew_type::Quantile;
  else if (typ_e == "Nonparametric" || typ_e == "nonparametric" || typ_e == "n")
    return skew_type::Nonparametric;
  else 
    return skew_type::Pearson;
}  // end calc_skew_type



////////////////////////////////////////////////////////////
//' Calculate the skewness of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{t_series} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{typ_e} A \emph{string} specifying the type of skewness (see
//'   Details). (The default is the \code{typ_e = "pearson"}.)
//'
//' @param \code{al_pha} The confidence level for calculating the quantiles.
//'   (the default is \code{al_pha = 0.25}).
//'
//' @return A single-row matrix with the skewness of the columns of
//'   \code{t_series}.
//'
//' @details The function \code{calc_skew()} calculates the skewness of the
//'   columns of a \emph{time series} or a \emph{matrix} of data using
//'   \code{RcppArmadillo} \code{C++} code.
//'
//'   If \code{typ_e = "pearson"} (the default) then \code{calc_skew()}
//'   calculates the Pearson skewness using the third moment of the data.
//'
//'   If \code{typ_e = "quantile"} then it calculates the skewness using the
//'   differences between the quantiles of the data.
//'
//'   If \code{typ_e = "nonparametric"} then it calculates the skewness as the
//'   difference between the mean of the data minus its median, divided by the
//'   standard deviation.
//'   
//'   The code examples below compare the function \code{calc_skew()} with the
//'   skewness calculated using \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Calculate VTI returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[ ,"VTI", drop=FALSE])
//' # Calculate the Pearson skewness
//' HighFreq::calc_skew(re_turns)
//' # Compare HighFreq::calc_skew() with Pearson skewness
//' calc_skewr <- function(x) {
//'   x <- (x-mean(x)); nr <- NROW(x);
//'   nr*sum(x^3)/(var(x))^1.5/(nr-1)/(nr-2)
//' }  # end calc_skewr
//' all.equal(HighFreq::calc_skew(re_turns), 
//'   calc_skewr(re_turns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns),
//'   Rcode=calc_skewr(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile skewness
//' HighFreq::calc_skew(re_turns, typ_e = "quantile", al_pha = 0.1)
//' # Compare HighFreq::calc_skew() with quantile skewness
//' calc_skewq <- function(x) {
//'   	quantile_s <- quantile(x, c(0.25, 0.5, 0.75), type=5)
//'   	(quantile_s[3] + quantile_s[1] - 2*quantile_s[2])/(quantile_s[3] - quantile_s[1])
//' }  # end calc_skewq
//' all.equal(drop(HighFreq::calc_skew(re_turns, typ_e = "quantile")), 
//'   calc_skewq(re_turns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns, typ_e = "quantile"),
//'   Rcode=calc_skewq(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric skewness
//' HighFreq::calc_skew(re_turns, typ_e = "nonparametric")
//' # Compare HighFreq::calc_skew() with R nonparametric skewness
//' all.equal(drop(HighFreq::calc_skew(re_turns, typ_e = "nonparametric")), 
//'   (mean(re_turns)-median(re_turns))/sd(re_turns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns, typ_e = "nonparametric"),
//'   Rcode=(mean(re_turns)-median(re_turns))/sd(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_skew(arma::mat t_series,
                    const std::string& typ_e = "pearson", 
                    double al_pha = 0.25) {
  
  // switch statement for all the different types of skew
  switch(calc_skew_type(typ_e)) {
  case skew_type::Pearson: {  // Pearson
    double num_rows = t_series.n_rows;
    arma::mat mean_s = arma::mean(t_series);
    arma::mat var_s = arma::var(t_series);
    // De-mean the columns of t_series
    t_series.each_row() -= mean_s;
    return (num_rows/(num_rows-1)/(num_rows-2))*arma::sum(arma::pow(t_series, 3))/arma::pow(var_s, 1.5);
  }  // end pearson
  case skew_type::Quantile: {  // Quantile
    arma::vec prob_s = {al_pha, 0.5, 1.0 - al_pha};
    arma::mat quantile_s = quantile(t_series, prob_s);
    return (quantile_s.row(2) + quantile_s.row(0) - 2*quantile_s.row(1))/(quantile_s.row(2) - quantile_s.row(0));
  }  // end quantile
  case skew_type::Nonparametric: {  // Nonparametric
    return (arma::mean(t_series) - arma::median(t_series))/arma::stddev(t_series);
  }  // end nonparametric
  default : {
    cout << "Invalid typ_e" << endl;
    return 0;
  }  // end default
  }  // end switch
  
}  // end calc_skew



//' @export
// [[Rcpp::export]]
arma::mat roll_skew(const arma::mat& t_series,
                    const std::string& typ_e = "pearson", 
                    double al_pha = 0.25, 
                    const arma::uword& look_back = 11) {
  
  arma::uword num_rows = t_series.n_rows;
  arma::uword num_cols = t_series.n_cols;
  arma::mat rolling_skew(num_rows, num_cols);
  arma::mat sub_series;
  // Rcpp::List lm_list;
  
  // Warmup period
  // rolling_skew.rows(0, num_cols+1) = arma::zeros(num_cols+2, (num_cols + 1));
  
  for (arma::uword it = 1; it < look_back; it++) {
    sub_series = t_series.rows(0, it);
    arma::mat skew_ness = calc_skew(sub_series, typ_e, al_pha);
    rolling_skew.row(it) = conv_to< rowvec >::from(skew_ness);
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < num_rows; it++) {
    sub_series = t_series.rows(it-look_back+1, it);
    arma::mat skew_ness = calc_skew(sub_series, typ_e, al_pha);
    rolling_skew.row(it) = conv_to< rowvec >::from(skew_ness);
  }  // end for
  
  return rolling_skew;
}  // end roll_skew


