////////////////////////////
// Functions to test C++ syntax with Armadillo
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_arma.cpp")

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

// using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;


////////////////////////////////////////////////////////////
//' The function order_index() calculates the order index of a vector
//' using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::uvec order_index(const arma::vec& vectorv) {
  return arma::sort_index(vectorv);
}  // end order_index

// [[Rcpp::export]]
arma::uvec calc_endpoints(arma::uword length,  // The length of the sequence
                          arma::uword step = 1,  // The number of periods between neighboring end points
                          arma::uword stub = 0,  // The first non-zero end point
                          bool stubs = true) {  // Include stub intervals?
  
  // Number of initial end points
  arma::uword numpts = length / step + 3;
  // Define the end points
  arma::uvec endd;
  endd.zeros(numpts);
  // Define the last end point
  arma::uword lastp = length - 1;
  
  // Calculate the initial end points - including extra end points at the end
  if (stub == 0) {
    for (arma::uword it = 0; it < numpts; ++it) {
      endd(it) = it*step;
    }  // end for
  } else if ((stub > 0) & (stubs)) {
    for (arma::uword it = 1; it < numpts; ++it) {
      endd(it) = stub + (it-1)*step;
    }  // end for
  } else {
    for (arma::uword it = 0; it < numpts; ++it) {
      endd(it) = stub + it*step;
    }  // end for
  }  // end if
  // std::cout << "endd = " << arma::conv_to<rowvec>::from(endd) << std::endl;
  
  // arma::uvec endd = arma::regspace<uvec>(stub, step, lastp + step);
  // Find the index of the largest element of endd which is less than lastp
  arma::uword endpp = 0;
  for (arma::uword it = 0; endd(it) < lastp; ++it) {
    endpp++;
  }  // end for
  // std::cout << "endpp = " << endpp << std::endl;
  
  // Trim the initial end points at the end - remove extra end points at the end
  // Subset endd to include the smallest element of endd which is equal to or greater than lastp
  endd = endd.subvec(0, endpp);
  
  // Set the stub intervals at the end
  if (stubs) {
    // Include stub intervals
    // Set the last end point to lastp - last element of endd
    endd[endpp] = lastp;
  } else {
    // Do not include the last end point - no stub interval at the end
    // Exclude the last element greater than lastp
    if (endd[endpp] > lastp) {
      endd = endd.subvec(0, endpp-1);
    }  // end if
  }  // end if
  
  return endd;
  
  // Old code below
  // if ((stub == 0) & (remainp == 0)) {
  //   std::cout << "(stub == 0) & (remainp == 0)" << std::endl;
  //   // No stub interval
  //   endd = arma::regspace<uvec>(0, step, length);
  //   endd.back() = lastp;
  //   // endd.transform([](arma::uword val) {return (val - 1);});
  //   // endd.front() = 0;
  // } else if ((stub == 0) & (remainp > 0)) {
  //   std::cout << "(stub == 0) & (remainp > 0)" << std::endl;
  //   // Add stub interval at end
  //   endd = arma::regspace<uvec>(0, step, length + step);
  //   // endd.transform([](arma::uword val) {return (val - 1);});
  //   // endd.front() = 0;
  //   endd.back() = lastp;
  // } else if ((stub > 0) & (remainp == 0)) {
  //   std::cout << "(stub > 0) & (remainp == 0)" << std::endl;
  //   // Add initial stub interval equal to stub
  //   arma::uvec endm = arma::regspace<uvec>(stub, step, length);
  //   endm.back() = lastp;
  //   endd.zeros(numpts+1);
  //   endd.subvec(1, numpts) = endm;
  // } else if ((stub > 0) & (remainp > 0)) {
  //   std::cout << "(stub > 0) & (remainp > 0)" << std::endl;
  //   // Add initial stub interval equal to stub
  //   arma::uvec endm = arma::regspace<uvec>(stub, step, lastp + step);
  //   endm.back() = lastp;
  //   endd.zeros(numpts+1);
  //   endd.subvec(1, numpts) = endm;
  // }  // end if
  
}  // end calc_endpoints




//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& matrixv, const arma::uword& max_eigen) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  
  arma::eig_sym(eigen_val, eigen_vec, cov(matrixv));
  eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
  eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
  // arma::mat eigen_valmat = diagmat(eigen_val);
  
  return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  
}  // end calc_inv
// [[Rcpp::export]]
arma::mat diffit(const arma::mat& tseries, arma::sword lagg = 1, bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_rows-1);
  arma::uword ncols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    // Matrix difference without padding
    arma::mat diffmat = (tseries.rows(lagg, nrows) - tseries.rows(0, nrows - lagg));
    if (pad_zeros) {
      // Pad diffmat with zeros at the front
      return arma::join_cols(arma::zeros<mat>(lagg, ncols), diffmat);
    } else {
      // Don't pad the output
      return diffmat;
    }  // end if pad_zeros
  } else {
    // Negative lag
    // Matrix difference without padding
    arma::mat diffmat = (tseries.rows(0, nrows + lagg) - tseries.rows(-lagg, nrows));
    if (pad_zeros) {
      // Pad diffmat with zeros at the back
      return arma::join_cols(diffmat, arma::zeros<mat>(-lagg, ncols));
    } else {
      // Don't pad the output
      return diffmat;
    }  // end if pad_zeros
  }  // end if lagg
  
}  // end diffit



// [[Rcpp::export]]
arma::mat calc_var_ag(const arma::mat& pricev, 
                      arma::uword step = 1) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return arma::var(diffit(pricev, 1, false));
  else {
    // Allocate aggregations, end points, and variance.
    arma::uword nrows = pricev.n_rows;
    arma::mat aggs;
    arma::uvec endd;
    // The number of rows is equal to step so that loop works for stub=0
    arma::mat vars(step, pricev.n_cols);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub < step; stub++) {
      endd = calc_endpoints(nrows, step, stub, false);
      // endd = arma::regspace<uvec>(stub, step, nrows + step);
      // endd = endd.elem(find(endd < nrows));
      aggs = pricev.rows(endd);
      vars.row(stub) = arma::var(diffit(aggs, 1, false));
    }  // end for
    return mean(vars);
  }  // end if
  
}  // end calc_var_ag


// [[Rcpp::export]]
arma::mat calc_hurst(const arma::mat& tseries, 
                     const arma::vec& aggv) {
  
  // If only single agg value then calculate the Hurst exponent from a single data point
  if (aggv.n_rows == 1) {
    return 0.5*arma::log(calc_var_ag(tseries, aggv(0))/calc_var_ag(tseries, 1))/log(aggv(0));
  }  // end if
  
  // Allocate the objects
  arma::uword nrows = aggv.n_rows;
  // cout << "nrows: " << nrows << endl;
  // cout << "aggv: " << aggv << endl;
  arma::mat volv(nrows, tseries.n_cols, fill::zeros);
  
  // Calculate the log volatility at the agg points
  for (arma::uword it = 0; it < nrows; it++) {
    volv.row(it) = 0.5*arma::log(calc_var_ag(tseries, aggv(it)));
  }  // end for
  
  // Calculate the log of the agg points
  arma::mat agglog = arma::log(aggv);
  // cout << "agglog: " << agglog << endl;
  
  // Calculate the Hurst exponent from the regression slopes
  arma::mat varagg = arma::var(agglog);
  // cout << "varagg: " << varagg << endl;
  return (arma::cov(volv, agglog))/varagg(0);
  
}  // end calc_hurst



//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(const arma::vec& respv,  // Response vector
                   const arma::mat& predm) { // Predictor matrix
  
  // Add column for intercept to the predictor matrix
  arma::uword nrows = predm.n_rows;
  // arma::mat predm = arma::join_rows(ones(nrows), predm);
  arma::uword ncols = predm.n_cols;
  arma::uword degf = (nrows - ncols);
  
  // Calculate alpha and beta coefficients for the model response ~ predictor
  arma::colvec coeff = arma::solve(predm, respv);
  // Calculate residuals
  arma::colvec residuals = respv - predm*coeff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (nrows-1)*arma::var(respv);
  double res_sumsq = arma::dot(residuals, residuals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate R-squared and F-statistic
  double rsquared = exp_sumsq/tot_sumsq;
  double fstat = (exp_sumsq*degf)/(res_sumsq*(ncols-1));
  // arma::rowvec stats=join_horiz(rsquared, fstat);
  Rcpp::NumericVector stats(2);
  stats(0) = rsquared;
  stats(1) = fstat;
  stats.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // Calculate standard errors of betac coefficients
  arma::colvec stderrv = arma::sqrt(res_sumsq/degf*arma::diagvec(arma::pinv(arma::trans(predm)*predm)));
  // Calculate t-values and p-values of betac coefficients
  arma::colvec tvals = coeff/stderrv;
  arma::colvec pvals = 2*Rcpp::pt(-Rcpp::abs(Rcpp::wrap(tvals)), degf);
  Rcpp::NumericMatrix coeffmat = Rcpp::wrap(arma::join_rows(arma::join_rows(arma::join_rows(coeff, stderrv), tvals), pvals));
  Rcpp::colnames(coeffmat) = Rcpp::CharacterVector::create("coeff", "stderr", "tvals", "pvals");
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = coeffmat,
                            // Named("residuals") = residuals,
                            Rcpp::Named("zscore") = residuals(nrows-1)/arma::stddev(residuals),
                            Rcpp::Named("stats") = stats);
  
}  // end calc_lm



//' @export
// [[Rcpp::export]]
double calc_alpha(arma::mat& returns, 
                  const arma::mat& indeks, 
                  const std::string& typev = "jensen",
                  const double& betav = 1.0) {
  // Initialize
  double alpha = 0;

  // Calculate alpha depending on typev
  if (typev == "jensen") {
    // Mean returns by columns
    Rcpp::NumericMatrix coeff = calc_lm(returns, indeks)["coefficients"];
    return coeff(0, 0);
  } else if (typev == "wilcoxon") {
    returns = (returns - indeks);
    // arma::uword nrows = (returns.n_rows);
    arma::uvec ranks = (arma::sort_index(arma::sort_index(returns)) + 1);
    return  dot(sign(returns), ranks);
  } else if (typev == "kruskal_wallis") {
    arma::uword nrows = (returns.n_rows);
    arma::mat combined = join_cols(returns, indeks);
    arma::uvec ranks = (arma::sort_index(arma::sort_index(combined)) + 1);
    // Apply regularized inverse to unit vector
    // weights = calc_inv(returns, max_eigen)*arma::ones(returns.n_cols);
    return  (0.0 + sum(ranks.subvec(0, nrows-1)) - sum(ranks.subvec(nrows, 2*nrows-1)));
  } else if (typev == "rank") {
    // Mean returns by columns
    // arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    // arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    // sd_cols.replace(0, 1);
    // meancols = meancols/sd_cols;
    // alpha equal to ranks of Sharpe
    // weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    // weights = (weights - arma::mean(weights));
    return alpha;
  } else if (typev == "rankrob") {
    // Mean returns by columns
    // arma::vec meancols = arma::trans(arma::median(returns, 0));
    // meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, max_eigen);
    // weights = calc_inv(returns, max_eigen)*meancols;
    // weights = calc_inv(returns, max_eigen)*meancols;
    // // Standard deviation by columns
    // arma::vec sd_cols = meancols;
    // for (arma::uword it=0; it < returns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((returns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // meancols = meancols/sd_cols;
    // alpha equal to ranks of Sharpe
    // weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    // weights = conv_to< vec >::from(arma::sort_index(meancols));
    // probv;
    // weights = (weights - arma::mean(weights));
    return alpha;
  } else {
    cout << "Warning: Incorrect typev argument: " << typev << endl;
    return alpha;
  }  // end if
  
  return alpha;
  
}  // end calc_alpha


// Calculate the top and bottom quantiles of the columns of a matrix
//' @export
// [[Rcpp::export]]
arma::vec calc_quant(const arma::mat& returns, 
                     const double& probv = 0.1) {
  
  arma::vec probs = {probv, 1-probv};
  // return arma::quantile(returns, probs, 0);
  return conv_to< colvec >::from(arma::sum(arma::quantile(returns, probs, 0), 0));
  
}  // end calc_quant



// Calculate the top and bottom quantiles of a vector,
// and return a vector of zeros, ones, and minus ones, 
// with the top quantile elements equal to 1, and the 
// bottom equal to -1.
//' @export
// [[Rcpp::export]]
arma::vec calc_top_bottom(const arma::vec& returns, 
                          const double& probv = 0.1) {
  
  // arma::uword ndataem = returns.n_elem;
  arma::vec weights = zeros(returns.n_elem);
  arma::vec probs = {probv, 1-probv};
  arma::vec quantiles = arma::quantile(returns, probs);
  
  // Bottom quantile
  arma::uvec indeks = find(returns <= quantiles(0));
  weights(indeks).fill(-1);
  
  // Top quantile
  // arma::rowvec quantilev = arma::quantile(returns, probv);
  indeks = find(returns >= quantiles(1));
  // weights.elem(indeks).ones();
  weights(indeks).fill(1);
  
  return weights;
  
}  // end calc_top_bottom



// Calculate the top and bottom quantiles of the columns of a matrix.
// Add the quantiles for each column, and rank the vector of sums.
// Return a vector of zeros, ones, and minus ones, with the top 
// ranks equal to 1, and the bottom equal to -1.
//' @export
// [[Rcpp::export]]
arma::vec calc_top_bottom_columns(const arma::vec& returns, 
                      const double& probv = 0.1) {
  
  arma::uword ncols = returns.n_cols;
  arma::uword threshold = round(probv*ncols);
  arma::vec probs = {probv, 1-probv};
  arma::vec weights = zeros(ncols);

  // return arma::quantile(returns, probs, 0);
  arma::rowvec sum_quant = arma::sum(arma::quantile(returns, probs, 0), 0);
  arma::uvec ranks = (arma::sort_index(arma::sort_index(sum_quant)));
  // arma::uvec indeks = find((ranks <= threshold) || (ranks >= (ncols - threshold)));
  arma::uvec indeks = find(ranks >= (ncols - threshold));
  // weights.elem(indeks).ones();
  weights(indeks).fill(1);
  indeks = find(ranks <= threshold);
  weights(indeks).fill(-1);
  
  // indeks = find((ranks > threshold) || (ranks < (ncols - threshold)));
  // weights.elem(indeks) = 0;
  // weights.head(threshold) = 1;
  // weights.tail(threshold) = 1;
  // return conv_to< colvec >::from(weights);
  return weights;
  
}  // end calc_top_bottom_columns



//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& returns, 
                       const std::string& typev = "quantilev",
                       int max_eigen = 1,
                       const double& probv = 0.1,
                       const double& alpha = 0.0,
                       const bool scalit = true) {
  // Initialize
  arma::vec weights(returns.n_cols);
  if (max_eigen == 1)  max_eigen = returns.n_cols;
  
  // Calculate weights depending on typev
  if (typev == "max_sharpe") {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, max_eigen);
    weights = calc_inv(returns, max_eigen)*meancols;
  } else if (typev == "quantilev") {
    // Sum of quantiles by columns
    arma::vec probs = {probv, 1-probv};
    weights = conv_to< vec >::from(arma::sum(arma::quantile(returns, probs, 0), 0));
    // Weights equal to ranks
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(weights)));
    weights = (weights - arma::mean(weights));
  } else if (typev == "min_var") {
    // Apply regularized inverse to unit vector
    weights = calc_inv(returns, max_eigen)*arma::ones(returns.n_cols);
  } else if (typev == "min_varpca") {
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, cov(returns));
    weights = eigen_vec.col(0);
  } else if (typev == "rank") {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
  } else if (typev == "rankrob") {
    // Median returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Apply regularized inverse
    // arma::mat inverse = calc_inv(returns, max_eigen);
    // weights = calc_inv(returns, max_eigen)*meancols;
    // weights = calc_inv(returns, max_eigen)*meancols;
    // // Standard deviation by columns
    // arma::vec sd_cols = meancols;
    // for (arma::uword it=0; it < returns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((returns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to< vec >::from(arma::sort_index(arma::sort_index(meancols)));
    // probv;
    weights = (weights - arma::mean(weights));
  } else {
    cout << "Warning: Incorrect typev argument: " << typev << endl;
    return arma::ones(returns.n_cols);
  }  // end if
  
  if (scalit == TRUE) {
    // Returns of equally weighted portfolio
    // arma::vec meanrows = arma::mean(returns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = returns*weights;
    // Scale weights to equally weighted portfolio and return them
    // Return weights/sqrt(sum(square(weights)));
    // Return weights/sum(weights);
    return weights*arma::stddev(arma::mean(returns, 1))/arma::stddev(returns*weights);
  }  // end if
  
  return weights;
}  // end calc_weights


