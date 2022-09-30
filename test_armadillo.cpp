////////////////////////////
// Functions to test C++ syntax with Armadillo
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/test_armadillo.cpp")

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Use STL
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


////////////////////////////////////////////////
// Tests misc

// Example of transform
//' @export
// [[Rcpp::export]]
arma::mat floor_it(arma::mat& data, double minval) {
  
  arma::mat mydata = data;
  mydata.transform([&minval](double x) {return max(x, minval);});
  return mydata;
  
}  // end floor_it


////////////////////////////////////////////////
// Included to facilitate Tests - remove later - don't migrate


arma::uvec calc_endpoints(arma::uword length,  // The length of the sequence
                          arma::uword step = 1,  // The number of periods between neighboring end points
                          arma::uword stub = 0,  // The first non-zero end point
                          bool stubs = true) {  // Include stub intervals?
  
  // Number of initial end points
  arma::uword numpts = length / step + 3;
  // Define the end points
  arma::uvec endp;
  endp.zeros(numpts);
  // Define the last end point
  int lastp = length - 1;
  
  // Calculate the initial end points - including extra end points at the end
  if (stub == 0) {
    for (arma::uword it = 0; it < numpts; ++it) {
      endp[it] = it*step;
    }  // end for
  } else if ((stub > 0) & (stubs)) {
    for (arma::uword it = 1; it < numpts; ++it) {
      endp[it] = stub + (it-1)*step;
    }  // end for
  } else {
    for (arma::uword it = 0; it < numpts; ++it) {
      endp[it] = stub + it*step;
    }  // end for
  }  // end if
  // std::cout << "endp = " << arma::conv_to<rowvec>::from(endp) << std::endl;
  
  // arma::uvec endp = arma::regspace<uvec>(stub, step, lastp + step);
  // Find the index of the largest element of endp which is less than lastp
  arma::uword endpp = 0;
  for (arma::uword it = 0; endp[it] < lastp; ++it) {
    endpp++;
  }  // end for
  // std::cout << "endpp = " << endpp << std::endl;
  
  // Trim the initial end points at the end - remove extra end points at the end
  // Subset endp to include the smallest element of endp which is equal to or greater than lastp
  endp = endp.subvec(0, endpp);
  
  // Set the stub intervals at the end
  if (stubs) {
    // Include stub intervals
    // Set the last end point to lastp - last element of endp
    endp[endpp] = lastp;
  } else {
    // Do not include the last end point - no stub interval at the end
    // Exclude the last element greater than lastp
    if (endp[endpp] > lastp) {
      endp = endp.subvec(0, endpp-1);
    }  // end if
  }  // end if
  
  return endp;
  
}  // end calc_endpoints


arma::uvec calc_startpoints(arma::uvec endp, arma::uword look_back) {
  
  arma::uword numpts = endp.n_elem;
  arma::uvec startp = arma::join_cols(arma::zeros<uvec>(look_back), 
                                      endp.subvec(0, numpts - look_back - 1) + 1);
  
  return startp;
  
}  // end calc_startpoints



enum methodenum {moment, least_squares, quantile, nonparametric, regular, sharpem, 
                 maxsharpe, maxsharpemed, minvarlin, minvarquad, kellym, robustm, 
                 sumsq, sumone, voltarget, voleqw};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
methodenum calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return methodenum::moment;
  else if (method == "least_squares" || method == "l")
    return methodenum::least_squares;
  else if (method == "quantile" || method == "q")
    return methodenum::quantile;
  else if (method == "nonparametric" || method == "n")
    return methodenum::nonparametric;
  else if (method == "regular")
    return methodenum::regular;
  else if (method == "sharpem")
    return methodenum::sharpem;
  else if (method == "maxsharpe")
    return methodenum::maxsharpe;
  else if (method == "maxsharpemed")
    return methodenum::maxsharpemed;
  else if (method == "minvarlin")
    return methodenum::minvarlin;
  else if (method == "minvarquad")
    return methodenum::minvarquad;
  else if (method == "kellym")
    return methodenum::kellym;
  else if (method == "robustm")
    return methodenum::robustm;
  else if (method == "sumsq")
    return methodenum::sumsq;
  else if (method == "voltarget")
    return methodenum::voltarget;
  else 
    return methodenum::moment;
}  // end calc_method


// [[Rcpp::export]]
arma::mat calc_mean(const arma::mat& tseries,
                    std::string method = "moment", 
                    double confl = 0.75) {
  
  // Apply different calculation methods for location
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::mean(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(0) + quantiles.row(1));
  }  // end quantile
  case methodenum::nonparametric: {  // nonparametric
    return arma::median(tseries);
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_mean



arma::mat calc_var(const arma::mat& tseries,
                   std::string method = "moment", 
                   double confl = 0.75) {
  
  double ncols = tseries.n_cols;
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(ncols);
  }  // end if
  
  // Apply different calculation methods for dispersion
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::var(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(1) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // MAD
    arma::mat medians = arma::median(tseries);
    arma::mat mads(1, ncols);
    // Loop over columns of tseries
    for (arma::uword it = 0; it < ncols; it++) {
      mads.col(it) = arma::median(arma::abs(tseries.col(it) - arma::as_scalar(medians.col(it))));
    }  // end for
    // tseries.each_row() -= arma::median(tseries, 0);
    // return 1.4826*arma::median(arma::abs(tseries), 0);
    return 1.4826*mads;
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(ncols);
  }  // end default
  }  // end switch
  
}  // end calc_var


arma::mat calc_skew(const arma::mat& tseries,
                    std::string method = "moment", 
                    double confl = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Apply different calculation methods for skew
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    double nrows = tseries.n_rows;
    double ncols = tseries.n_cols;
    arma::mat means = arma::mean(tseries);
    arma::mat vars = arma::var(tseries);
    arma::mat skewness(1, ncols);
    // De-mean the columns of tseries
    // tseries.each_row() -= means;
    // return arma::sum(arma::pow(tseries, 3))/arma::pow(vars, 1.5)/nrows;
    for (arma::uword it = 0; it < ncols; it++) {
      skewness.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(means.col(it)), 3))/arma::pow(vars.col(it), 1.5)/nrows;
    }  // end for
    return skewness;
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, 0.5, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(2) + quantiles.row(0) - 2*quantiles.row(1))/(quantiles.row(2) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // nonparametric
    return (arma::mean(tseries) - arma::median(tseries))/arma::stddev(tseries);
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_skew


arma::mat calc_kurtosis(const arma::mat& tseries,
                        std::string method = "moment", 
                        double confl = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Apply different calculation methods for kurtosis
  switch(calc_method(method)) {
  case methodenum::moment: {  // Fourth moment
    double nrows = tseries.n_rows;
    double ncols = tseries.n_cols;
    arma::mat means = arma::mean(tseries);
    arma::mat vars = arma::var(tseries);
    arma::mat kurtosis(1, ncols);
    // Don't de-mean the columns of tseries because that requires copying the matrix of data, so it's time-consuming
    // Loop over columns of tseries
    for (arma::uword it = 0; it < ncols; it++) {
      kurtosis.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(means.col(it)), 4))/arma::pow(vars.col(it), 2)/nrows;
    }  // end for
    // tseries.each_row() -= means;
    // return arma::sum(arma::pow(tseries, 4))/arma::pow(vars, 2)/nrows;
    return kurtosis;
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, 0.25, 0.75, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(3) - quantiles.row(0))/(quantiles.row(2) - quantiles.row(1));
  }  // end quantile
  case methodenum::nonparametric: {  // nonparametric
    return (arma::mean(tseries) - arma::median(tseries))/arma::stddev(tseries);
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_kurtosis


arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword dimax = 0) {
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD
  arma::svd(svdu, svdval, svdv, tseries);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > eigen_thresh*arma::sum(svdval));
  
  if (dimax == 0) {
    // Set dimax
    dimax = svdnum - 1;
  } else {
    // Adjust dimax
    dimax = std::min(dimax - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, dimax);
  svdu = svdu.cols(0, dimax);
  svdv = svdv.cols(0, dimax);
  
  // Calculate the regularized inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_inv


// [[Rcpp::export]]
void parselist(Rcpp::List controlv) {
  // Add intercept column to the predictor matrix?
  bool intercept = as<int>(controlv["intercept"]);
  // Type of the regression model
  std::string method = as<std::string>(controlv["method"]);
  // Threshold level for discarding small singular values
  double eigen_thresh = as<double>(controlv["eigen_thresh"]);
  // Dimension reduction
  arma::uword dimax = as<int>(controlv["dimax"]);
  // Confidence level for calculating the quantiles of returns
  double confl = as<double>(controlv["confl"]);
  // Return shrinkage intensity
  double alpha = as<double>(controlv["alpha"]);
  cout << "intercept = " << intercept << endl;
  cout << "method = " << method << endl;
  cout << "eigen_thresh = " << eigen_thresh << endl;
  cout << "dimax = " << dimax << endl;
  cout << "confl = " << confl << endl;
  cout << "alpha = " << alpha << endl;
}  // end parselist


////////////////////////////////////////////////
// Test versions to be migrated to package HighFreq

// [[Rcpp::export]]
Rcpp::List calc_control(bool intercept = true,  // Add intercept column to the predictor matrix?
                        std::string method = "least_squares",  // Type of the regression model
                        double eigen_thresh = 1e-5, // Threshold level for discarding small singular values
                        arma::uword dimax = 0, // Dimension reduction
                        double confl = 0.1, // Confidence level for calculating the quantiles of returns
                        double alpha = 0.0) {  // Return shrinkage intensity
  
  Rcpp::List controlv = Rcpp::List::create(Named("intercept") = intercept,
                                           Named("method") = method,
                                           Named("eigen_thresh") = eigen_thresh,
                                           Named("dimax") = dimax,
                                           Named("confl") = confl,
                                           Named("alpha") = alpha);
  
  return controlv;
  
}  // end calc_control


// [[Rcpp::export]]
arma::mat calc_reg(const arma::mat& response, 
                   const arma::mat& predictor,
                   Rcpp::List controlv) {

  // Unpack control list
  // Type of the regression model
  std::string method = as<std::string>(controlv["method"]);
  // Add intercept column to the predictor matrix?
  bool intercept = as<int>(controlv["intercept"]);
  // Threshold level for discarding small singular values
  double eigen_thresh = as<double>(controlv["eigen_thresh"]);
  // Dimension reduction
  arma::uword dimax = as<int>(controlv["dimax"]);
  // Confidence level for calculating the quantiles of returns
  double confl = as<double>(controlv["confl"]);
  // Return shrinkage intensity
  double alpha = as<double>(controlv["alpha"]);
  
  // Add column for intercept to predictor matrix
  arma::uword nrows = predictor.n_rows;
  arma::mat predictori = predictor;
  if (intercept)
    predictori = join_rows(ones(nrows), predictor);
  
  arma::uword ncols = predictori.n_cols;
  arma::uword deg_free = (nrows - ncols);
  arma::vec coeff;
  arma::vec tvals;
  arma::mat reg_data = arma::zeros<mat>(2*ncols+1, 1);
  
  // Apply different calculation methods for weights
  switch(calc_method(method)) {
  case methodenum::least_squares: {
    // Calculate regression coefficients for the model response ~ predictor
    coeff = arma::solve(predictori, response);
    break;
  }  // end least_squares
  case methodenum::regular: {
    // Calculate shrinkage regression coefficients
    coeff = calc_inv(predictori, eigen_thresh, dimax)*response;
    break;
  }  // end regular
  case methodenum::quantile: {
    // Not implemented yet
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return reg_data;
  }  // end default
  }  // end switch
  
  // Calculate residuals
  arma::mat residuals = response - predictori*coeff;
  
  // Calculate TSS, RSS, and ESS
  // double tot_sumsq = (nrows-1)*arma::var(response);
  double res_sumsq = arma::dot(residuals, residuals);
  // double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate standard errors of the beta coefficients
  arma::mat stderrv = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(predictori)*predictori)));
  // Calculate t-values of the beta coefficients
  tvals = coeff/stderrv;
  
  // Calculate z-score
  mat zscore = residuals(nrows-1, 0)/arma::stddev(residuals);
  
  // Combine regression data
  reg_data.rows(0, ncols-1) = coeff;
  reg_data.rows(ncols, 2*ncols-1) = tvals;
  reg_data.row(2*ncols) = zscore;
  
  return reg_data.t();
  
}  // end calc_reg


// [[Rcpp::export]]
arma::mat roll_reg(const arma::mat& response, 
                   const arma::mat& predictor,
                   Rcpp::List controlv, 
                   arma::uvec startp = 0, 
                   arma::uvec endp = 0, 
                   arma::uword step = 1, 
                   arma::uword look_back = 1, 
                   arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword nrows = predictor.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate regression matrix
  arma::mat responsi;
  arma::mat predicti;
  arma::uword numpts = endpts.n_elem;
  arma::uword ncols = predictor.n_cols;
  // Add intercept column to the predictor matrix?
  bool intercept = as<int>(controlv["intercept"]);
  if (intercept) ncols += 1;
  arma::mat reg_stats(numpts, (2*ncols + 1), fill::zeros);
  
  // Perform loop over the endpts
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate regression coefficients
    if (endpts(ep) > startpts(ep)) {
      responsi = response.rows(startpts(ep), endpts(ep));
      predicti = predictor.rows(startpts(ep), endpts(ep));
      reg_stats.row(ep) = calc_reg(responsi, predicti, controlv);
    }  // end if
  }  // end for
  
  // Warmup period
  // reg_stats.rows(0, ncols+1) = zeros(ncols+2, (ncols + 1));
  // for (arma::uword it = (ncols+2); it < look_back; it++) {
  //   responsi = response.rows(0, it);
  //   predicti = predictor.rows(0, it);
  //   reg_data = calc_reg(responsi, predicti);
  //   reg_stats.row(it) = arma::conv_to<rowvec>::from(reg_data);
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = look_back; it < nrows; it++) {
  //   responsi = response.rows(it-look_back+1, it);
  //   predicti = predictor.rows(it-look_back+1, it);
  //   reg_data = calc_reg(responsi, predicti, method, eigen_thresh, dimax, confl, alpha);
  //   reg_stats.row(it) = arma::conv_to<rowvec>::from(reg_data);
  // }  // end for
  
  return reg_stats;
  
}  // end roll_reg






////////////////////////////////////////////////
// Old stuff - can be deleted later

// Calculate the maximum of two vectors
//' @export
// [[Rcpp::export]]
arma::colvec max2(arma::colvec tseries, arma::colvec tseries2) {
  
  arma::uword nrows = tseries.n_rows;
  arma::colvec maxs = tseries;
  
  // Perform loop over rows
  for (arma::uword it = 0; it < nrows; it++) {
    maxs.row(it) = arma::max(tseries.row(it), tseries2.row(it));
  }  // end for
  
  return maxs;
  
}  // end max2

