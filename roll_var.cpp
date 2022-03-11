// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


//' @export
// [[Rcpp::export]]
arma::mat lagit(arma::mat& matrixv, int lagg=1) {
  int nrows = (matrixv.n_rows-1);
  
  if (lagg > 0)
    return arma::join_cols(arma::repelem(matrixv.row(0), lagg, 1), 
                           matrixv.rows(0, nrows-lagg));
  else
    return arma::join_cols(matrixv.rows(-lagg, nrows), 
                           arma::repelem(matrixv.row(nrows), -lagg, 1));
  
}  // end lagit



// The function variance_rcpp() calculates the variance of a vector using Rcpp.
//' @export
// [[Rcpp::export]]
double variance_rcpp(NumericVector& vectorv) {
  return sum(pow(vectorv - sum(vectorv)/vectorv.size(), 2))/(vectorv.size()-1);
}  // end variance_rcpp


// The function variance_arma() calculates the variance of a vector using 
// RcppArmadillo.
//' @export
// [[Rcpp::export]]
double variance_arma(arma::vec& vectorv) {
  return arma::var(vectorv);
}  // end variance_arma


//' @export
// [[Rcpp::export]]
double vec_prod(const arma::colvec& vec1, const arma::colvec& vec2) {
  return arma::dot(vec1, vec2);
}  // end vec_prod


//' @export
// [[Rcpp::export]]
double calc_mat(arma::mat ohlc) {
  return arma::sum(arma::var(ohlc));
}  // end calc_mat


//' @export
// [[Rcpp::export]]
double calc_matp(arma::mat& ohlc, 
                 arma::colvec lag_close=0, 
                 std::string calc_method="yang_zhang", 
                 arma::colvec indeks=0, 
                 bool scalit=true) {
  return arma::sum(arma::var(ohlc));
}  // end calc_matp


//' @export
// [[Rcpp::export]]
double variance_ohlc(arma::mat& ohlc) {
  int nrows = ohlc.n_rows;
  double coeff = 0.34/(1.34 + (nrows+1)/(nrows-1));
  arma::colvec clo_se = ohlc.col(3);
  arma::colvec open_close = arma::join_cols(clo_se.subvec(0, 0), clo_se.subvec(0, clo_se.n_elem-2));
  open_close = ohlc.col(0) - open_close;
  
  return arma::var(open_close) + coeff*arma::var(ohlc.col(0) - clo_se) +
    (1 - coeff)*(arma::dot((ohlc.col(1) - clo_se), (ohlc.col(1) - ohlc.col(0))) +
    arma::dot((ohlc.col(2) - clo_se), (ohlc.col(2) - ohlc.col(0))))/nrows;
  
}  // end variance_ohlc


// The below function with an explicit loop instead of a lag is actually 
// a bit slower than the above function with a lag
// double variance_ohlc_plus(arma::mat& ohlc) {
//   int nrows = ohlc.n_rows;
//   double coeff = 0.34/(1.34 + (nrows+1)/(nrows-1));
//   arma::colvec openp = ohlc.col(0);
//   arma::colvec clo_se = ohlc.col(3);
//   arma::colvec open_close(nrows);
//   // open_close = ohlc.col(0) - open_close;
// 
//   // So instead copy using an explicit loop:
//   for (int it = 1; it < nrows; it++) {
//     open_close(it) = (openp(it) - clo_se(it-1));
//   }  // end for
//   open_close(0) = (openp(0) - clo_se(0));
//   
//   return arma::var(open_close) + coeff*arma::var(ohlc.col(0) - clo_se) +
//     (1 - coeff)*(arma::dot((ohlc.col(1) - clo_se), (ohlc.col(1) - ohlc.col(0))) +
//     arma::dot((ohlc.col(2) - clo_se), (ohlc.col(2) - ohlc.col(0))))/nrows;
//   
// }  // end variance_ohlc_plus



//' @export
// [[Rcpp::export]]
double calc_var_ohlc(arma::mat& ohlc, 
                     std::string calc_method="yang_zhang", 
                     arma::colvec lag_close=0, 
                     arma::colvec indeks=0, 
                     bool scalit=true) {
  
  int nrows = ohlc.n_rows;
  double coeff = 0.34/(1.34 + (nrows+1)/(nrows-1));

  if (!scalit || (indeks.n_rows == 1)) {
    indeks = arma::ones(nrows);
    // cout << "ohlc.n_rows = " << nrows << endl;
    // cout << "indeks.n_rows = " << indeks.n_rows << endl;
  }  // end if
  
  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::colvec clo_se = ohlc.col(3);
  arma::colvec open_close(clo_se.n_rows);
  if (lag_close.n_rows == 1) {
    open_close = arma::join_cols(clo_se.subvec(0, 0), clo_se.subvec(0, clo_se.n_elem-2));
    open_close = (ohlc.col(0) - open_close)/indeks;
  } else {
    open_close = (ohlc.col(0) - lag_close)/indeks;
  }  // end if
  arma::colvec close_open = (clo_se - ohlc.col(0))/indeks;
  arma::colvec close_high = (clo_se - ohlc.col(1))/indeks;
  arma::colvec close_low = (clo_se - ohlc.col(2))/indeks;
  arma::colvec high_low = (ohlc.col(1) - ohlc.col(2))/indeks;
  arma::colvec high_open = (ohlc.col(1) - ohlc.col(0))/indeks;
  arma::colvec low_open = (ohlc.col(2) - ohlc.col(0))/indeks;
  
  cout << "ohlc.n_rows = " << nrows << endl;
  // cout << "ohlc.n_cols = " << ohlc.n_cols << endl;
  cout << "open_close.n_rows = " << open_close.n_rows << endl;
  // cout << "open_close.n_cols = " << open_close.n_cols << endl;
  // cout << "indeks.n_rows = " << indeks.n_rows << endl;
  // cout << "indeks.n_cols = " << indeks.n_cols << endl;
  // cout << "ohlc last row = " << ohlc.row(nrows-1) << endl;
  // cout << "sum clo_se = " << arma::var(arma::diff(clo_se)) << endl;

  if (calc_method == "close") {
    // cout << "Calc method is Close" << endl;
    return arma::var(arma::diff(clo_se));
  } else if (calc_method == "rogers_satchell") {
    // cout << "Calc method is Rogers-Satchell" << endl;
    return -(arma::dot(close_high, high_open) +
             arma::dot(close_low, low_open))/nrows;
  } else if (calc_method == "garman_klass") {
    // cout << "Calc method is Garman-Klass" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/nrows;
  } else if (calc_method == "garman_klass_yz") {
    // cout << "Calc method is Garman-Klass-YZ" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/nrows + 
            arma::var(open_close);
  } else if (calc_method == "yang_zhang") {
    // cout << "Calc method is Yang-Zhang" << endl;
    return arma::var(open_close) + coeff*arma::var(close_open) +
      (coeff-1)*(arma::dot(close_high, high_open) + 
      arma::dot(close_low, low_open))/nrows;
  } else {
    cout << "Wrong calc method!" << endl;
    return 1;
  }  // end if
  
  // cout << "Calc method is " << calc_method << endl;
  
}  // end calc_var_ohlc




// The function roll_var_rcpp() calculates a vector of variance estimates over a
// rolling look-back interval for a vector, using Rcpp.
//' @export
// [[Rcpp::export]]
NumericVector roll_var_rcpp(NumericVector vectorv, int look_back=11) {
  int nrows = vectorv.size();
  NumericVector sub_vec(look_back);
  NumericVector var_vec(nrows);
  NumericVector roll_mean(nrows);
  
  var_vec[0] = 0;
  for (int it = 1; it < vectorv.size(); it++) {
    sub_vec = vectorv[Range(std::max(it-look_back+1, 0), it)];
    var_vec[it] = variance_rcpp(sub_vec);
  }  // end for
  
  return var_vec;
  
}  // end roll_var_rcpp


// The function roll_var_vec() calculates a vector of variance estimates over a
// rolling look-back interval for a vector, using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::vec roll_var_vec(arma::vec& vectorv, arma::uword look_back=11) {
  arma::uword nrows = vectorv.n_elem;
  arma::vec var_vec = arma::zeros(nrows);

  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    var_vec(it) = arma::var(vectorv.subvec(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < nrows; it++) {
    var_vec(it) = arma::var(vectorv.subvec(it-look_back+1, it));
  }  // end for
  
  return var_vec;
  
}  // end roll_var_vec


//' @export
// [[Rcpp::export]]
arma::mat roll_var_mat(arma::mat& matrixv, arma::uword look_back=11) {
  arma::uword nrows = matrixv.n_rows;
  arma::mat var_mat = arma::zeros(nrows, matrixv.n_cols);
  
  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    var_mat.row(it) = arma::var(matrixv.rows(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < nrows; it++) {
    var_mat.row(it) = arma::var(matrixv.rows(it-look_back+1, it));
  }  // end for
  
  return var_mat;
  
}  // end roll_var_mat


//' @export
// [[Rcpp::export]]
arma::vec roll_var_ohlc(arma::mat& ohlc, 
                        std::string calc_method="yang_zhang", 
                        arma::colvec indeks=0, 
                        bool scalit=true, 
                        arma::uword look_back=11) {
  
  arma::uword nrows = ohlc.n_rows;
  arma::vec var_vec = arma::zeros(nrows);
  arma::colvec clo_se = ohlc.col(3);
  arma::colvec lag_close = lagit(clo_se);
  // std::string method = calc_method;
  // bool do_scale = scalit;

  if (!scalit || (indeks.n_rows == 1)) {
    indeks = arma::ones(nrows);
    // cout << "ohlc.n_rows = " << nrows << endl;
    // cout << "indeks.n_rows = " << indeks.n_rows << endl;
  }  // end if
  // 
  // cout << "indeks.n_rows = " << indeks.n_rows << endl;
  // cout << "indeks.rows(0, 3) = " << indeks.rows(0, 3) << endl;
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    cout << "it = " << it << endl;
    // cout << "lag_close = " << lag_close.subvec(0, it) << endl;
    // cout << "indeks = " << indeks.subvec(0, it) << endl;
    // var_vec(it) = sum(indeks.subvec(0, it));
    arma::mat sub_ohlc = ohlc.rows(0, it);
    arma::colvec sub_close = lag_close.rows(0, it);
    arma::colvec sub_index = indeks.subvec(0, it);
    // cout << "sum(sub_ohlc) = " << sum(sub_ohlc) << endl;
    // cout << "sum(sub_close) = " << sum(sub_close) << endl;
    // cout << "sum(sub_index) = " << sum(sub_index) << endl;
    // var_vec(it) = it;
    // var_vec(it) = calc_var_ohlc(ohlc=ohlc, lag_close=lag_close, calc_method=calc_method, indeks=indeks, scalit=scalit);
    var_vec(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scalit);
    // var_vec(it) = calc_var_ohlc(ohlc=ohlc.rows(0, it), lag_close=lag_close.rows(0, it), calc_method=calc_method, indeks=indeks.subvec(0, it), scalit=scalit);
    // arma::mat sub_ohlc = ohlc.rows(0, it);
    // var_vec(it) = calc_matp(sub_ohlc, sub_close, calc_method, sub_index, scalit);
    // cout << "it = " << it << endl;
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < nrows; it++) {
    // cout << "it = " << it << endl;
    arma::mat sub_ohlc = ohlc.rows(it-look_back+1, it);
    arma::colvec sub_close = lag_close.rows(it-look_back+1, it);
    arma::colvec sub_index = indeks.subvec(it-look_back+1, it);
    var_vec(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scalit);
    // var_vec(it) = sum(indeks.subvec(it-look_back+1, it));
  }  // end for
  
  return var_vec;
  
}  // end roll_var_ohlc


