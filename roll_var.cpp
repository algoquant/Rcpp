// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


//' @export
// [[Rcpp::export]]
arma::mat lag_it(arma::mat& mat_rix, int lagg=1) {
  int num_rows = (mat_rix.n_rows-1);
  
  if (lagg > 0)
    return arma::join_cols(arma::repelem(mat_rix.row(0), lagg, 1), 
                           mat_rix.rows(0, num_rows-lagg));
  else
    return arma::join_cols(mat_rix.rows(-lagg, num_rows), 
                           arma::repelem(mat_rix.row(num_rows), -lagg, 1));
  
}  // end lag_it



// The function variance_rcpp() calculates the variance of a vector using Rcpp.
//' @export
// [[Rcpp::export]]
double variance_rcpp(NumericVector& vec_tor) {
  return sum(pow(vec_tor - sum(vec_tor)/vec_tor.size(), 2))/(vec_tor.size()-1);
}  // end variance_rcpp


// The function variance_arma() calculates the variance of a vector using 
// RcppArmadillo.
//' @export
// [[Rcpp::export]]
double variance_arma(arma::vec& vec_tor) {
  return arma::var(vec_tor);
}  // end variance_arma


//' @export
// [[Rcpp::export]]
double vec_prod(const arma::colvec& vec1, const arma::colvec& vec2) {
  return arma::dot(vec1, vec2);
}  // end vec_prod


//' @export
// [[Rcpp::export]]
double calc_mat(arma::mat oh_lc) {
  return arma::sum(arma::var(oh_lc));
}  // end calc_mat


//' @export
// [[Rcpp::export]]
double calc_matp(arma::mat& oh_lc, 
                 arma::colvec lag_close=0, 
                 std::string calc_method="yang_zhang", 
                 arma::colvec in_dex=0, 
                 bool scal_e=true) {
  return arma::sum(arma::var(oh_lc));
}  // end calc_matp


//' @export
// [[Rcpp::export]]
double variance_ohlc(arma::mat& oh_lc) {
  int num_rows = oh_lc.n_rows;
  double co_eff = 0.34/(1.34 + (num_rows+1)/(num_rows-1));
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec open_close = arma::join_cols(clo_se.subvec(0, 0), clo_se.subvec(0, clo_se.n_elem-2));
  open_close = oh_lc.col(0) - open_close;
  
  return arma::var(open_close) + co_eff*arma::var(oh_lc.col(0) - clo_se) +
    (1 - co_eff)*(arma::dot((oh_lc.col(1) - clo_se), (oh_lc.col(1) - oh_lc.col(0))) +
    arma::dot((oh_lc.col(2) - clo_se), (oh_lc.col(2) - oh_lc.col(0))))/num_rows;
  
}  // end variance_ohlc


// The below function with an explicit loop instead of a lag is actually 
// a bit slower than the above function with a lag
// double variance_ohlc_plus(arma::mat& oh_lc) {
//   int num_rows = oh_lc.n_rows;
//   double co_eff = 0.34/(1.34 + (num_rows+1)/(num_rows-1));
//   arma::colvec op_en = oh_lc.col(0);
//   arma::colvec clo_se = oh_lc.col(3);
//   arma::colvec open_close(num_rows);
//   // open_close = oh_lc.col(0) - open_close;
// 
//   // So instead copy using an explicit loop:
//   for (int it = 1; it < num_rows; it++) {
//     open_close(it) = (op_en(it) - clo_se(it-1));
//   }  // end for
//   open_close(0) = (op_en(0) - clo_se(0));
//   
//   return arma::var(open_close) + co_eff*arma::var(oh_lc.col(0) - clo_se) +
//     (1 - co_eff)*(arma::dot((oh_lc.col(1) - clo_se), (oh_lc.col(1) - oh_lc.col(0))) +
//     arma::dot((oh_lc.col(2) - clo_se), (oh_lc.col(2) - oh_lc.col(0))))/num_rows;
//   
// }  // end variance_ohlc_plus



//' @export
// [[Rcpp::export]]
double calc_var_ohlc(arma::mat& oh_lc, 
                     std::string calc_method="yang_zhang", 
                     arma::colvec lag_close=0, 
                     arma::colvec in_dex=0, 
                     bool scal_e=true) {
  
  int num_rows = oh_lc.n_rows;
  double co_eff = 0.34/(1.34 + (num_rows+1)/(num_rows-1));

  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
    // cout << "oh_lc.n_rows = " << num_rows << endl;
    // cout << "in_dex.n_rows = " << in_dex.n_rows << endl;
  }  // end if
  
  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec open_close(clo_se.n_rows);
  if (lag_close.n_rows == 1) {
    open_close = arma::join_cols(clo_se.subvec(0, 0), clo_se.subvec(0, clo_se.n_elem-2));
    open_close = (oh_lc.col(0) - open_close)/in_dex;
  } else {
    open_close = (oh_lc.col(0) - lag_close)/in_dex;
  }  // end if
  arma::colvec close_open = (clo_se - oh_lc.col(0))/in_dex;
  arma::colvec close_high = (clo_se - oh_lc.col(1))/in_dex;
  arma::colvec close_low = (clo_se - oh_lc.col(2))/in_dex;
  arma::colvec high_low = (oh_lc.col(1) - oh_lc.col(2))/in_dex;
  arma::colvec high_open = (oh_lc.col(1) - oh_lc.col(0))/in_dex;
  arma::colvec low_open = (oh_lc.col(2) - oh_lc.col(0))/in_dex;
  
  cout << "oh_lc.n_rows = " << num_rows << endl;
  // cout << "oh_lc.n_cols = " << oh_lc.n_cols << endl;
  cout << "open_close.n_rows = " << open_close.n_rows << endl;
  // cout << "open_close.n_cols = " << open_close.n_cols << endl;
  // cout << "in_dex.n_rows = " << in_dex.n_rows << endl;
  // cout << "in_dex.n_cols = " << in_dex.n_cols << endl;
  // cout << "oh_lc last row = " << oh_lc.row(num_rows-1) << endl;
  // cout << "sum clo_se = " << arma::var(arma::diff(clo_se)) << endl;

  if (calc_method == "close") {
    // cout << "Calc method is Close" << endl;
    return arma::var(arma::diff(clo_se));
  } else if (calc_method == "rogers_satchell") {
    // cout << "Calc method is Rogers-Satchell" << endl;
    return -(arma::dot(close_high, high_open) +
             arma::dot(close_low, low_open))/num_rows;
  } else if (calc_method == "garman_klass") {
    // cout << "Calc method is Garman-Klass" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows;
  } else if (calc_method == "garman_klass_yz") {
    // cout << "Calc method is Garman-Klass-YZ" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows + 
            arma::var(open_close);
  } else if (calc_method == "yang_zhang") {
    // cout << "Calc method is Yang-Zhang" << endl;
    return arma::var(open_close) + co_eff*arma::var(close_open) +
      (co_eff-1)*(arma::dot(close_high, high_open) + 
      arma::dot(close_low, low_open))/num_rows;
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
NumericVector roll_var_rcpp(NumericVector vec_tor, int look_back=11) {
  int len_gth = vec_tor.size();
  NumericVector sub_vec(look_back);
  NumericVector var_vec(len_gth);
  NumericVector roll_mean(len_gth);
  
  var_vec[0] = 0;
  for (int it = 1; it < vec_tor.size(); it++) {
    sub_vec = vec_tor[Range(std::max(it-look_back+1, 0), it)];
    var_vec[it] = variance_rcpp(sub_vec);
  }  // end for
  
  return var_vec;
  
}  // end roll_var_rcpp


// The function roll_var_vec() calculates a vector of variance estimates over a
// rolling look-back interval for a vector, using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::vec roll_var_vec(arma::vec& vec_tor, arma::uword look_back=11) {
  arma::uword len_gth = vec_tor.n_elem;
  arma::vec var_vec = arma::zeros(len_gth);

  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    var_vec(it) = arma::var(vec_tor.subvec(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < len_gth; it++) {
    var_vec(it) = arma::var(vec_tor.subvec(it-look_back+1, it));
  }  // end for
  
  return var_vec;
  
}  // end roll_var_vec


//' @export
// [[Rcpp::export]]
arma::mat roll_var_mat(arma::mat& mat_rix, arma::uword look_back=11) {
  arma::uword num_rows = mat_rix.n_rows;
  arma::mat var_mat = arma::zeros(num_rows, mat_rix.n_cols);
  
  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    var_mat.row(it) = arma::var(mat_rix.rows(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < num_rows; it++) {
    var_mat.row(it) = arma::var(mat_rix.rows(it-look_back+1, it));
  }  // end for
  
  return var_mat;
  
}  // end roll_var_mat


//' @export
// [[Rcpp::export]]
arma::vec roll_var_ohlc(arma::mat& oh_lc, 
                        std::string calc_method="yang_zhang", 
                        arma::colvec in_dex=0, 
                        bool scal_e=true, 
                        arma::uword look_back=11) {
  
  arma::uword num_rows = oh_lc.n_rows;
  arma::vec var_vec = arma::zeros(num_rows);
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec lag_close = lag_it(clo_se);
  // std::string meth_od = calc_method;
  // bool do_scale = scal_e;

  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
    // cout << "oh_lc.n_rows = " << num_rows << endl;
    // cout << "in_dex.n_rows = " << in_dex.n_rows << endl;
  }  // end if
  // 
  // cout << "in_dex.n_rows = " << in_dex.n_rows << endl;
  // cout << "in_dex.rows(0, 3) = " << in_dex.rows(0, 3) << endl;
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    cout << "it = " << it << endl;
    // cout << "lag_close = " << lag_close.subvec(0, it) << endl;
    // cout << "in_dex = " << in_dex.subvec(0, it) << endl;
    // var_vec(it) = sum(in_dex.subvec(0, it));
    arma::mat sub_ohlc = oh_lc.rows(0, it);
    arma::colvec sub_close = lag_close.rows(0, it);
    arma::colvec sub_index = in_dex.subvec(0, it);
    // cout << "sum(sub_ohlc) = " << sum(sub_ohlc) << endl;
    // cout << "sum(sub_close) = " << sum(sub_close) << endl;
    // cout << "sum(sub_index) = " << sum(sub_index) << endl;
    // var_vec(it) = it;
    // var_vec(it) = calc_var_ohlc(oh_lc=oh_lc, lag_close=lag_close, calc_method=calc_method, in_dex=in_dex, scal_e=scal_e);
    var_vec(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
    // var_vec(it) = calc_var_ohlc(oh_lc=oh_lc.rows(0, it), lag_close=lag_close.rows(0, it), calc_method=calc_method, in_dex=in_dex.subvec(0, it), scal_e=scal_e);
    // arma::mat sub_ohlc = oh_lc.rows(0, it);
    // var_vec(it) = calc_matp(sub_ohlc, sub_close, calc_method, sub_index, scal_e);
    // cout << "it = " << it << endl;
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < num_rows; it++) {
    // cout << "it = " << it << endl;
    arma::mat sub_ohlc = oh_lc.rows(it-look_back+1, it);
    arma::colvec sub_close = lag_close.rows(it-look_back+1, it);
    arma::colvec sub_index = in_dex.subvec(it-look_back+1, it);
    var_vec(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
    // var_vec(it) = sum(in_dex.subvec(it-look_back+1, it));
  }  // end for
  
  return var_vec;
  
}  // end roll_var_ohlc


