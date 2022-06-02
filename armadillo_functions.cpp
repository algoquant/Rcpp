// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/armadillo_functions.cpp")

////////////////////////////
// RcppArmadillo functions miscellaneous
////////////////////////////

// The function which_method() switches through options
//' @export
// [[Rcpp::export]]
void which_method(const std::string& calc_method = "yang_zhang") {
  
  // cout << "Calc method is " << calc_method << endl;
  
  if (calc_method == "close") 
    cout << "Calc method is Close" << endl;
  else if (calc_method == "garman_klass") 
    cout << "Calc method is Garman Klass" << endl;
  else if (calc_method == "garman_klass_yz") 
    cout << "Calc method is Garman Klass YZ" << endl;
  else if (calc_method == "yang_zhang") 
    cout << "Calc method is Yang Zhang" << endl;
  else 
    cout << "Wrong calc method!" << endl;
  
  // return "Done";
}  // end which_method



// The function countfun() counts the number of times it was called.
// It creates a static integer variable countv, and advances it every 
// time it is called.
// countv is static so it remains alive outside the scope of countfun(),
// between calls to countfun().
// The function countfun() can also be defined as a static function as 
// an illustration - it doesn't have to be static.
// A static function in C is only visible to those functions in the same source file. 
// Making a function static limits its scope to functions from the same source file. 
// https://www.cprogramming.com/tutorial/statickeyword.html
// [[Rcpp::export]]
void countfun(int initv=1) {
  
  // countv is static so it's initialized only once the first time countfun() is called.
  static int countv = initv;
  
  countv++;
  
  std::cout << "The function countfun() was called " << countv << " times." << std::endl;
  
}  // end countfun



// First sort the datav and then unsort it back to original
//' @export
// [[Rcpp::export]]
arma::vec sort_back(const arma::vec& datav) {
  
  // Reverse sort index
  arma::uvec indeks = arma::sort_index(arma::sort_index(datav));
  // Sort the datav
  arma::vec sortedv = arma::sort(datav);
  // Reverse sort the datav
  sortedv = sortedv.elem(indeks);
  
  return sortedv;
}  // end sort_back




// Calculates the ranks of a vector of data
//' @export
// [[Rcpp::export]]
arma::uvec calc_ranks(const arma::vec& datav) {
  return (arma::sort_index(arma::sort_index(datav)) + 1);
}  // end calc_ranks



// Calculates the de-meaned ranks of a vector of data
//' @export
// [[Rcpp::export]]
arma::vec calc_ranksm(const arma::vec& datav) {
  // Ranks
  arma::vec ranks = conv_to< vec >::from(arma::sort_index(arma::sort_index(datav)));
  return (ranks - arma::mean(ranks));
}  // end calc_ranksm




////////////////////////////
// RcppArmadillo functions for manipulating vectors and matrices
////////////////////////////

// Conditional sums using Rcpp

// The function sum_na() performs sum(is.na(vectorv)) using Rcpp.
// Adapted from:
// https://stackoverflow.com/questions/46892399/fast-checking-of-missing-values-in-rcpp
//' @export
// [[Rcpp::export]]
int sum_na(const NumericVector& vectorv) {
  int count_ = 0;
  for (int i = 0; i < vectorv.size(); i++)
    count_ += NumericVector::is_na(vectorv[i]);
  return count_;
}  // end sum_na


// The function sum_na_stl() performs sum(is.na(vectorv)) using Rcpp and STL.
//' @export
// [[Rcpp::export]]
int sum_na_stl(const NumericVector& vectorv) {
  return std::count_if(vectorv.begin(), vectorv.end(), NumericVector::is_na);
}  // end sum_na_stl


// The function sum_if_cpp() performs sum(ivectorv<findv) using Rcpp.
//' @export
// [[Rcpp::export]]
int sum_if_cpp(NumericVector& vectorv, double findv) {
  int count_ = 0;
  for (int i = 0; i < vectorv.size(); i++)
    count_ += (vectorv[i]<findv);
  return count_;
}  // end sum_if_cpp


// The function sum_if() performs sum(ivectorv<findv) using Rcpp.
//' @export
// [[Rcpp::export]]
int sum_if(NumericVector& vectorv, double findv) {
  int count_ = 0;
  NumericVector::iterator it;
  for (it = vectorv.begin(); it != vectorv.end(); it++) {
    if (*it < findv) count_++;
  }  // end for
  return count_;
}  // end sum_if


// The function sum_if_stl() performs sum(ivectorv<findv) using Rcpp and STL.
// Warning: bind2nd is deprecated 
//' @export
// [[Rcpp::export]]
int sum_if_stl(const NumericVector& vectorv, double findv) {
  return std::count_if(vectorv.begin(), vectorv.end(), bind2nd(std::less<double>(), findv));
}  // end sum_if_stl


// The function whichv() performs which() using Rcpp with Sugar and RcppArmadillo
//' @export
// [[Rcpp::export]]
arma::uvec whichv(arma::uvec logicalv) {
  // arma::uvec logicalv;
  // logicalv = arma::find(logicalv);
  return arma::find(logicalv);
}  // end whichv


// The function whichv2() performs which() using Rcpp and arma::find() from
// RcppArmadillo
//' @export
// [[Rcpp::export]]
arma::uvec whichv2(LogicalVector& vectorv) {
  arma::uvec logicalv;
  logicalv = arma::find(Rcpp::as<uvec>(vectorv));
  return logicalv+1;
}  // end whichv2


// The function whichv3() performs which() using Rcpp and STL
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector whichv3(const LogicalVector& vectorv) {
  // arma::uvec logicalv;
  // logicalv = arma::find(Rcpp::as<uvec>(vectorv));
  int nrow = vectorv.size();
  std::vector<int> indeks;
  indeks.reserve(nrow);
  for (int i=0; i<nrow; i++) {
    if (vectorv[i]) indeks.push_back(i+1);
  }  // end for
  return Rcpp::wrap(indeks);
  // return logicalv;
}  // end whichv3


// The function whichv32() performs which() using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector whichv32(const LogicalVector& vectorv) {
  int nrow = vectorv.size();
  IntegerVector indeks(sum(vectorv));
  int j = 0;
  for (int i=0; i<nrow; i++) {
    if (vectorv[i]) 
      indeks(j++) = i+1;
  }  // end for
  return indeks;
}  // end whichv32


// The function whichv34() performs which() using Rcpp.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector whichv34(LogicalVector& vectorv) {
  // int nrow = vectorv.size();
  IntegerVector indeks(sum(vectorv));
  int i=0, j=0;
  LogicalVector::iterator it;
  for (it = vectorv.begin(); it != vectorv.end(); it++, i++) {
    if (*it) 
      indeks(j++) = i+1;
  }  // end for
  return indeks;
}  // end whichv34


// The function whichv35() performs which() using Rcpp and STL.
// Currently it doesn't work.
// int whichv35(LogicalVector& vectorv) {
//   return std::find(std::begin(vectorv), std::end(vectorv), TRUE);
// }  // end whichv35


// The function whichv33() performs which() using RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::uvec whichv33(arma::uvec& vectorv) {
  int nrow = vectorv.size();
  arma::uvec indeks(accu(vectorv));
  int j = 0;
  for (int i=0; i<nrow; i++) {
    if (vectorv[i]) 
      indeks(j++) = i+1;
  }  // end for
  return indeks;
}  // end whichv33


// The function whichv4() performs which() using Rcpp
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector whichv4(LogicalVector& vectorv) {
  // arma::uvec logicalv;
  // logicalv = arma::find(Rcpp::as<uvec>(vectorv));
  int nrow = vectorv.size();
  IntegerVector indeks(nrow);
  int j=0;
  for (int i=0; i<nrow; i++) {
    if (vectorv[i]) {
      indeks(j)=i+1;
      j++;
    }  // end if
  }  // end for
  vector<int> sub_vector(indeks.begin(), indeks.begin() + j);
  return wrap(sub_vector);
  // return logicalv;
}  // end whichv4


// The function whichv5() calls R function which().
// It's very slow.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector whichv5(LogicalVector& vectorv) {
  // Obtain environment containing the function
  Rcpp::Environment base("package:base"); 
  // Define function which() callable from Rcpp
  Rcpp::Function which = base["which"];
  return which(vectorv);
  // return logicalv;
}  // end whichv5


// The function tapply_arma() performs aggregations over a vector using a
// factor.
// It produces the same result as the R code: 
//    tapply(X=vectorv, INDEX=ratio, FUN=NROW)
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::vec tapply_arma(const arma::vec& vectorv, const arma::vec& ratio) {
  Function whichv3("whichv3");
  // int nrow = vectorv.size();
  arma::vec uniq_ue = arma::unique(ratio);
  int n_unique = uniq_ue.size();
  arma::vec agg_s(n_unique);
  arma::uvec indeks;
  arma::vec sub_vec;
  for (int i=0; i<n_unique; i++) {
    indeks = find(ratio == uniq_ue(i));
    sub_vec = vectorv(indeks);
    agg_s(i) = sub_vec.n_elem;
  }  // end for
  return agg_s;
  // return find_unique(vectorv);
}  // end tapply_arma



// The function cbind_rcpp() binds the columns of two matrices.
// It uses Rcpp.
//' @export
// [[Rcpp::export]]
NumericMatrix cbind_rcpp(NumericMatrix matrixv1, NumericMatrix matrixv2) {
  int ncol1 = matrixv1.ncol();
  int ncol2 = matrixv2.ncol();
  NumericMatrix output = Rcpp::no_init_matrix(matrixv1.nrow(), ncol1 + ncol2);
  for (int j = 0; j < ncol1 + ncol2; j++) {
    if (j < ncol1) {
      output(_, j) = matrixv1(_, j);
    } else {
      output(_, j) = matrixv2(_, j - ncol1);
    }
  }
  return output;
}  // end cbind_rcpp


// The function cbind_arma() binds the columns of two matrices.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::mat cbind_arma(const arma::mat& matrixv1, const arma::mat& matrixv2) {
  return arma::join_rows(matrixv1, matrixv2);
}  // end cbind_arma


// The function sub_mat() selects a sub-matrix using RcppArmadillo
//' @export
// [[Rcpp::export]]
arma::mat sub_mat(const arma::mat& matrixv, arma::uvec row_num, arma::uvec col_num) {
  // arma::mat sub_matrix = matrixv.cols(col_num);
  arma::mat sub_matrix = matrixv.submat(row_num, col_num);
  return sub_matrix;
}  // end sub_mat


// The function sub_mat_cast() selects a sub-matrix using RcppArmadillo with Rcpp type arguments
// It casts Rcpp arguments into RcppArmadillo structures
//' @export
// [[Rcpp::export]]
arma::mat sub_mat_cast(NumericMatrix& matrixv, IntegerVector& row_num, IntegerVector col_num=0) {
  // arma::mat sub_matrix = matrixv.cols(row_num);
  // arma::mat sub_matrix = matrixv.cols(col_num);
  arma::mat sub_matrix = as<mat>(matrixv);
  // arma::uvec urow_num = as<uvec>(row_num);
  // arma::uvec ucol_num = as<uvec>(col_num);
  // sub_matrix = sub_matrix.submat(as<uvec>(row_num), as<uvec>(col_num));
  sub_matrix = sub_matrix.submat(as<uvec>(row_num), as<uvec>(col_num));
  return sub_matrix;
  // return Rcpp::wrap(sub_matrix);
}  // end sub_mat_cast


// The function sub_assign() assigns values to selected vector elements using
// RcppArmadillo.
// Adapted from:
// http://gallery.rcpp.org/articles/armadillo-subsetting/index.html
//' @export
// [[Rcpp::export]]
arma::rowvec sub_assign(arma::rowvec vectorv, arma::uvec indeks, arma::vec datav) {
  vectorv.elem(indeks) = datav;
  return vectorv;
}  // end sub_assign


// The function find_assign_vec() assigns values to selected vector elements
// using RcppArmadillo.
// Adapted from:
// http://gallery.rcpp.org/articles/armadillo-subsetting/index.html
//' @export
// [[Rcpp::export]]
arma::rowvec find_assign_vec(arma::rowvec& vectorv, double findv, double datav) {
  arma::uvec indeks = find(vectorv == findv);
  vectorv(indeks).fill(datav);
  return vectorv;
}  // end find_assign_vec


// The function find_assign_vec_point() assigns values to selected vector
// elements using RcppArmadillo.
// It accepts a pointer to the vector, assigns values in place without copying
// the input vector, and returns an integer with the number of assignments
// performed. 
// Adapted from:
// http://gallery.rcpp.org/articles/armadillo-subsetting/index.html
//' @export
// [[Rcpp::export]]
double find_assign_vec_point(arma::rowvec& vectorv, double findv, double datav) {
  arma::uvec indeks = find(vectorv > findv);
  vectorv(indeks).fill(datav);
  return arma::accu(indeks);
}  // end find_assign_vec_point


// The function find_assign_mat() assigns values to selected matrix elements
// using RcppArmadillo.
// Adapted from:
// http://gallery.rcpp.org/articles/armadillo-subsetting/index.html
//' @export
// [[Rcpp::export]]
arma::mat find_assign_mat(arma::mat matrixv, double findv, double datav) {
  arma::uvec indeks = find(matrixv >= findv);
  matrixv.elem(indeks).fill(datav);
  return matrixv;
}  // end find_assign_mat


// The function compare_col() selects a matrix column and compares it 
// to an input value, using Rcpp with Sugar.
//' @export
// [[Rcpp::export]]
LogicalVector compare_col(NumericMatrix& matrixv, double findv, int col_num=0) {
  // NumericVector colnum = matrixv(_, col_num);
  // LogicalVector bar;
  // bar = (colnum > findv);
  return (matrixv(_, col_num) > findv);
}  // end compare_col


// The function compare_col_arma() selects a matrix column and compares 
// it to an input value, using Rcpp with Sugar and RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::uvec compare_col_arma(NumericMatrix& matrixv, double findv, int col_num=0) {
  NumericVector colnum = matrixv(_, col_num);
  LogicalVector bar;
  bar = (colnum > findv);
  return Rcpp::as<uvec>(bar);
}  // end compare_col_arma


// The function compare_col_armaa() selects a matrix column and compares 
// it to an input value, using Rcpp with Sugar and RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::uvec compare_col_armaa(const arma::mat& matrixv, double findv, int col_num=0) {
  arma::vec colnum = matrixv.cols(as<uvec>(wrap(col_num)));
  // NumericVector colnum = matrixv(_, col_num);
  // LogicalVector bar;
  // arma::uvec bar;
  // bar = (colnum > findv);
  // for (int i=0; i<colnum.n_rows; i++) {
  //   bar[i] = (colnum[i] > findv);
  // }  // end for
  return (colnum > findv);
}  // end compare_col_armaa


// The function which_col() selects a matrix column and performs 
// which() on it.
// It uses Rcpp with Sugar and RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::uvec which_col(NumericMatrix& matrixv, double findv, int col_num=0) {
  NumericVector colnum = matrixv(_, col_num);
  LogicalVector vectorv = (colnum > findv);
  arma::uvec whichv;
  whichv = arma::find(Rcpp::as<uvec>(vectorv));
  return whichv;
}  // end which_col


// The function find_extract_mat() extracts matrix elements using RcppArmadillo.
// Adapted from:
// http://gallery.rcpp.org/articles/armadillo-subsetting/index.html
//' @export
// [[Rcpp::export]]
arma::vec find_extract_mat(arma::mat matrixv, double findv) {
  arma::uvec indeks = find(matrixv >= findv);
  return matrixv.elem(indeks);
}  // end find_extract_mat


// The function find_sub_mat() selects a matrix column, performs 
// find() on it, and then selects the corresponding matrix rows.
// It uses Rcpp with Sugar and RcppArmadillo function find().
//' @export
// [[Rcpp::export]]
NumericMatrix find_sub_mat(NumericMatrix& matrixv, double findv, int col_num=0) {
  NumericVector colnum = matrixv(_, col_num);
  LogicalVector vectorv = (colnum > findv);
  arma::uvec whichv;
  whichv = arma::find(Rcpp::as<uvec>(vectorv));
  arma::mat sub_matrix = as<mat>(matrixv);
  sub_matrix = sub_matrix.rows(whichv);
  return wrap(sub_matrix);
}  // end find_sub_mat


// The function select_sub_mat() selects a matrix column, performs 
// which on it, and then selects the corresponding matrix rows.
// It uses Rcpp with Sugar and RcppArmadillo, but without function find().
//' @export
// [[Rcpp::export]]
NumericMatrix select_sub_mat(NumericMatrix& matrixv, double findv, int col_num=0) {
  NumericVector colnum = matrixv(_, col_num);
  LogicalVector vectorv = (colnum > findv);
  int nrow = vectorv.size();
  // perform which
  std::vector<int> indeks;
  indeks.reserve(nrow);
  for (int i=0; i<nrow; i++) {
    if (vectorv[i]) indeks.push_back(i);
  }  // end for
  // perform which
  arma::mat sub_matrix = as<mat>(matrixv);
  sub_matrix = sub_matrix.rows(as<uvec>(wrap(indeks)));
  return wrap(sub_matrix);
}  // end select_sub_mat


// The function apply_agg() performs aggregations over a matrix using its
// first column as a factor.
// It produces the same result as the R code: 
//      sapply(X=unique(matrixv[, 1]), FUN=function(matrixv[, -1]))
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::vec apply_agg(const arma::mat& matrixv) {
  // Function whichv3("whichv3");
  // int nrow = vectorv.size();
  arma::vec uniq_ue = arma::unique(matrixv.col(0));
  int n_unique = uniq_ue.size();
  arma::vec agg_s(n_unique);
  arma::uvec indeks;
  arma::mat sub_mat = matrixv.cols(1, matrixv.n_cols-1);
  for (int i=0; i<n_unique; i++) {
    indeks = find(matrixv.col(0) == uniq_ue(i));
    agg_s(i) = as_scalar(sum(prod(sub_mat.rows(indeks), 1)));
    // agg_s(i) = as_scalar(accu(sub_mat.rows(indeks)));
  }  // end for
  return agg_s;
  // return find_unique(matrixv);
}  // end apply_agg


//' @export
// [[Rcpp::export]]
double  agg_mat(const arma::mat& matrixv) {
  return  arma::as_scalar(sum(prod(matrixv, 1)));
  // return find_unique(matrixv);
}  // end agg_mat


// The function variance() calculates the variance of a vector using Rcpp
//' @export
// [[Rcpp::export]]
double variance(NumericVector vectorv) {
  return sum(pow(vectorv - sum(vectorv)/vectorv.size(), 2));
}  // end variance



// The function roll_rows() rolls over the rows of a matrix
//' @export
// [[Rcpp::export]]
arma::vec roll_rows(arma::mat& ohlc, 
                    arma::uword look_back=11) {
  
  arma::uword nrows = ohlc.n_rows;
  arma::vec var_vec = arma::zeros(nrows);
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    arma::mat sub_ohlc = ohlc.rows(0, it);
    var_vec(it) = sum(sum(sub_ohlc));
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < nrows; it++) {
    arma::mat sub_ohlc = ohlc.rows(it-look_back+1, it);
    var_vec(it) = sum(sum(sub_ohlc));
  }  // end for
  
  return var_vec;
  
}  // end roll_rows



// The function roll_cols() rolls over the columns of a matrix
//' @export
// [[Rcpp::export]]
arma::vec roll_cols(arma::mat& ohlc, 
                    arma::uword look_back=11) {
  
  arma::uword ncols = ohlc.n_cols;
  arma::vec var_vec = arma::zeros(ncols);
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    arma::mat sub_ohlc = ohlc.cols(0, it);
    var_vec(it) = sum(sum(sub_ohlc));
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < ncols; it++) {
    arma::mat sub_ohlc = ohlc.cols(it-look_back+1, it);
    var_vec(it) = sum(sum(sub_ohlc));
  }  // end for
  
  return var_vec;
  
}  // end roll_cols




////////////////////////////
// RcppArmadillo functions for matrix algebra
////////////////////////////


// The function demean_arma() calculates a matrix with de-meaned columns.
// It accepts a pointer to a matrix and returns a copy of the de-meaned matrix.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::mat demean_arma(const arma::mat& matrixv) {
  // de-mean response and explanatory variables
  arma::mat mat_demean(matrixv.n_rows, matrixv.n_cols);
  for (uword i = 0; i < matrixv.n_cols; i++) {
    mat_demean.col(i) = matrixv.col(i) - arma::mean(matrixv.col(i));
    // mat_demean(i) = arma::mean(matrixv.col(i));
  }  // end for
  return mat_demean;
}  // end demean_arma


// The function demean_mat() calculates a matrix with de-meaned columns.
// It accepts a pointer to a matrix and operates on the matrix in place.
// It returns the number of columns of the input matrix.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
int demean_mat(arma::mat& matrixv) {
  // de-mean response and explanatory variables
  // arma::mat mat_demean(matrixv.n_cols);
  for (uword i = 0; i < matrixv.n_cols; i++) {
    matrixv.col(i) -= arma::mean(matrixv.col(i));
    // mat_demean(i) = arma::mean(matrixv.col(i));
  }  // end for
  return matrixv.n_cols;
}  // end demean_mat


// The function outer_vec() calculates the outer product of two vectors.
// It accepts pointers to the two vectors and returns a matrix.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::mat outer_vec(const arma::vec& vec1, const arma::vec& vec2) {
  return vec1*vec2.t();
}  // end outer_vec


// The function inner_vec() calculates the inner (dot) product of two vectors.
// It accepts pointers to the two vectors and returns a double.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
double inner_vec(const arma::vec& vec1, const arma::vec& vec2) {
  return arma::dot(vec1, vec2);
}  // end inner_vec


// The function mat_inner_vec() calculates the product of a matrix times a vector.
// It accepts pointers to the matrix and vector, and returns a vector.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::vec mat_inner_vec(const arma::vec& vectorv, const arma::mat& matrixv) {
  return matrixv * vectorv;
}  // end mat_inner_vec


// The function inner_mat() calculates the inner (dot) product of a matrix
// with two vectors.
// It accepts pointers to the matrix and vectors, and returns a double.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
double inner_mat(const arma::vec& vectorv2, const arma::mat& matrixv, const arma::vec& vectorv1) {
  return arma::as_scalar(trans(vectorv2) * (matrixv * vectorv1));
}  // end inner_mat


// [[Rcpp::export]]
void mult_rows(arma::vec& weightv, arma::mat& matrixv) {
  
  matrixv.each_row() %= weightv.t();
  
  // matrixv.each_row() %= weightv.t();
  // return matrixv;
  
}  // end mult_rows



// Already in HighFreq.cpp
////////////////////////////////////////////////////////////
//' Multiply the rows or columns of a \emph{matrix} times a \emph{vector},
//' element-wise and in place (without copying).
//' 
//' @param \code{vector} A \emph{numeric} \emph{vector}.
//' 
//' @param \code{matrix} A \emph{numeric} \emph{matrix}.
//' 
//' @param \code{byrow} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the rows of \code{matrix} by \code{vector}, otherwise multiply the columns
//'   (the default is \code{byrow = TRUE}.)
//' 
//' @return Void (no return value).
//' 
//' @details
//'   The function \code{mult_mat_ref()} multiplies the rows or columns of a
//'   \emph{matrix} times a \emph{vector}, element-wise and in place (without
//'   copying).
//'
//'   It accepts a \emph{pointer} to the argument \code{matrix}, and replaces
//'   the old \code{matrix} values with the new values. It performs the
//'   calculation in place, without copying the \emph{matrix} in memory, which
//'   can significantly increase the computation speed for large matrices.
//'
//'   If \code{byrow = TRUE} (the default), then function \code{mult_mat_ref()}
//'   multiplies the rows of the argument \code{matrix} times the argument
//'   \code{vector}.
//'   Otherwise it multiplies the columns of \code{matrix}.
//' 
//'   In \code{R}, \emph{matrix} multiplication is performed by columns.
//'   Performing multiplication by rows is often required, for example when
//'   multiplying stock returns by portfolio weights.
//'   But performing multiplication by rows requires explicit loops in \code{R},
//'   or it requires \emph{matrix} transpose.  And both are slow.
//'
//'   The function \code{mult_mat_ref()} uses \code{RcppArmadillo} \code{C++}
//'   code, so when multiplying large \emph{matrix} columns it's several times
//'   faster than vectorized \code{R} code, and it's even much faster compared
//'   to \code{R} when multiplying the \emph{matrix} rows.
//' 
//'   The function \code{mult_mat_ref()} performs loops over the \emph{matrix} rows
//'   and columns using the \emph{Armadillo} operators \code{each_row()} and
//'   \code{each_col()}, instead of performing explicit \code{for()} loops (both
//'   methods are equally fast).
//'   
//' @examples
//' \dontrun{
//' # Create vector and matrix data
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' 
//' # Multiply the matrix rows using R
//' matrixr <- t(vectorv*t(matrixv))
//' # Multiply the matrix rows using C++
//' HighFreq::mult_mat_ref(vectorv, matrixv, byrow=TRUE)
//' all.equal(matrixr, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat_ref(vectorv, matrixv, byrow=TRUE),
//'     Rcode=t(vectorv*t(matrixv)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//'     
//' # Multiply the matrix columns using R
//' matrixr <- vectorv*matrixv
//' # Multiply the matrix columns using C++
//' HighFreq::mult_mat_ref(vectorv, matrixv, byrow=FALSE)
//' all.equal(matrixr, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat_ref(vectorv, matrixv, byrow=FALSE),
//'     Rcode=vectorv*matrixv,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
void mult_mat_ref(arma::vec vector,
                  arma::mat matrix,
                  bool byrow = true) {
  
  arma::uword nelem = vector.n_elem;
  arma::uword nrows = matrix.n_rows;
  arma::uword ncols = matrix.n_cols;
  
  if (byrow && (nelem == ncols)) {
    // Multiply every row of matrix by vector
    matrix.each_row() %= vector.t();
  } else if (!byrow && (nelem == nrows)) {
    // Multiply every column of matrix by vector
    matrix.each_col() %= vector;
  } else {
    // Do nothing
    cout << "Nothing done: Vector length is neither equal to the number of columns nor to the rows of the matrix!" << std::endl;
  }  // end if
  
  // return matrix;
  
}  // end mult_mat_ref


////////////////////////////////////////////////////////////
//' Multiply the columns or rows of a matrix times a vector, element-wise.
//' 
//' @param vectorv A numeric \emph{vector}.
//' @param matrixv A numeric \emph{matrix}.
//' @param by_col A \emph{Boolean} argument: if \code{TRUE} then multiply the
//'   columns, else multiply the rows. (The default is \code{TRUE})
//' 
//' @return A single \emph{integer} value, equal to either the number of matrix
//'   columns or the number of rows.
//' 
//' @details The function \code{mult_vec_mat()} multiplies the columns or rows
//'   of a matrix times a vector, element-wise.
//'
//'   It accepts pointers to the matrix and vector, and performs the calculation
//'   in place, without copying the matrix in memory.
//'
//'   If the number of vector elements is equal to the number of matrix columns,
//'   then it multiplies the columns by the vector, and returns the number of
//'   columns. If the number of vector elements is equal to the number of rows,
//'   then it multiplies the rows, and returns the number of rows.
//'
//'   If the matrix is square and if \code{by_col} is \code{TRUE} then it
//'   multiplies the columns, else it multiplies the rows.
//'   
//'   The function \code{mult_vec_mat()} uses \emph{RcppArmadillo}, so when
//'   multiplying large matrix columns it's several times faster than vectorized
//'   \emph{R} code, and it's even much faster when multiply the matrix rows.
//'   
//' @examples
//' \dontrun{
//' # Multiply matrix columns
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' prod_uct <- vectorv*matrixv
//' mult_vec_mat(vectorv, matrixv)
//' all.equal(prod_uct, matrixv)
//' summary(microbenchmark(
//'     r_cpp=mult_vec_mat(vectorv, matrixv),
//'     r_code=vectorv*matrixv,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' 
//' # Multiply matrix rows
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' prod_uct <- t(vectorv*t(matrixv))
//' mult_vec_mat(vectorv, matrixv, by_col=FALSE)
//' all.equal(prod_uct, matrixv)
//' library(microbenchmark)
//' summary(microbenchmark(
//'     r_cpp=mult_vec_mat(vectorv, matrixv, by_col=FALSE),
//'     r_code=t(vectorv*t(matrixv)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
uword mult_vec_mat(const arma::vec& vectorv,
                   arma::mat& matrixv,
                   const bool& by_col=true) {
  uword nelem = vectorv.n_elem;
  uword nrows = matrixv.n_rows;
  uword ncols = matrixv.n_cols;
  
  if ((ncols == nrows) && (nelem == nrows)) {
    if (by_col) {
      // Multiply each column of matrixv by vectorv
      matrixv.each_col() %= vectorv;
      return nrows;
    } else {
      // Multiply each row of matrixv by vectorv
      matrixv.each_row() %= conv_to< rowvec >::from(vectorv);
      return ncols;
    }
  } else if (nelem == nrows) {
    // Multiply each column of matrixv by vectorv
    matrixv.each_col() %= vectorv;
    return nrows;
  } else if (nelem == ncols) {
    // Multiply each row of matrixv by vectorv
    matrixv.each_row() %= conv_to< rowvec >::from(vectorv);
    return ncols;
  } else 
    return NA_INTEGER;
}  // end mult_vec_mat


// Below are other examples of functions which multiply a matrix 
// times a vector, but are less efficient and/or less compact.

// The function mult_vec_mat2() calculates the product of matrix rows times a
// vector. It copies the vector into matrix rows so it's less efficient.
// It uses RcppArmadillo.
// 
//' @export
// [[Rcpp::export]]
int mult_vec_mat2(arma::rowvec& vectorv, 
                  arma::mat& matrixv) {
  // unsigned int mnrows = matrixv.n_rows;
  uword nrows = matrixv.n_rows;
  arma::mat vec_mat = repmat(vectorv, nrows, 1);
  matrixv = matrixv % vec_mat;
  return matrixv.n_cols;
}  // end mult_vec_mat2



// The function mult_vec_mat3() calculates the product of a matrix times a vector.
// It accepts pointers to the matrix and vector, and returns an integer.
// It accepts a pointer to a matrix and operates on the matrix in place.
// It performs a for() loop explicitly and is about as fast as using the
// implicit loops in each_col() and each_row().
// It uses RcppArmadillo.
// It returns the number of columns of the input matrix.
// 
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat3(arma::vec& vectorv, 
                          arma::mat& matrixv) {
  for (uword it = 0; it < matrixv.n_cols; it++) {
    matrixv.col(it) *= vectorv(it);
  }  // end for
  return matrixv.n_cols;
}  // end mult_vec_mat3


//'   It performs an implicit loop over the matrix rows and columns using the
//'   \emph{Armadillo} operators \code{each_row()} and \code{each_col()},
//'   instead of performing explicit \code{for()} loops (both methods are
//'   equally fast).
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat31(arma::vec& vectorv, 
                           arma::mat& matrixv) {
  matrixv.each_row() %= conv_to< rowvec >::from(vectorv);
  return matrixv.n_cols;
}  // end mult_vec_mat31


//'   It performs two explicit \code{for()} loops.
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat32(arma::vec& vectorv, 
                           arma::mat& matrixv) {
  for (uword it = 0; it < matrixv.n_cols; it++) {
    for (uword jt = 0; jt < matrixv.n_rows; jt++) {
      matrixv(jt, it) *= vectorv(it);
    }  // end for
  }  // end for
  return matrixv.n_cols;
}  // end mult_vec_mat32


//'   It performs an implicit loop over the matrix rows and columns using the
//'   \emph{Armadillo} operators \code{each_row()} and \code{each_col()},
//'   instead of performing explicit \code{for()} loops (both methods are
//'   equally fast).
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat4(arma::vec& vectorv, 
                          arma::mat& matrixv) {
  matrixv.each_col() %= vectorv;
  return matrixv.n_cols;
}  // end mult_vec_mat4

//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat5(arma::vec& vectorv, 
                          arma::mat& matrixv) {
  for (uword it = 0; it < matrixv.n_cols; it++) {
    matrixv.col(it) %= vectorv;
  }  // end for
  return matrixv.n_cols;
}  // end mult_vec_mat5


//'   It performs two explicit \code{for()} loops.
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat6(arma::vec& vectorv, 
                          arma::mat& matrixv) {
  for (uword it = 0; it < matrixv.n_rows; it++) {
    for (uword jt = 0; jt < matrixv.n_cols; jt++) {
      matrixv(it, jt) *= vectorv(it);
    }  // end for
  }  // end for
  return matrixv.n_cols;
}  // end mult_vec_mat6


// The function mult_vec_mat2() calculates the product of matrix columns times a
// vector. It copies the vector into matrix columns so it's less efficient.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
int mult_vec_mat_copy(arma::vec& vectorv, arma::mat& matrixv) {
  // uword nrows = matrixv.n_rows;
  uword ncols = matrixv.n_cols;
  arma::mat vec_mat = repmat(vectorv, 1, ncols);
  matrixv = matrixv % vec_mat;
  return matrixv.n_cols;
}  // end mult_vec_mat_copy



// The function mult_vec2_mat_copy() calculates the product of a matrix times
// two vectors.
// It multiplies the rows of the matrix by vectorv1, and the columns by
// vectorv2.
// It produces the same result as the R code: t(t(vectorv2*matrixv)*vectorv1)
// It accepts pointers to the matrix and vectors, assigns the matrix values in
// place without copying the input vector, and returns a double.
// It uses RcppArmadillo.
// Adapted from:
// https://stackoverflow.com/questions/24933290/elementwise-matrix-multiplication-r-versus-rcpp-how-to-speed-this-code-up
//' @export
// [[Rcpp::export]]
int mult_vec2_mat_copy(arma::vec& vectorv2, arma::mat& matrixv, arma::vec& vectorv1) {
  uword nrows = matrixv.n_rows;
  uword ncols = matrixv.n_cols;
  arma::mat vec_mat = repmat(vectorv2, 1, ncols);
  matrixv = matrixv % vec_mat;
  vec_mat = repmat(trans(vectorv1), nrows, 1);
  matrixv = matrixv % vec_mat;
  return matrixv.n_cols;
}  // end mult_vec2_mat_copy


// The function mult_vec2_mat_rcpp() calculates the product of a matrix times
// two vectors.
// It multiplies the rows of the matrix by vectorv1, and the columns by
// vectorv2.
// It produces the same result as the R code: t(t(vectorv2*matrixv)*vectorv1)
// It accepts pointers to the matrix and vectors, assigns the matrix values in
// place without copying the input vector, and returns a double.
// It uses Rcpp.
//' @export
// [[Rcpp::export]]
int mult_vec2_mat_rcpp(NumericVector& vectorv2, NumericMatrix& matrixv, NumericVector& vectorv1) {
  uword nrows = matrixv.nrow();
  uword ncols = matrixv.ncol();
  if (!(nrows == vectorv2.size())) stop("vectorv2 length not equal to number of rows of matrixv");
  if (!(ncols == vectorv1.size())) stop("vectorv1 length not equal to number of columns of matrixv");
  for (uword i = 0; i < nrows; i++) {
    matrixv(i, _) = vectorv1 * matrixv(i, _);
  }  // end for
  for (uword i = 0; i < ncols; i++) {
    matrixv(_, i) = vectorv2 * matrixv(_, i);
  }  // end for
  return ncols;
}  // end mult_vec2_mat_rcpp


// The function mult_vec2_mat_rcpp2() is similar to mat2vec_rcpp_by(), 
// but uses a simple loop without Rcpp, and is much slower.
//' @export
// [[Rcpp::export]]
int mult_vec2_mat_rcpp2(NumericVector& vectorv2, NumericMatrix& matrixv, NumericVector& vectorv1) {
  for (int i = 0; i < matrixv.nrow(); i++) {
    for (int j = 0; j < matrixv.ncol(); j++) {
      matrixv(i, j) = vectorv1(j) * vectorv2(i) * matrixv(i, j);
    }  // end for
  }  // end for
  return matrixv.ncol();
}  // end mult_vec2_mat_rcpp2


// The function get_cor() calculates the correlation of the matrix returns.
//' @export
// [[Rcpp::export]]
arma::mat get_cor(const arma::mat& returns) {
  return arma::cor(returns);
}  // end get_cor


// homework
// The function get_eigenvals() calculates the eigen_values 
// of the matrix covmat.
//' @export
// [[Rcpp::export]]
arma::vec get_eigenvals(const arma::mat& covmat) {
  arma::vec eigen_vals = arma::eig_sym(covmat);
  return eigen_vals;
}  // end get_eigenvals


// The function get_eigen() calculates the eigen decomposition 
// of the matrix returns.
//' @export
// [[Rcpp::export]]
List get_eigen(const arma::mat& returns) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  arma::eig_sym(eigen_val, eigen_vec, cor(returns));
  return List::create(Named("eigval") = eigen_val,
                      Named("eigvec") = eigen_vec);
}  // end get_eigen


// homework
// The function get_pca() calculates the PCA 
// of the matrix returns.
//' @export
// [[Rcpp::export]]
List get_pca(const arma::mat& returns) {
  arma::mat coeff;
  arma::mat sco_re;
  arma::vec la_tent;
  arma::vec t_squared;
  
  arma::princomp(coeff, sco_re, la_tent, t_squared, returns);
  
  return List::create(Named("coefficients") = coeff,
                      Named("score") = sco_re,
                      Named("latent") = la_tent,
                      Named("tsquared") = t_squared);
  
}  // end get_pca


// The function invspd_rcpp() calculates the inverse of symmetric positive
// definite matrix.
// It uses Rcpp and RcppArmadillo.
//' @export
// [[Rcpp::export]]
SEXP invspd_rcpp(SEXP X_) {
  arma::mat X    = Rcpp::as<arma::mat>(X_);
  arma::mat Xinv = arma::inv_sympd(X);
  return(Rcpp::wrap(Xinv));
}  // end invspd_rcpp


// The function invspd_arma() calculates the inverse of symmetric positive 
// definite matrix.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
arma::mat invspd_arma(const arma::mat& matrixv) {
  arma::mat mat_inv = arma::inv_sympd(matrixv);
  return mat_inv;
}  // end invspd_arma


// The function inv_mat() calculates the inverse of symmetric positive
// definite matrix.
// It accepts a pointer to a matrix and operates on the matrix in place.
// It returns the number of columns of the input matrix.
// It uses RcppArmadillo.
//' @export
// [[Rcpp::export]]
double inv_mat(arma::mat& matrixv) {
  matrixv = arma::inv_sympd(matrixv);
  return matrixv.n_cols;
}  // end inv_mat


// The function lm_arma() performs multivariate linear regression, and 
// calculates the alpha and beta coefficients and their t-values and p-values, 
// and the R-squared and F-statistic.
// It uses RcppArmadillo.
// Adapted from:
// http://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
//' @export
// [[Rcpp::export]]
Rcpp::List lm_arma(const arma::colvec& response, const arma::mat& design) {
  // add column for intercept to explanatory matrix
  arma::mat designp = join_rows(ones(design.n_rows), design);
  int nrows = design.n_rows, ncols = designp.n_cols;
  int deg_free = (nrows - ncols);
  
  // fit the model response ~ design, and calculate alpha and beta coefficients
  arma::colvec coeff = arma::solve(designp, response);
  // calculate residuals
  arma::colvec residuals = response - designp*coeff;
  
  // calculate TSS, RSS, and ESS
  double tot_sumsq = (nrows-1)*arma::var(response);
  double res_sumsq = arma::dot(residuals, residuals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // calculate R-squared and F-statistic
  double rsquared = exp_sumsq/tot_sumsq;
  double fstat = (exp_sumsq*deg_free)/(res_sumsq*(ncols-1));
  // arma::rowvec stats=join_horiz(rsquared, fstat);
  Rcpp::NumericVector stats(2);
  stats(0) = rsquared;
  stats(1) = fstat;
  stats.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // calculate standard errors of beta coefficients
  arma::colvec stderr = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(designp)*designp)));
  // calculate t-values and p-values of beta coefficients
  arma::colvec beta_tvals = coeff/stderr;
  arma::colvec beta_pvals = 2*Rcpp::pt(-abs(wrap(beta_tvals)), deg_free);
  NumericMatrix coeffmat = wrap(join_rows(join_rows(join_rows(coeff, stderr), beta_tvals), beta_pvals));
  Rcpp::colnames(coeffmat) = Rcpp::CharacterVector::create("coeff", "stderr", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeffmat,
                            // Named("residuals") = residuals,
                            Named("z_score") = residuals(nrows-1)/arma::stddev(residuals),
                            Named("stats") = stats);
}  // end lm_arma



////////////////////////////
// RcppArmadillo test functions
////////////////////////////


// The function test_rcpp() is for testing some Rcpp code snippets.
//' @export
// [[Rcpp::export]]
LogicalVector test_rcpp(NumericVector& vectorv2, NumericMatrix& matrixv, NumericVector& vectorv1) {
  Rcpp::IntegerVector stats(4);
  stats(0) = matrixv.nrow();
  stats(1) = matrixv.ncol();
  stats(2) = vectorv1.size();
  stats(3) = vectorv2.size();
  stats.attr("names") = Rcpp::CharacterVector::create("nrows", "ncols", "v_siz1", "v_siz2");
  // uword nrows = matrixv.n_row();
  uword ncols = matrixv.ncol();
  // uword v_siz1 = vectorv1.size();
  // uword v_siz2 = vectorv2.size();
  // return stats;
  return !(ncols == vectorv2.size());
}  // end test_rcpp


// The function test_arma() is for testing some RcppArmadillo code snippets.
//' @export
// [[Rcpp::export]]
LogicalVector test_arma(arma::mat& matrixv) {
  // Rcout << "Num rows: " << matrixv.n_rows << std::endl;
  // Rcout << "Num cols: " << matrixv.n_cols << std::endl;
  matrixv.print("This is the input matrix:");
  matrixv(0, 1) = 111;
  matrixv.print("This is the output matrix:");
  return !(matrixv.n_cols == matrixv.n_rows);
}  // end test_arma


//' @export
// [[Rcpp::export]]
arma::uvec test_more_arma(const LogicalVector& vectorv) {
  Function whichv3("whichv3");
  // Rcout << "Num rows: " << matrixv.n_rows << std::endl;
  // Rcout << "Num cols: " << matrixv.n_cols << std::endl;
  // matrixv.print("This is the input matrix:");
  // matrixv(0, 1) = 111;
  // matrixv.print("This is the output matrix:");
  arma::uvec indeks = Rcpp::as<uvec>(whichv3(vectorv));
  return indeks;
}  // end test_more_arma

