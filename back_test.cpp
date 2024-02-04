////////////////////////////
// Functions to backtest trading strategies.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/back_test.cpp")

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


// The function bollinger_strat() simulates a Bollinger strategy 
// and returns its PnLs, positions, and z-scores.
// The strategy buys 1 share if the z-score is greater than the threshold, 
// and sells if the z-score is less than minus the threshold.
// The z-score is equal to the current price minus the last fill price, divided by the volatility.  
// The variance is the sum of the squared differences of the current price minus the EMA price.  
// It uses Rcpp.
//' @export
 // [[Rcpp::export]]
 arma::mat bollinger_strat(const arma::mat& pricev, // Time series of prices
                           double lambdaf, // Decay factor which multiplies the past values
                           double threshv = 1.0, // Threshold level
                           double varin = 100.0) { // Initial variance
   
   arma::uword nrows = pricev.n_rows;
   double lambda1 = 1-lambdaf;
   double vars = varin;
   double pricefill = pricev(0);
   double pricema = pricev(0);
   double retv;
   arma::mat zscores = arma::zeros(nrows, 1);
   arma::mat posv = arma::zeros(nrows, 1);
   arma::mat pnlv = arma::zeros(nrows, 1);
   // arma::mat pnlv = arma::zeros(nrows, 1);
   
   // Calculate the positions in a loop
   for (arma::uword it = 1; it < (nrows-1); it++) {
     retv = pricev(it) - pricev(it-1);
     // Update the pnls
     pnlv(it) = posv(it)*retv;
     // Update the z-score
     // zscores(it) = (pricev(it) - pricefill) / sqrt(vars);
     zscores(it) = (pricev(it) - pricema) / sqrt(vars);
     // Update the variance
     // vars = lambdaf*vars + lambda1*pow(pricev(it) - pricema, 2);
     vars = lambdaf*vars + lambda1*pow(retv, 2);
     // Update the EMA price
     pricema = lambdaf*pricema + lambda1*pricev(it);
     if ((zscores(it) > threshv) && (true)) {
       // Sell 1 share
       posv(it+1) = posv(it) - 1;
       // Update the fill price
       pricefill = pricev(it);
     } else if ((zscores(it) < (-threshv)) && (true)) {
       // Buy 1 share
       posv(it+1) = posv(it) + 1;
       // Update the fill price
       pricefill = pricev(it);
     } else {
       // Do nothing
       posv(it+1) = posv(it);
     }  // end if
   }  // end for
   
   return arma::join_rows(pnlv, posv, zscores);
   
 }  // end bollinger_strat



// The function contrastrat() simulates a contrarian (mean reverting) 
// strategy using z-scores and thresholds. 
// It returns the strategy PnLs, positions, and z-scores.
// The strategy sells short 1 share if the z-score is greater than the 
// threshold, and buys 1 share if the z-score is less than minus the 
// threshold.
// The z-score is equal to the current price minus the reference price, 
// divided by the trailing volatility of prices.
// The reference price can be chosen as the EMA price or the trade fill
// price.
// The strategy also filters (scrubs) bad prices. 
// It ignores prices which have an absolute z-score greater than the 
// threshold level threshbad.
// It uses Rcpp.
//' @export
 // [[Rcpp::export]]
 arma::mat contrastrat(const arma::mat& pricev, // Time series of prices
                       double lambdaf, // Decay factor which multiplies the past values
                       double threshv = 1.0, // Threshold level
                       double threshd = 1.0, // Threshold double down
                       double threshbad = 10.0, // Threshold bad
                       double varin = 1.0) { // Initial variance
   
   bool isvalid;
   arma::uword nrows = pricev.n_rows;
   double lambda1 = 1-lambdaf;
   double vars = varin;
   double pricefill = pricev(0);
   double pricema = pricev(0);
   // double retv;
   arma::mat zscores = arma::zeros(nrows, 1);
   arma::mat posv = arma::zeros(nrows, 1);
   arma::mat pnlv = arma::zeros(nrows, 1);
   arma::mat retv = arma::zeros(nrows, 1);
   // arma::mat pnlv = arma::zeros(nrows, 1);
   
   // Calculate the positions in a loop
   for (arma::uword it = 1; it < (nrows-1); it++) {
     retv(it) = pricev(it) - pricev(it-1);
     // Update the pnls
     pnlv(it) = posv(it)*retv(it);
     // Calculate the z-score using the EMA price
     zscores(it) = (pricev(it) - pricema) / sqrt(vars);
     // Check if the price is valid
     // isvalid = (fabs(zscores(it)) < threshbad);
     if (fabs(zscores(it)) < threshbad) {
       // Calculate the z-score using the fill price
       // zscores(it) = (pricev(it) - pricefill) / sqrt(vars);
       // Calculate the z-score using the EMA price
       zscores(it) = (pricev(it) - pricema) / sqrt(vars);
       // Update the variance
       vars = lambdaf*vars + lambda1*pow(pricev(it) - pricema, 2);
       // vars = lambdaf*vars + lambda1*pow(retv, 2);
       // Update the EMA price
       pricema = lambdaf*pricema + lambda1*pricev(it);
       // if (isvalid && (((zscores(it) > threshv) && (posv(it) >= 0)) || (zscores(it) > threshd))) {
       if (zscores(it) > threshv) {
         // Sell 1 share
         // posv(it+1) = posv(it) - 1;
         posv(it+1) = -1;
         // Update the fill price
         pricefill = pricev(it);
         // } else if (isvalid && (((zscores(it) < (-threshv)) && (posv(it) <= 0)) || (zscores(it) < (-threshd)))) {
       } else if (zscores(it) < (-threshv)) {
         // Buy 1 share
         // posv(it+1) = posv(it) + 1;
         posv(it+1) = 1;
         // Update the fill price
         pricefill = pricev(it);
       } else {
         // Do nothing - carry over the current position
         posv(it+1) = posv(it);
       }  // end if
     } else {
       // Price is invalid - carry over the position
       retv(it) = retv(it-1);
       posv(it+1) = posv(it);
     }  // end if
   }  // end for
   
   return arma::join_rows(pnlv, posv, zscores, retv);
   
 }  // end contrastrat



// The function bollinger_brackets() simulates a Bollinger price brackets strategy 
// and returns its positions.
// The Bollinger brackets strategy buys shares at the buy price and sells at the sell price.
// The buy price is equal to the previous trade fill price minus the standard deviation of returns.
// The sell price is equal to the previous fill price plus the standard deviation.
// It uses Rcpp.
//' @export
 // [[Rcpp::export]]
 arma::mat bollinger_brackets(const arma::mat& retv, // Time series of returns
                              double lambdaf, // Decay factor which multiplies the past values 
                              double varf, // Variance factor
                              double varin) { // Initial variance
   
   arma::uword nrows = retv.n_rows;
   arma::mat pricec = arma::zeros(nrows, 1);
   arma::mat vars = arma::zeros(nrows, 1);
   arma::mat posv = arma::zeros(nrows, 1);
   // arma::mat pnlv = arma::zeros(nrows, 1);
   double lambda1 = 1-lambdaf;
   double volv = sqrt(varin);
   double pricesell = volv;
   double pricebuy = -volv;
   
   pricec(0) = retv(0);
   vars(0) = varin;
   
   // Calculate the positions in a loop
   for (arma::uword it = 1; it < (nrows-1); it++) {
     // Update the pnls
     // pnlv(it) = posv(it)*retv(it);
     // Update the prices
     pricec(it) = pricec(it-1) + retv(it);
     // Update the variance
     vars(it) = lambdaf*vars(it-1) + varf*lambda1*pow(retv(it) - retv(it-1), 2);
     volv = sqrt(vars(it));
     if ((pricec(it) > pricesell) && (pricec(it-1) > pricesell)) {
       // Sell 1 share
       posv(it+1) = posv(it) - 1;
       // Update the buy and sell prices
       pricesell = pricec(it) + volv;
       pricebuy = pricec(it) - volv;
     } else if ((pricec(it) < pricebuy) && (pricec(it-1) < pricebuy)) {
       // Buy 1 share
       posv(it+1) = posv(it) + 1;
       // Update the buy and sell prices
       pricesell = pricec(it) + volv;
       pricebuy = pricec(it) - volv;
     } else {
       // Do nothing
       posv(it+1) = posv(it);
     }  // end if
   }  // end for
   
   return posv;
   
 }  // end bollinger_brackets


// The function revert_to_open() calculates the positions of a mean-reversion strategy.
// It bets on prices reverting to the open price.
// It uses Rcpp.
//' @export
 // [[Rcpp::export]]
 arma::mat revert_to_open(const arma::mat& pricev, // Time series of prices
                          double lambdaf, // Decay factor which multiplies the past values 
                          double varf, // Variance factor
                          double varin) { // Initial variance
   
   arma::uword nrows = pricev.n_rows;
   arma::mat retv = arma::zeros(nrows, 1);
   // arma::mat cumax = arma::zeros(nrows, 1);
   // arma::mat cumin = arma::zeros(nrows, 1);
   // arma::mat vars = arma::zeros(nrows, 1);
   arma::mat posv = arma::zeros(nrows, 1);
   // arma::mat pnlv = arma::zeros(nrows, 1);
   // double lambda1 = 1-lambdaf;
   // double volv = sqrt(varin);
   double pricinit = pricev(0);
   // double pricesell = volv;
   // double pricebuy = -volv;
   
   // cumax(0) = pricinit;
   // cumin(0) = pricinit;
   // vars(0) = varin;
   
   // Calculate the positions in a loop
   for (arma::uword it = 1; it < (nrows-1); it++) {
     // Update the pnls
     // pnlv(it) = posv(it)*retv(it);
     // Update the returns
     retv(it) = pricev(it) - pricev(it-1);
     // Update the variance
     // vars(it) = lambdaf*vars(it-1) + varf*lambda1*pow(retv(it) - retv(it-1), 2);
     // volv = sqrt(vars(it));
     if (pricev(it) > pricinit) {
       // Unwind positive position
       if (posv(it) > 0) posv(it) = 0;
       if (retv(it) > 0) {
         // Add to negative position
         posv(it) = posv(it-1) - retv(it);
       } else {
         // Do nothing
         posv(it) = posv(it-1);
       }  // end if
     } else if (pricev(it) <= pricinit) {
       // Unwind negative position
       if (posv(it) < 0) posv(it) = 0;
       if (retv(it) < 0) {
         // Add to positive position
         posv(it) = posv(it-1) - retv(it);
       } else {
         // Do nothing
         posv(it) = posv(it-1);
       }  // end if
     }  // end if
   }  // end for
   
   return posv;
   
 }  // end revert_to_open


