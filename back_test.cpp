////////////////////////////
// Functions to backtest trading strategies.
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/Rcpp/back_test.cpp")

#include "RcppArmadillo.h"
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


// The function trend_flip() simulates a trending strategy 
// and returns its PnLs, positions, and variance.
// The strategy follows the price trend, and reverses its 
// position when the trend reverses.
// It waits for a clear trend reversal based on the price 
// volatility.
// The strategy buys shares if the price gain above the 
// recent minimum price is greater than the volatility. 
// It sells shares if the price drawdown below the recent 
// maximum price is greater than the volatility.
// It also scrubs isolated price spikes based on the return 
// (price difference).
// It ignores prices which have an absolute return greater 
// than the threshold level threshbad.
// It uses the Boolean variable isvalid to prevent scrubbing 
// more than one isolated price spike at a time.
// It uses Rcpp.
//' @export
// [[Rcpp::export]]
arma::mat trend_flip(const arma::mat& pricev, // Time series of prices
                     double lambdaf, // Decay factor which multiplies the past values
                     int poslimit = 2, // Position limit
                     double volf = 0.2, // Volatility floor in dollars
                     double threshbad = 0.1) { // Threshold bad

  bool isvalid = true; // Check if the price is valid
  arma::uword nrows = pricev.n_rows;
  double lambda1 = 1-lambdaf;
  // double pnlreal = 0;
  // double pnlunreal = 0;
  double pricec = pricev(0);
  double pricemax = pricec;
  double pricemin = pricec;
  double pricedrawd = 0;
  double pricegain = 0;
  double zerov = 0; // Zero value to make max and min work
  // double pricefill = pricec;
  double pricema = pricec; // Moving average price
  arma::mat vars = arma::zeros(nrows, 1);
  vars(0) = pow(pricev(1) - pricec, 2);
  double volv;
  double retv;
  // arma::mat zscores = arma::zeros(nrows, 1);
  arma::mat posv = arma::zeros(nrows, 1);
  arma::mat pnlv = arma::zeros(nrows, 1);
  // arma::mat pnlv = arma::zeros(nrows, 1);
  
  // Calculate the positions in a loop
  for (arma::uword it = 1; it < (nrows-1); it++) {
    pricec = pricev(it);
    retv = pricec - pricev(it-1);
    // Check if the price is valid
    if ((fabs(retv) > threshbad) && isvalid) {
      isvalid = false;
    } else {
      // The price is valid
      isvalid = true;
      // Update the pnls
      pnlv(it) = posv(it)*retv;
      // Update the z-score
      volv = std::max(sqrt(vars(it)), volf);
      // zscores(it) = (pricec - pricefill) / volv;
      // zscores(it) = (pricec - pricema) / volv;
      // Update the variance
      vars(it+1) = lambdaf*vars(it) + lambda1*pow(pricec - pricema, 2);
      // vars = lambdaf*vars + lambda1*pow(retv, 2);
      // Update the EMA price
      pricema = lambdaf*pricema + lambda1*pricec;
      // Update the prices
      pricemax = std::max(pricemax, pricec);
      pricemin = std::min(pricemin, pricec);
      pricedrawd = std::max(pricemax - pricec, zerov);
      pricegain = std::max(pricec - pricemin, zerov);
      if ((pricedrawd > volv) && (posv(it) >= 0)) {
        // Buy 1 share
        posv(it+1) = -poslimit;
        pricemin = pricec;
        // Update the fill price
        // pricefill = pricec;
      } else if ((pricegain > volv) && (posv(it) <= 0)) {
        // Sell 1 share
        posv(it+1) = poslimit;
        pricemax = pricec;
        // Update the fill price
        // pricefill = pricec;
      } else {
        // Do nothing
        posv(it+1) = posv(it);
      }  // end if
    }  // end if the price is valid
  }  // end for
  
  return arma::join_rows(pnlv, posv, vars);
  
}  // end trend_flip



// The function bollinger_double() simulates a Bollinger double down 
// strategy.
// The strategy sells short 1 share if the z-score is greater than 
// the threshold, and buys 1 share if the z-score is less than minus 
// the threshold.
// The strategy doubles down - it keeps buying stocks as prices continue 
// dropping, and it keeps selling as prices continue rising. 
// But it applies a brake on consecutive trades in the same direction 
// (consecutive buys or sells).  
// It submits consecutive buys or sells only under two conditions.  
// First, if the number of consecutive trades does not exceed the limit tradelimit.  
// Second, if the absolute z-score exceeds the double down threshold threshd, 
// even if it exceeds the max number of consecutive trades.
// The z-score is equal to the current price minus the reference price, 
// divided by the trailing volatility of prices.
// The reference price can be chosen as the EMA price or the trade 
// fill price.
// It returns the strategy PnLs, positions, and z-scores.
// It uses Rcpp.
//' @export
// [[Rcpp::export]]
arma::mat bollinger_double(const arma::mat& pricev, // Time series of prices
                           double lambdaf, // Decay factor which multiplies the past values
                           double threshv = 1.0, // Threshold level
                           double threshd = 2.0, // Threshold double down
                           double maxd = 1.0, // Max number of double downs
                           double volf = 1.0) { // Volatility floor in dollars
  
  // bool isvalid;
  arma::uword nrows = pricev.n_rows;
  double lambda1 = 1-lambdaf;
  double vars = pow(volf, 2); // Variance
  double volv; // Volatility
  double pricefill = pricev(0);
  double pricema = pricev(0); // Moving average price
  double nbuys = 0.0;
  double nsells = 0.0;
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
    // zscores(it) = (pricev(it) - pricema) / sqrt(vars);
    // Check if the price is valid
    // isvalid = (fabs(zscores(it)) < threshbad);
    // if (fabs(zscores(it)) < threshbad) {
    // Calculate the z-score using the fill price
    volv = std::max(sqrt(vars), volf);
    zscores(it) = (pricev(it) - pricefill) / sqrt(vars);
    // Calculate the z-score using the EMA price
    // zscores(it) = (pricev(it) - pricema) / sqrt(vars);
    // Update the variance
    vars = lambdaf*vars + lambda1*pow(pricev(it) - pricema, 2);
    // vars = lambdaf*vars + lambda1*pow(retv, 2);
    // Update the EMA price
    pricema = lambdaf*pricema + lambda1*pricev(it);
    // if ((((zscores(it) > threshv) && (posv(it) >= 0)) || (zscores(it) > threshd))) {
    if ((((zscores(it) > threshv) && (nsells < maxd)) || (zscores(it) > threshd))) {
      // if (zscores(it) > threshv) {
      // Sell 1 share
      posv(it+1) = posv(it) - 1;
      nsells += 1;
      nbuys = 0;
      // posv(it+1) = -1;
      // Update the fill price
      pricefill = pricev(it);
    } else if ((((zscores(it) < (-threshv)) && (nbuys < maxd)) || (zscores(it) < (-threshd)))) {
      // } else if (zscores(it) < (-threshv)) {
      // Buy 1 share
      posv(it+1) = posv(it) + 1;
      nbuys += 1;
      nsells = 0;
      // posv(it+1) = 1;
      // Update the fill price
      pricefill = pricev(it);
    } else {
      // Do nothing - carry over the current position
      posv(it+1) = posv(it);
    }  // end if
    // } else {
    //   // Price is invalid - carry over the position
    //   retv(it) = retv(it-1);
    //   posv(it+1) = posv(it);
    // }  // end if
  }  // end for
  
  // return arma::join_rows(pnlv, posv, zscores, retv);
  return arma::join_rows(pnlv, posv);
  
}  // end bollinger_double



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
                              double volf) { // Volatility floor in dollars
     
   arma::uword nrows = retv.n_rows;
   arma::mat pricec = arma::zeros(nrows, 1);
   arma::mat vars = arma::zeros(nrows, 1);
   arma::mat posv = arma::zeros(nrows, 1);
   // arma::mat pnlv = arma::zeros(nrows, 1);
   double lambda1 = 1-lambdaf;
   double volv = volf;
   double pricesell = volv;
   double pricebuy = -volv;
   
   pricec(0) = retv(0);
   vars(0) = pow(volf, 2);
   
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
                         double volv = 0.1) { // Volatility

  arma::uword nrows = pricev.n_rows;
  double zscore = 0.0;
  double zscorep = 0.0; // Previous z-score
  double zscorem = 0.0; // Extreme z-score
  double pricinit = pricev(0);
  arma::mat pnlv = arma::zeros(nrows, 1);
  arma::mat posv = arma::zeros(nrows, 1);

  // Calculate the positions in a loop
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the z-scores
    zscorep = zscore;
    zscore = (pricev(it) - pricinit)/volv;
    // zscore = sign(zscore)*pow(zscore, 0.1);
    // zscore = pow(zscore, 3);
    // Check if z-score has changed sign
    if ((zscore*zscorep) > 0) {
      // The z-score has not changed sign
      // Update the extreme z-score
      if (zscore > 0) {
        zscorem = std::max(zscore, zscorem);
      } else {
        zscorem = std::min(zscore, zscorem);
      }  // end if
    } else {
      // The z-score has changed sign
      zscorem = 0.0;
    }  // end if
    posv(it) = -zscorem;
    // Update the pnls
    pnlv(it) = posv(it-1)*(pricev(it) - pricev(it-1));
  }  // end for
  
  return arma::join_rows(pnlv, posv);
  
}  // end revert_to_open


