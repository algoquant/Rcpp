
////////////////////////////////////////////////////////////
// The function trend_follow() simulates a trending following 
// strategy and returns its PnLs, positions, and variance.
// The strategy follows the price trend, and reverses its 
// position when the trend reverses.
// It waits for a clear trend reversal based on the price 
// volatility.
// The strategy buys shares if the price gain above the 
// reference price is greater than the price volatility. 
// It sells shares if the price drawdown below the reference
// price is greater than the volatility.
// It's written in Rcpp.
//' @export
// [[Rcpp::export]]
arma::mat trend_follow(const arma::mat& pricev, // Time series of prices
                       Rcpp::Nullable<double> lambdaf = R_NilValue) { // Decay factor which multiplies the past values
  
  arma::uword nrows = pricev.n_rows;
  double pricinit = pricev(0); // Initial price
  double priceref = pricev(0); // Reference price
  double pricema = pricev(0); // EMA price
  bool use_ema = lambdaf.isNotNull();
  double lambdafv = use_ema ? Rcpp::as<double>(lambdaf) : 0.0;
  double lambda1 = use_ema ? (1.0 - lambdafv) : 0.0;
  double pricec = pricinit; // Current price
  // arma::mat varv = arma::zeros(nrows, 1); // Price variance
  // varv(0) = pow(pricev(1) - pricec, 2);
  // double volv; // Price volatility
  double retp; // Stock return
  // arma::mat zscores = arma::zeros(nrows, 1);
  arma::mat posv = arma::zeros(nrows, 1);
  arma::mat pnlv = arma::zeros(nrows, 1);
  // arma::mat pnlv = arma::zeros(nrows, 1);
  
  // Calculate the PnLs in a loop
  for (arma::uword it = 1; it < (nrows-1); it++) {
    pricec = pricev(it);
    retp = pricec - pricev(it-1);
    pnlv(it) = posv(it-1) * retp;
    
    if (use_ema) {
      pricema = lambda1 * pricema + lambdafv * pricec;
      priceref = pricema;
    } // end if
    
    // varv(it) = lambda1 * varv(it-1) + lambdaf * pow(pricec - pricema, 2);
    // volv = sqrt(varv(it));
    // zscores(it) = (pricec - pricema) / volv;
    posv(it) = arma::sign(pricec - priceref);
  }  // end for
  
  // return arma::join_rows(pnlv, posv, varv);
  return arma::join_rows(pnlv, posv);
  
}  // end trend_follow



