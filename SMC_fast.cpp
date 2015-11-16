// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

List SMC_cpp(int N, NumericVector y, double sv, double sw) {
  int T = y.size();
  int index;
  IntegerVector seq_one_to_N = seq_len(N);
  NumericMatrix X(T, N);
  NumericMatrix X_updated(T, N);
  NumericMatrix X_temp(T, N);
  NumericMatrix A(T, N);
  NumericMatrix W(T, N);
  NumericVector x(N);
  NumericVector xsq(N);
  NumericVector newx(N);
  NumericVector weights(N);
  NumericVector ancestors(N);
  NumericVector noise(N);
  
  // time t=1
  x = rnorm(N, 0, sv);
  // find x squared (need this for the weights)
  for(int i=0; i<N; i++){
    xsq[i] = pow(x[i], 2);
    weights[i] = R::dnorm(y[0], 0.05*xsq[i], sw, 0);
  }
  X(0, _) = x;
  X_updated(0, _) = x;
  W(0, _) = weights;
  A(0, _) = seq_len(N);
  
  // time t>1
  for(int t=1; t<T; t++){
    ancestors = RcppArmadillo::sample(seq_one_to_N, N, true, weights);
    
    // update previous x values according to their ancestors
    for(int j=0; j<N; j++){
      index = ancestors[j]-1;
      newx[j] = x[index];
      // update X_updated as well
      for(int k=0; k<t; k++){
        X_temp(k, j) = X_updated(k, index);
      }
    }
    x = newx;
    X_updated = X_temp;
    
    // sample new x values
    noise = rnorm(N, 0, sv);
    for(int i=0; i<N; i++){
      x[i] = 0.5*x[i] + 25*x[i]/(1+pow(x[i], 2)) + 8*cos(1.2*t) + noise[i];
      xsq[i] = pow(x[i], 2);
      weights[i] = R::dnorm(y[t], 0.05*xsq[i], sw, 0) + 1.0e-15;
    }
    
    X(t, _) = x;
    X_updated(t, _) = x;
    W(t, _) = weights;
    A(t, _) = ancestors;
  }
  
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("W") = W,
                            Rcpp::Named("X") = X, 
                            Rcpp::Named("X_updated") = X_updated);
}
