// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

List SMC_cpp(int N, NumericVector y, double sv, double sw) {
  int T = y.size();
  int index;
  double maxlogweight = 0;
  IntegerVector seq_one_to_N = seq_len(N);
  NumericMatrix X(T, N);
  NumericMatrix X_updated(T, N);
  NumericMatrix X_temp(T, N);
  NumericMatrix prevX(T, N);
  NumericMatrix A(T, N);
  NumericMatrix W(T, N);
  NumericVector x(N);
  NumericVector xsq(N);
  NumericVector prevx(N);
  NumericVector logweights(N);
  NumericVector weights(N);
  NumericVector weights_shifted(N);
  NumericVector ancestors(N);
  NumericVector noise(N);
  NumericVector W_rowsums(N);
  NumericVector W_rowmeans(N);
  NumericVector filtering_means_unnormalised(T);
  NumericVector filtering_means(T);
  double loglikelihood = 0.0;
  
  // time t=1
  x = rnorm(N, 0, sv);
  // find x squared (need this for the weights)
  for(int i=0; i<N; i++){
    xsq[i] = pow(x[i], 2);
    logweights[i] = R::dnorm(y[0], 0.05*xsq[i], sw, 1);
  }
  
  maxlogweight = max(logweights); 
  for(int i=0; i<N; i++){
    weights[i] = exp(logweights[i]); 
    weights_shifted[i] = exp(logweights[i] - maxlogweight) + 1.0e-15; 
    W_rowsums[0] += weights[i];
    filtering_means_unnormalised[0] += weights[i] * x[i];
  }
  
  
  X(0, _) = x;
  X_updated(0, _) = x;
  W(0, _) = logweights;
  A(0, _) = seq_len(N);
  
  for(int t=1; t<T; t++){
    ancestors = RcppArmadillo::sample(seq_one_to_N, N, true, weights_shifted);
    
    // update previous x values according to their ancestors
    // essentially prevx <- x
    //std::copy(x.begin(), x.end(), prevx.begin());
    prevx = x;
    prevX = X_updated;
    for(int j=0; j<N; j++){
      index = ancestors[j]-1;
      x[j] = prevx[index];
      // update X_updated as well
      for(int k=0; k<t; k++){
        X_updated(k, j) = prevX(k, index);
      }
    }
    
    // sample new x values
    noise = rnorm(N, 0, sv);
    for(int i=0; i<N; i++){
      x[i] = 0.5*x[i] + 25*x[i]/(1+pow(x[i], 2)) + 8*cos(1.2*(t+1)) + noise[i];
      xsq[i] = pow(x[i], 2);
      logweights[i] = R::dnorm(y[t], 0.05*xsq[i], sw, 1);
    }
    
    
    // shift logweights before taking exponential
    maxlogweight = max(logweights); 
    for(int i=0; i<N; i++){
      weights[i] = exp(logweights[i]); 
      weights_shifted[i] = exp(logweights[i] - maxlogweight) + 1.0e-15; 
      W_rowsums[t] += weights[i];
      filtering_means_unnormalised[t] += weights[i] * x[i];
    }
    
    X(t, _) = x;
    X_updated(t, _) = x;
    W(t, _) = logweights;
    A(t, _) = ancestors;
  }
  
  for(int t=0; t<T; t++){
    W_rowmeans[t] = W_rowsums[t] / N;
    filtering_means[t] = filtering_means_unnormalised[t] / W_rowsums[t]; 
    loglikelihood += log(W_rowmeans[t]);
  }
  NumericVector final_weights = W(T-1, _);
  
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("logW") = W,
                            Rcpp::Named("X") = X, 
                            Rcpp::Named("X_updated") = X_updated, 
                            Rcpp::Named("filtering_means") = filtering_means,
                            Rcpp::Named("final_weights") = exp(final_weights),
                            Rcpp::Named("marginal.LL") = loglikelihood);
}
