library(Rcpp)
sourceCpp("SMC_fast.cpp")

PIMH_fast <- function(N, y, sv, sw, iters) {
  chain <- matrix(NA, iters+1, length(y))
  marginal.ll <- numeric(iters+1)
  accept <- 0
  u <- runif(iters)
  filtering_means <- matrix(NA, iters+1, length(y))

  # step 1 (i = 0)
  s <- SMC_cpp(N, y, sv, sw)
  index <- sample(1:N, size = 1, prob = s$final_weights)
  x.star <- s$X_updated[,index]
  marginal.ll[1] <- s$marginal.LL
  chain[1,] <- x.star
  filtering_means[1,] <- s$filtering_means

  # step 2
  for (i in 1:iters) {
    s <- SMC_cpp(N, y, sv, sw)
    index <- sample(1:N, size = 1, prob = s$final_weights)
    x.star <- s$X_updated[,index]
    p.star <- s$marginal.LL
    #a <- min(1, exp(p.star)/exp(marginal.ll[i]))
    
    if (u[i] <= exp(p.star-marginal.ll[i])) {
      chain[i+1,] <- x.star
      marginal.ll[i+1] <- p.star
      filtering_means[i+1,] <- s$filtering_means
      accept <- accept + 1
    }
    else {
      chain[i+1,] <- chain[i,]
      marginal.ll[i+1] <- marginal.ll[i]
      filtering_means[i+1,] <- filtering_means[i,]
    }
  }
  
  accept_prob <- accept/iters
  return(list(chain = chain, marginal.ll = marginal.ll, accept = accept_prob, filtering_means=filtering_means))
}

