library(Rcpp)
sourceCpp("SMC_fast.cpp")
SMC_cpp(200, y = c(2, 2.5 ,3 ,3), 10, 10)


PIMH_fast <- function(N, y, sv, sw, iters) {
  chain <- matrix(NA, iters+1, length(y))
  marginal.ll <- numeric(iters+1)
  accept <- 0
  u <- runif(iters)

  # step 1 (i = 0)
  s <- SMC_cpp(N, y, sv, sw)
  index <- sample(1:N, size = 1, prob = s$final_weights)
  x.star <- s$X_updated[,index]
  marginal.ll[1] <- s$marginal.LL
  chain[1,] <- x.star

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
      accept <- accept + 1
    }
    else {
      chain[i+1,] <- chain[i,]
      marginal.ll[i+1] <- marginal.ll[i]
    }
  }
  
  accept_prob <- accept/iters
  return(list(chain = chain, marginal.ll = marginal.ll, accept = accept_prob))
}

PIMH_fast(N = 200, y = data$Y, sv = sqrt(10), sw = sqrt(10), iters = 1000)
