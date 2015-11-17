library(MCMCpack)

PMMH_fast <- function(N, y, iters) {
  theta1 = 10
  theta2 = 10
  theta_chain = matrix(NA, iters+1, 2)
  theta_chain[1, ] = c(theta1, theta2)
  
  chain <- matrix(NA, iters+1, length(y))
  marginal.ll <- numeric(iters+1)
  accept <- 0
  u <- runif(iters)
  
  # step 1 (i = 0)
  s <- SMC_cpp(N, y, theta1, theta2)
  index <- sample(1:N, size = 1, prob = s$final_weights + 1.0e-16)
  x.star <- s$X_updated[,index]
  marginal.ll[1] <- s$marginal.LL
  chain[1,] <- x.star
  
  # step 2
  for (i in 1:iters) {
    theta1_proposal = max(0, rnorm(1, theta1, 0.15))
    theta2_proposal = max(0, rnorm(1, theta2, 0.08))
    
    s <- SMC_cpp(N, y, theta1_proposal, theta2_proposal)
    index <- sample(1:N, size = 1, prob = s$final_weights)
    x.star <- s$X_updated[,index]
    p.star <- s$marginal.LL
    if(any(c(theta1, theta2, theta1_proposal, theta2_proposal) < 1.0e-3)){
      a = exp(p.star - marginal.ll[i])
    }
    else{
      a = exp(p.star - marginal.ll[i]) * dinvgamma(theta1_proposal**2, 0.01, 0.01) * dinvgamma(theta2_proposal**2, 0.01, 0.01) / dinvgamma(theta1**2, 0.01, 0.01) / dinvgamma(theta2**2, 0.01, 0.01) * abs(theta1_proposal / theta1 * theta2_proposal / theta2)
    }
    
    if (u[i] <= a) {
      theta1 = theta1_proposal
      theta2 = theta2_proposal
      theta_chain[i+1, ] = c(theta1, theta2)
      chain[i+1,] <- x.star
      marginal.ll[i+1] <- p.star
      accept <- accept + 1
    }
    else {
      theta_chain[i+1, ] = c(theta1, theta2)
      chain[i+1,] <- chain[i,]
      marginal.ll[i+1] <- marginal.ll[i]
    }
  }
  
  accept_prob <- accept/iters
  return(list(chain = chain, theta_chain = theta_chain, marginal.ll = marginal.ll, accept = accept_prob))
}

