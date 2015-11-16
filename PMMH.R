libs <- list("MCMCpack", "MASS", "mvtnorm")
lapply(libs, function(x) library(x, character.only = T))
PMMH <- function(N, y, dg_gen, rf_gen, rmu=rf, iters) {
  chain <- matrix(NA, iters+1, length(y))
  marginal.ll <- numeric(iters+1)
  accept <- 0
  filtering_means <- matrix(NA, iters+1, length(y))
  theta <- matrix(NA, nrow = iters+1, ncol = 2)
  u <- runif(iters)
  # step 1 (i = 0)
  theta[1,] <- c(10,10)
  rf <- rf_gen(theta[1,1])
  dg <- dg_gen(theta[1,2])
  s <- SMC_EK(N, y, dg, rf, rmu)
  index <- sample(1:N, size = 1, prob = s$W)
  x.star <- s$X_updated[,index]
  marginal.ll[1] <- s$Marginal.LL
  chain[1,] <- x.star
  filtering_means[1,] <- s$means
  # step 2
  for (i in 1:iters) {
    theta[i+1,] <- mvrnorm(1, mu = theta[i,], Sigma = diag(c(0.15^2, 0.08^2)))
    rf <- rf_gen(abs(theta[i+1,1]))
    dg <- dg_gen(abs(theta[i+1,2]))
    s <- SMC_EK(N, y, dg, rf, rmu)
    index <- sample(1:N, size = 1, prob = s$W)
    x.star <- s$X_updated[,index]
    p.star <- s$Marginal.LL
    a <- exp(p.star - marginal.ll[i]) * dinvgamma(theta[i+1,1]^2,0.01,0.01) * dinvgamma(theta[i+1,2]^2,0.01,0.01) / (dinvgamma(theta[i,1]^2,0.01,0.01) * dinvgamma(theta[i,1]^2,0.01,0.01))
    
    if (u[i] <= a) {
      chain[i+1,] <- x.star
      marginal.ll[i+1] <- p.star
      filtering_means[i+1,] <- s$means
      accept <- accept + 1
    } else {
      chain[i+1,] <- chain[i,]
      marginal.ll[i+1] <- marginal.ll[i]
      filtering_means[i+1,] <- filtering_means[i,]
      theta[i+1,] <- theta[i,]
    }
    }
  
  accept_prob <- accept/iters
  return(list(chain = chain, marginal.ll = marginal.ll, accept = accept_prob, means = filtering_means, theta = theta))
}

## Example on page 280
sv = sqrt(10)
sw = sqrt(10)
rmu <- function(n) rnorm(n, mean = 0, sd = sqrt(5))
rf <- function(x, t) rnorm(length(x), mean = 0.5*x + 25*x/(1+x**2) + 8*cos(1.2*t), sd = sv)
dg <- function(x,y) dnorm(y, mean = 0.05*x^2, sd = sw, log = TRUE)

rf_gen <- function(sv) {
  return(function(x, t) rnorm(length(x), mean = 0.5*x + 25*x/(1+x**2) + 8*cos(1.2*t), sd = sv))
}
dg_gen <- function(sw) {
  return(function(x,y) dnorm(y, mean = 0.05*x^2, sd = sw, log = TRUE))
}



p <- PMMH(N = 500, y = data$Y, dg_gen, rf_gen, rmu, iters = 1000)

N = 200
y = data$Y
iters = 1000
