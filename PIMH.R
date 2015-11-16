## Some functions to implement SMC

# function to normalise logweights
normalise_weights <- function(logweights) {
  weights <- exp(logweights - max(logweights))
  weights <- weights / sum(weights)
  return(weights)
}

# multinomial resampling
multinomial_resampling <- function(normalised_weights) {
  N <- length(normalised_weights)
  return(sample(1:N, N, replace = TRUE, prob = normalised_weights))
  # returns vector of ancestors
}


## like the SMC sampler, but working in log space and introducing some numerical stability tricks
#' SMC sampler
#' 
#' @param N the number of particles
#' @param y the vector of observations
#' @param dg function which evaluates the log-density of g
#' @param rf function to sample from the proposal distribution q = f (time t>=2)
#' @param rmu function to sample from the proposal distribution q = mu (time t=1), defaults to rf
#' @return vector of x_{1:T} which approximates the posterior density p_{\theta}(x_{1:T} | y_{1:T})

SMC_EK = function(N, y, dg, rf, rmu=rf){
  T = length(y)
  # matrices for storing vectors of X, W, A
  # (however, X can be saved in two ways: 
  #   1) keeping all the previous generated values, or
  #   2) updating the values for dead particles, i.e. containing repeating values in the previous positions
  X = matrix(NA, T, N)
  W = matrix(NA, T, N)
  A = matrix(NA, T, N)
  X_updated = matrix(NA, T, N)
  filtering_means <- rep(NA,T)
  # time t = 1
  xparticles <- rmu(N)
  logweights <- dg(xparticles, y[1])
  weights <- exp(logweights - max(logweights))
  
  X[1, ] = xparticles
  W[1, ] = weights
  A[1, ] = 1:N
  X_updated[1, ] = xparticles
  filtering_means[1] <- sum(W[1, ]/sum(W[1,]) * X_updated[1,])
  # time t > 1
  for(t in 2:T){
    # indexes from multinomial resampling (note: no need to normalise weights here)
    ancestors <- sample(N, replace=TRUE, prob=weights)
    # sample the next particles x_t, given the previous ones x_{t-1} after the resampling step
    xparticles <- rf(xparticles[ancestors], t)
    # compute weights for the latter
    logweights <- dg(xparticles, y[t])
    weights <- exp(logweights - max(logweights))
    
    X[t, ] = xparticles
    W[t, ] = weights
    A[t, ] = ancestors
    X_updated = X_updated[, ancestors]
    X_updated[t, ] = xparticles
    filtering_means[t] <- sum(W[t, ]/sum(W[t,]) * X_updated[t,])
  }
  marginal.ll <- sum(log(rowMeans(W)))
  # perhaps not necessary to return everything
  return(list("X" = X, "W" = W[nrow(W),], "A" = A[-1,], "X_updated" = X_updated, "Marginal.LL" = marginal.ll, means = filtering_means))
}


PIMH <- function(N, y, dg, rf, rmu=rf, iters) {
  chain <- matrix(NA, iters+1, length(y))
  marginal.ll <- numeric(iters+1)
  accept <- 0
  filtering_means <- matrix(NA, iters+1, length(y))
  
  # step 1 (i = 0)
  s <- SMC_EK(N, y, dg, rf, rmu)
  index <- sample(1:N, size = 1, prob = s$W)
  x.star <- s$X_updated[,index]
  marginal.ll[1] <- s$Marginal.LL
  chain[1,] <- x.star
  filtering_means[1,] <- s$means
  
  # step 2
  for (i in 1:iters) {
    s <- SMC_EK(N, y, dg, rf, rmu)
    index <- sample(1:N, size = 1, prob = s$W)
    x.star <- s$X_updated[,index]
    p.star <- s$Marginal.LL
    a <- min(1, exp(p.star)/exp(marginal.ll[i]))
    u <- runif(1)
    
    if (u <= a) {
      chain[i+1,] <- x.star
      marginal.ll[i+1] <- p.star
      filtering_means[i+1,] <- s$means
      accept <- accept + 1
    }
    else {
      chain[i+1,] <- chain[i,]
      marginal.ll[i+1] <- marginal.ll[i]
      filtering_means[i+1,] <- filtering_means[i,]
    }
  }
  
  accept_prob <- accept/iters
  return(list(chain = chain, marginal.ll = marginal.ll, accept = accept_prob, means = filtering_means))
}

## Example on page 280
sv = sqrt(10)
sw = sqrt(10)
rmu <- function(n) rnorm(n, mean = 0, sd = sqrt(5))
rf <- function(x, t) rnorm(length(x), mean = 0.5*x + 25*x/(1+x**2) + 8*cos(1.2*t), sd = sv)
dg <- function(x,y) dnorm(y, mean = 0.05*x^2, sd = sw, log = TRUE)

set.seed(0)
data <- generate_model2(100, sv = sv, sw = sw)
s <- SMC_EK(200, data$Y, dg, rf, rmu)
s$means

system.time(pimh <- PIMH(N = 200, y = data$Y, dg, rf, rmu, iters = 1000))
# 3 sec: 100 data, N = 200, iters = 100
# 30 sec: 100 data, N = 200, iters = 1000
# 300 sec: 100 data, N = 200, iters = 10000
pimh$marginal.ll
pimh$accept

pimh2 <- PIMH(N = 200, y = data$Y, dg, rf, rmu, iters = 1000)
dim(pimh2$means)
