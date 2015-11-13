## Some functions to implement SMC

#' SMC sampler
#' 
#' @param N the number of particles
#' @param y the vector of observations
#' @param dg function which evaluates the density of g
#' @param rf function to sample from the proposal distribution q = f (time t>=2)
#' @param rmu function to sample from the proposal distribution q = mu (time t=1), defaults to rf
#' @return vector of x_{1:T} which approximates the posterior density p_{\theta}(x_{1:T} | y_{1:T})

SMC = function(N, y, dg, rf, rmu=rf){
  T = length(y)
  # matrices for storing vectors of X, W, A
  # (however, X can be saved in two ways: 
  #   1) keeping all the previous generated values, or
  #   2) updating the values for dead particles, i.e. containing repeating values in the previous positions
  X = matrix(NA, T, N)
  W = matrix(NA, T, N)
  A = matrix(NA, T, N)
  X_updated = matrix(NA, T, N)
  
  # time t = 1
  xparticles <- rmu(N)
  weights <- dg(xparticles, y[1])
  
  X[1, ] = xparticles
  W[1, ] = weights
  A[1, ] = 1:N
  X_updated[1, ] = xparticles
  
  # time t > 1
  for(t in 2:T){
    # indexes from multinomial resampling (note: no need to normalise weights here)
    ancestors <- sample(N, replace=TRUE, prob=weights)
    # sample the next particles x_t, given the previous ones x_{t-1} after the resampling step
    xparticles <- rf(xparticles[ancestors])
    # compute weights for the latter
    weights <- dg(xparticles, y[t])
    
    X[t, ] = xparticles
    W[t, ] = weights
    A[t, ] = ancestors
    X_updated = X_updated[, ancestors]
    X_updated[t, ] = xparticles
  }
  # perhaps not necessary to return everything
  return(list("X" = X, "W" = W, "A" = A, "X_updated" = X_updated))
}


# toy example

rmu <- function(n) rnorm(n, mean = 0, sd = 1)
rf <- function(x) rnorm(length(x), mean = x, sd = 1)
dg <- function(x,y) dnorm(y, mean = x, sd = 1)

# note that at the first time point, x is generated from N(0, 1) and 
# only these values are kept which are close to y[1] = 2
SMC(N = 5, y = c(2, 2.5, 3), dg, rf, rmu)
