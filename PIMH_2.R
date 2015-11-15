# Particle independent MH sampler
PIMH <- function(iter, N, y, dg, rf, rmu=rf) {
  # Initialise the variables
  x <- matrix(NA, nrow = length(y), ncol = iter)
  p <- rep(NA, iter)
  
  # Run SMC algorithm for initial values
  init <- SMC(N, y, dg, rf, rmu)
  # Select one of the particles based on weights
  x[,1] <- init$X_updated[,sample(1:N, size = 1, prob = init$W)]
  p[1] <- init$Marginal.LL
  
  # Loop through the process
  for(i in 2:iter) {
    # Run an SM algorithm for proposal states
    init <- SMC(N, y, dg, rf, rmu)
    x.star <- init$X_updated[,sample(1:N, size = 1, prob = init$W)]
    p.star <- init$Marginal.LL
    
    # Accept with MH probability or reject
    if(runif(1) <= p.star/p[i-1]) {
      x[,i] <- x.star
      p[i] <- p.star
    } else {
      x[,i] <- x[i-1]
      p[i] <- p[i-1]
    }
  }
  return(list("X" = x, "p" = p))
}
th = 1
rmu <- function(n) rnorm(n, mean = 0, sd = 5)
rf <- function(x, t) rnorm(length(x), mean = x, sd = th)
dg <- function(x,y) dnorm(y, mean = x, sd = 1)

PIMH(iter = 1000,N = 5, y = c(2, 2.5, 3, 3), dg, rf, rmu)


# Particle marginal same as above but theta sampled at each step
PMMH <- function(iter, N, y, dg, rf, rmu=rf, rth) {
  x <- matrix(NA, nrow = length(y), ncol = iter)
  p <- rep(NA, iter)
  th <- rth()
  # Run SMC algorithm for initial values
  init <- SMC(N, y, dg, rf, rmu)
  x[,1] <- init$X_updated[,sample(1:N, size = 1, prob = init$W)]
  p[1] <- init$Marginal.LL
  for(i in 2:iter) {
    th <- rth()
    init <- SMC(N, y, dg, rf, rmu)
    x.star <- init$X_updated[,sample(1:N, size = 1, prob = init$W)]
    p.star <- init$Marginal.LL
    if(runif(1) <= p.star/p[i-1]) {
      x[,i] <- x.star
      p[i] <- p.star
    } else {
      x[,i] <- x[i-1]
      p[i] <- p[i-1]
    }
  }
  return(list("X" = x, "p" = p))
}

rth <- function(){rnorm(1,10,1)}

PMMH(iter = 1000,N = 5, y = c(2, 2.5, 3, 3), dg, rf, rmu, rth)

