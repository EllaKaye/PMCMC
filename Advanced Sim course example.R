##
browseURL("https://dl.dropboxusercontent.com/u/1430652/exercise6.hmm.html")

## Particle Filter
# Simulate from model

rho = 0.95
datalength <- 100
Y <- rep(0, datalength)
X <- rep(0, datalength)
X[1] <- rnorm(1)
Y[1] <- rnorm(1, mean = X[1])

for (t in 2:datalength) {
  X[t] <- rnorm(1, mean = rho * X[t-1])
  Y[t] <- rnorm(1, mean = X[t])
}

# Use Kalman filter to calculate the filtering means, E[X_t | Y_{0:t}]
# Use package FKF
library(FKF)
kf.res <- fkf(a0 = 0, P0 = matrix(1), dt = matrix(0), ct = matrix(0), Tt = matrix(rho), Zt = matrix(1), HHt = matrix(1), GGt = matrix(1), yt = matrix(Y, nrow = 1), check.input = TRUE)
kalman.means <- data.frame(mean = as.numeric(kf.res$att), time = 1:datalength)

## Sequential Importance Sampling
# Set up proposal distributions (prior proposal)
rmu <- function(n) rnorm(n, mean = 0, sd = 1)
rf <- function(x) rnorm(length(x), mean = rho * x, sd = 1)
dg <- function(x,y) dnorm(y, mean = x, sd = 1, log = TRUE)

# normalise weights
# remove the maximum from the log weights, to increase numerical stability
normalise_weights <- function(logweights) {
  weights <- exp(logweights - max(logweights))
  weights <- weights / sum(weights)
  return(weights)
}

## SIS - propose new particles using the proposal distribution 
# and computing the new incremental weights.
# here only interested in the filtering means.

SIS <- function(N, datalength) {
  SIS_means <- rep(0, datalength) # filtering means
  
  # step t = 1
  xparticles <- rmu(N)
  logweights <- dg(xparticles, Y[1])
  weights <- normalise_weights(logweights)
  SIS_means[1] <- sum(weights * xparticles)
  
  # step t > 1
  for (t in 2:datalength) {
    xparticles <- rf(xparticles)
    logweights <- logweights + dg(xparticles, Y[t])
    weights <- normalise_weights(logweights)
    SIS_means[t] <- sum(weights * xparticles)
  }
  return(SIS_means)
}

# Run SIS with various number of particles N, 100 times independently
# store results in big data frame.

nrepeats <- 100
Ns <- c(128, 256, 512, 1024)
library(foreach)
sis.df <- foreach (index_repeat = 1:nrepeats, .combine = rbind) %:%
  foreach (N = Ns, .combine = rbind) %do% {
    SIS_means <- SIS(N, datalength)
    data.frame(time = 1:datalength, SIS_means = SIS_means, index_repeat = index_repeat, N=N)
  }

# We want to compare the results with ones obtained with the Kalman filter. 
# To do so, we merge the data frames by time indices. 
# Then we use dplyr to compute the mean squared error for each choice of N and each time. 
# The results are then plotted in a grid.

sis.df <- merge(sis.df, kalman.means, by = "time")
sis.df$index_repeat <- factor(sis.df$index_repeat)
library(dplyr)
library(ggplot2)
MSE.df <- sis.df %>% mutate(SE = (SIS_means - mean)**2) %>%
  group_by(time, N) %>% summarize(MSE = mean(SE))
ggplot(MSE.df, aes(x = time, y = MSE)) + geom_line()+ facet_grid(N ~ .)

# multinomial resampling
multinomial_resampling <- function(normalised_weights) {
  N <- length(normalised_weights)
  return(sample(1:N, N, replace = TRUE, prob = normalised_weights))
  # returns vector of ancestors
}

# Sequential Monte Carlo
SMC <- function(N, datalength) {
  SMC_means <- rep(0, datalength) # filtering means

  # step t = 1
  xparticles <- rmu(N)
  logweights <- dg(xparticles, Y[1])
  weights <- normalise_weights(logweights)
  SMC_means[1] <- sum(weights * xparticles)
  
  # resampling
  xparticles <- xparticles[multinomial_resampling(weights)]
  
  # step t > 1
  for (t in 2:datalength) {
    xparticles <- rf(xparticles)
    logweights <- logweights + dg(xparticles, Y[t])
    weights <- normalise_weights(logweights)
    SMC_means[t] <- sum(weights * xparticles)

    # resampling
    xparticles <- xparticles[multinomial_resampling(weights)]
    
  }
  return(SMC_means)
}

library(dplyr)
smc.df <- foreach (index_repeat = 1:nrepeats, .combine = rbind) %:%
  foreach (N = Ns, .combine = rbind) %do% {
    SMC_means <- SMC(N, datalength)
    data.frame(time = 1:datalength, SMC_means = SMC_means, index_repeat = index_repeat, N = N)
  }
smc.df <- merge(smc.df, kalman.means, by = "time")
smc.df$index_repeat <- factor(smc.df$index_repeat)
MSE.df <- smc.df %>% mutate(SE = (SMC_means - mean)**2) %>%
  group_by(time, N) %>% summarize(MSE = mean(SE))
ggplot(MSE.df, aes(x = time, y = MSE)) + geom_line()+ facet_grid(N ~ .)
