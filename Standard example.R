# consider the model on page 280 of the PMCMC paper

# simulate N observations from the model
model_sim <- function(N, sv, sw) {
  
  X <- numeric(N)
  Y <- numeric(N)
  
  # t = 1
  X[1] <- rnorm(1, mean = 0, sd = 5)
  Y[1] <- X[1]^2/20 + rnorm(1, mean = 0, sd = sw)
  
  # t > 1
  for (t in 2:N) {
    X[t] = X[t-1]/2 + 25 * X[t-1]/(1 + X[t-1]^2) + 8 * cos(1.2 * t) + rnorm(1, mean = 0, sd = sv)
    Y[t] <- X[t]^2/20 + rnorm(1, mean = 0, sd = sw)
  }
  
  return(df)
}