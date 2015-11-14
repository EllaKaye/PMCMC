# simulate N observations from the model
generate_model1 <- function(N, rho=0.95) {
  
  X <- numeric(N)
  Y <- numeric(N)
  
  # t = 1
  X[1] <- rnorm(1, mean = 0)
  Y[1] <- rnorm(1, mean = X[1])
  
  # t > 1
  for (t in 2:N) {
    X[t] <- rnorm(1, mean = rho * X[t-1])
    Y[t] <- rnorm(1, mean = X[t])
  }
  
  df = data.frame(X=X, Y=Y)
  return(df)
}

# consider the model on page 280 of the PMCMC paper

# simulate N observations from the model
generate_model2 <- function(N, sv, sw) {
  
  X <- numeric(N)
  Y <- numeric(N)
  
  # t = 1
  X[1] <- rnorm(1, mean = 0, sd = sqrt(1))
  Y[1] <- X[1]^2/20 + rnorm(1, mean = 0, sd = sw)
  
  # t > 1
  for (t in 2:N) {
    X[t] = X[t-1]/2 + 25 * X[t-1]/(1 + X[t-1]^2) + 8 * cos(1.2 * t) + rnorm(1, mean = 0, sd = sv)
    Y[t] <- X[t]^2/20 + rnorm(1, mean = 0, sd = sw)
  }
  
  df = data.frame(X=X, Y=Y)
  return(df)
}
