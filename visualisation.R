library(reshape2)
library(ggplot2)
library(dplyr)

source("data_generation.R")
source("SMC_functions.R")

create_plot = function(df.m, df.discarded.m, maxt, ylim){
  if(!is.null(df.discarded.m)){
    p0 = ggplot(df.m, aes(t, value, group=variable)) + 
      geom_point(aes(size=n), data=df.discarded.m, col="grey") +
      geom_path(data=df.discarded.m, col="grey", linetype="dashed") + 
      geom_point(aes(size=n)) + 
      geom_path()
  }
  else{
    p0 = ggplot(df.m, aes(t, value, group=variable)) + 
      geom_point(aes(size=n)) + 
      geom_path()
  }
  p = p0 + theme_classic() + theme(legend.position="none") +
    scale_x_continuous(lim = c(1, maxt), breaks=1:maxt) + 
    scale_y_continuous(lim =  ylim) +
    ylab(expression(X[t]))
  return(p)
}

visualise_SMC = function(obj, maxt=5){
  X = obj$X
  A = obj$A
  df.discarded.m = data.frame(c())
  ylim = range(X)
  
  # when t = 1
  x = X[1, ]
  a = A[1, ]
  X_updated = rbind(x)
  df = as.data.frame(cbind(X_updated, t=1))
  df.m = melt(df, id.vars="t") %>%
    group_by(t, value) %>%
    mutate(n = n())
  
  p = create_plot(df.m, NULL, maxt, ylim)
  print(p)
  
  for(t in 2:maxt){
    x = X[t, ]
    a = A[t, ]
    X_updated = X_updated[, a, drop=FALSE]
    
    df.discarded = df[, c(setdiff(1:(ncol(df)-1), a), ncol(df))]
    temp = melt(df.discarded, id.vars="t") %>%
      group_by(t, value) %>%
      mutate(n = n(), 
             variable = as.numeric(variable)+ifelse(is.null(df.discarded.m$variable), 0, max(df.discarded.m$variable)))
    df.discarded.m = rbind(df.discarded.m, temp)
    
    df = data.frame(cbind(X_updated, t=1:(t-1)))
    df.m = melt(df, id.vars="t") %>%
      group_by(t, value) %>%
      mutate(n = n())
    
    p = create_plot(df.m, df.discarded.m, maxt, ylim)
    print(p)
    
    X_updated = rbind(X_updated, x)
    df = data.frame(cbind(X_updated, t=1:t))
    df.m = melt(df, id.vars="t") %>%
      group_by(t, value) %>%
      mutate(n = n())
    
    p = create_plot(df.m, df.discarded.m, maxt, ylim)
    print(p)
  }
}

# Visualise SMC for the first example

rmu <- function(n) rnorm(n, mean = 0, sd = 1)
rf <- function(x, t) rnorm(length(x), mean = 0.95*x, sd = 1)
dg <- function(x,y) dnorm(y, mean = x, sd = 1)

set.seed(0)
data <- generate_model1(100, rho = 0.95)
obj <- SMC(N = 5, data$Y, dg, rf, rmu)
visualise_SMC(obj, maxt=10)

# Visualise SMC for the second example

sv = 1
sw = 1
rmu <- function(n) rnorm(n, mean = 0, sd = sv)
rf <- function(x, t) rnorm(length(x), mean = 0.5*x + 25*x/(1+x**2) + 8*cos(1.2*t), sd = sv)
dg <- function(x,y) dnorm(y, mean = 0.05*x^2, sd = sw)

set.seed(0)
data <- generate_model2(100, sv = sv, sw = sw)
obj <- SMC(N = 5, data$Y, dg, rf, rmu)
visualise_SMC(obj, maxt=10)

