library(ggplot2)
library(gridExtra)
library(Rcpp)

source("data_generation.R")
source("PIMH_fast.R")
sourceCpp("SMC_fast.cpp")

# define functions needed for data generation
rmu <- function(n) rnorm(n, mean = 0, sd = sv)
rf <- function(x, t) rnorm(length(x), 0.5*x + 25*x/(1+x**2) + 8*cos(1.2*t), sd = sv)
dg <- function(x,y) dnorm(y, mean = 0.05*x^2, sd = sw)




create_df = function(sv, sw){
  set.seed(0)
  # generate data
  data <- generate_model2(100, sv = sv, sw = sw)
  # run PIMH sampler
  res = PIMH_fast(N = 200, y = data$Y, sv = sv, sw = sw, iters = 10000)
  # remove some burn-in
  mat = res$filtering_means[-c(1:1000), ]

  df = data.frame(t = 1:length(data$X), 
                  true_x = data$X, 
                  means = colMeans(mat))
  return(df)
}

df1 = create_df(sv = sqrt(10), sw = sqrt(10))
df2 = create_df(sv = sqrt(0.01), sw = sqrt(10))
df = rbind(cbind(df1, type="(a)"), 
           cbind(df2, type="(b)"))
save(df, file="data/meanX.RData")

p = ggplot(df) + 
  geom_ribbon(aes(x=t, ymin=pmin(true_x, means), ymax=pmax(true_x, means)), alpha=0.25, fill="red") + 
  geom_path(aes(t, true_x, col="True X")) + 
  geom_path(aes(x=t, y=means, col="PIMH mean")) + 
  facet_wrap(~ type, ncol=1, scales="free") + 
  scale_color_manual("", values = c("black", "red")) +
  ylab("value of X") +
  theme_bw()

pdf("fig/meanX.pdf", width = 8, height = 6)
p
dev.off()
