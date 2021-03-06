library(ggplot2)
library(Rcpp)

source("data_generation.R")
source("PIMH_fast.R")
sourceCpp("SMC_fast.cpp")

# set up N and T
num_particles <- c(10, 50, 100, 200, 500, 1000, 2000)
np <- length(num_particles)
Time <- c(10, 25, 50, 100)

create_plot = function(sv, sw, Time, num_particles, iters){
  res <- numeric(0)
  for (t in Time) {
    set.seed(1)
    data <- generate_model2(t, sv = sv, sw = sw)
    Y <- data$Y
    
    for (n in num_particles) {
      s <- PIMH_fast(n, Y, sv, sw, iters=iters) 
      res <- c(res, s$accept)
      cat("n =", n, ", t =", t, "\n")
    }
  }
  df <- data.frame(N = num_particles, T = as.factor(rep(Time, each = np)), accept = res)
  return(df)
}

# case 1, sv^2 = 10, sw^2 = 10
df1 = create_plot(sqrt(10), sqrt(10), Time, num_particles, iters=10000)

# case 2, sv^2 = 10, sw^2 = 1
df2 = create_plot(sqrt(10), sqrt(1), Time, num_particles, iters=10000)

# save the results
save(df1, df2, file="data/accept_PIMH.RData")

# join the two data frames and add variance column for facetting
df <- rbind(cbind(df1, variance = "(a)"),
            cbind(df2, variance = "(b)"))


# MAKE ANY CHANGES TO PLOT APPEARANCE AS YOU SEE FIT.
p.all <- ggplot(df, aes(N, accept)) + 
  geom_line(aes(group=T, colour=T)) + 
  geom_point(aes(group=T, colour=T, shape=T)) + 
  facet_grid(.~variance) + 
  ylim(0, 1) + ylab("Acceptance rate") +
  xlab("Number of particles") +
  theme_bw() +
  scale_color_brewer(palette="Spectral")
pdf("fig/accept_pimh.pdf", width = 8, height = 3.5)
p.all
dev.off()

