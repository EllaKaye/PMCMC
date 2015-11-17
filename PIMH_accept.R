library(ggplot2)

# set up N and T
source("data_generation.R")
num_particles <- c(10, 50, 100, 200, 500, 1000, 2000)
np <- length(num_particles)
Time <- c(10, 25, 50, 100)

create_plot = function(sv, sw, Time, num_particles){
  res <- numeric(0)
  for (t in Time) {
    set.seed(1)
    data <- generate_model2(t, sv = sv, sw = sw)
    Y <- data$Y
    
    for (n in num_particles) {
      s <- PIMH_fast(n, Y, sv, sw, iters=1000) 
      res <- c(res, s$accept)
      cat("n =", n, ", t =", t, "\n")
    }
  }
  df <- data.frame(N = num_particles, T = as.factor(rep(Time, each = np)), accept = res)
  return(df)
}

# case 1, sv^2 = 10, sw^2 = 1
df1 = create_plot(sqrt(10), sqrt(10), Time, num_particles)

# case 2, sv^2 = 10, sw^2 = 1
df2 = create_plot(sqrt(10), sqrt(1), Time, num_particles)


# join the two data frames and add variance column for facetting
df <- rbind(cbind(df1, variance = "sv^2 = 10, sw^2 = 10"),
            cbind(df2, variance = "sv^2 = 10, sw^2 = 1"))


# MAKE ANY CHANGES TO PLOT APPEARANCE AS YOU SEE FIT.
p.all <- ggplot(df, aes(N, accept)) + 
  geom_line(aes(group=T, colour=T)) + 
  geom_point(aes(group=T, colour=T)) + 
  facet_grid(.~variance) + 
  ylim(0, 1) + 
  theme_bw()
pdf("pimh_accept.pdf", width = 8, height = 4)
p.all
dev.off()




# last s from nested loops in 5000 iterations of 2000 particles
#plot(colMeans(s$means[1001:5001, ]), type = "l")
#lines(data$X, type = "l", col = "red")





