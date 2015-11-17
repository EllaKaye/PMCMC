### overnight plot
# set up N and T
num_particles <- c(10, 50, 100, 200, 500, 1000, 2000)
np <- length(num_particles)
Time <- c(10, 25, 50, 100)

# case 1, sv^2 = 10, sw^2 = 10
sv <- sqrt(10)
sw <- sqrt(10)
accept.10.10 <- numeric(0)
for (t in Time) {
  set.seed(1)
  data <- generate_model2(t, sv = sv, sw = sw)
  Y <- data$Y
  
  for (n in num_particles) {
    s <- PIMH(n, Y, dg, rf, rmu, iters=5000) # CHANGE TO WHATEVER FUNCTION YOU'RE USING - MAYBE UP ITERS
    accept.10.10 <- c(accept.10.10, s$accept)
    cat("n =", n, ", t =", t, "\n")
  }
}
accept.10.10

df.10.10 <- data.frame(N = num_particles, T = as.factor(rep(Time, each = np)), accept = accept.10.10)
p.10.10 <- ggplot(df.10.10, aes(N, accept)) + geom_line(aes(group=T, colour=T)) + geom_point(aes(group=T, colour=T))
p.10.10

# case 1, sv^2 = 10, sw^2 = 1
sw <- 1
accept.10.1 <- numeric(0)
for (t in Time) {
  set.seed(1)
  data <- generate_model2(t, sv = sv, sw = sw)
  Y <- data$Y
  
  for (n in num_particles) {
    s <- PIMH(n, Y, dg, rf, rmu, iters=5000) # CHANGE TO WHATEVER FUNCTION YOU'RE USING - MAYBE UP ITERS
    accept.10.1 <- c(accept.10.1, s$accept)
    cat("n =", n, ", t =", t, "\n")
  }
}
accept.10.1
df.10.1 <- data.frame(N = num_particles, T = as.factor(rep(Time, each = np)), accept = accept.10.1)
p.10.1 <- ggplot(df.10.1, aes(N, accept)) + geom_line(aes(group=T, colour=T)) + geom_point(aes(group=T, colour=T))
p.10.1


# join the two data frames and add variance column for facetting
df <- rbind(cbind(df.10.10, variance = "sv^2 = 10, sw^2 = 10"),
            cbind(df.10.1, variance = "sv^2 = 10, sw^2 = 1"))


# MAKE ANY CHANGES TO PLOT APPEARANCE AS YOU SEE FIT.
p.all <- ggplot(df, aes(N, accept)) + geom_line(aes(group=T, colour=T)) + geom_point(aes(group=T, colour=T)) + facet_grid(.~variance)
pdf("pimh_accept.pdf", width = 8, height = 4)
p.all
dev.off()




# last s from nested loops in 5000 iterations of 2000 particles
#plot(colMeans(s$means[1001:5001, ]), type = "l")
#lines(data$X, type = "l", col = "red")





