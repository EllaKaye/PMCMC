libs <- list("grid", "ggplot2")
lapply(libs,function(x) library(x, character.only = T))

data <- generate_model2(100, sv = sqrt(10), sw = 1)
setwd("./pmcmc/")
Rcpp::sourceCpp("SMC_fast.cpp")
source("PMMH_fast.R")
system.time(
  p <- PMMH_fast(N = 500, data$Y, iters = 10000)
)

df1 <- data.frame(theta_1 = p$theta_chain[,1], theta_2 = p$theta_chain[,2], index = 1:length(p$theta_chain[,1]))



hist_v <- ggplot(data = df1, aes(x = theta_1))+geom_histogram(fill = "white", colour = "black")+
          theme_classic() + ylab(expression(sigma[v])) + 
          theme(axis.title.y=element_text(angle=0), axis.title.x = element_blank(), axis.text = element_blank()) + 
          geom_vline(xintercept = sqrt(10), linetype = "longdash")+
          scale_y_continuous(expand = c(0,0))

hist_w <- ggplot(data = df1, aes(x = theta_2)) + 
  geom_histogram(fill = "white", colour = "black") + theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank())+
  geom_vline(xintercept = sqrt(10), linetype = "longdash")+
  scale_y_continuous(expand = c(0,0))

trace_v <- ggplot(data = df1, aes(y = theta_1, x = index)) + geom_line()+theme_classic()+
  theme(axis.title = element_blank(), axis.text = element_blank())

trace_w <- ggplot(data = df1, aes(y = theta_2, x = index)) + geom_line()+theme_classic()+
  theme(axis.title = element_blank(), axis.text = element_blank())

s_w <- ggplot(data = df1, aes(x = theta_1, y = theta_2)) + geom_point(shape = 4) + theme_classic() + 
ylab(expression(sigma[w])) + theme(axis.title.y=element_text(angle=0), axis.title.x = element_blank())


s_v <- ggplot(data = df1, aes(y = theta_1, x = theta_2)) + geom_point(shape = 4) + theme_classic() +
  theme(axis.title = element_blank(), axis.text.x = element_blank())


grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 4)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(hist_v, vp = vplayout(1, 1))  # key is to define vplayout
print(s_w, vp = vplayout(2, 1))
print(s_v, vp = vplayout(1, 2))
print(hist_w, vp = vplayout(2, 2))
print(trace_v, vp = vplayout(1, 3:4))
print(trace_w, vp = vplayout(2, 3:4))

