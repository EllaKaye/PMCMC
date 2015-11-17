#Set theta 2 to 1

data <- generate_model2(100, sv = 10, sw = 1)
system.time(
  p <- PMMH_fast(N = 5000, data$Y, iters = 100)
)

theta1 <- p$theta_chain[500:1000,1]
theta2 <- p$theta_chain[,2]

df1 <- data.frame(theta_1 = theta1, theta_2 = theta2[500:1000], index = 1:length(theta1))
qplot(theta_1, data = df1)
gp <- ggplot(data = df1, aes(x = theta_1))+geom_histogram(fill = "white", colour = "black")+theme_classic()
print(gp+geom_bar(fill = "white", colour = "black")+geom_density(stat = "step"))





library(grid)



hist_v <- ggplot(data = df1, aes(x = theta_1))+geom_histogram(fill = "white", colour = "black")+
          theme_classic() + ylab(expression(sigma[v])) + 
          theme(axis.title.y=element_text(angle=0), axis.title.x = element_blank(), axis.text.x = element_blank()) + 
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
  theme(axis.title = element_blank(), axis.text = element_blank())


grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 4)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(hist_v, vp = vplayout(1, 1))  # key is to define vplayout
print(s_w, vp = vplayout(2, 1))
print(s_v, vp = vplayout(1, 2))
print(hist_w, vp = vplayout(2, 2))
print(trace_v, vp = vplayout(1, 3:4))
print(trace_w, vp = vplayout(2, 3:4))

