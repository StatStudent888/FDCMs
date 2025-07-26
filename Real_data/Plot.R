# Read data
load(file="out_array_B_new.RData")

y_12 <- numeric(100)
y_113 <- numeric(100)
y_114 <- numeric(100)
y_57 <- numeric(100)
y_610 <- numeric(100)

for (i in 1:100){
  std_dev <- sqrt(diag(out_array_B_new[[i]][[1]]))
  u <- out_array_B_new[[i]][[1]]/ (std_dev %*% t(std_dev))
  y_12[i] <- u[1,2]
  y_113[i] <- u[1,13]
  y_114[i] <- u[1,14]
  y_57[i] <- u[5,7]
  y_610[i] <- u[6,10]
}

library(ggplot2)

data_1 <- data.frame(x = rep(seq(0.01, 1, by = 0.01),5), y = c(y_12,y_113,y_114,y_57,y_610),group = factor(rep(1:5, each = 100)),
                     Entry = factor(rep(c("red", "blue", "green","yellow",'black'), each = 100)))

# Figure2 in section 5.1
plot<-ggplot(data_1, aes(x = x, y = y, color = Entry)) +
  geom_line(linewidth=0.4, alpha=0.8) +
  labs(x = "Fama-French five-factor return vector", y = "Correlation") +
  scale_color_manual(values = unique(data_1$Entry), 
                     labels = c('corr(1,2)', 'corr(1,13)', 'corr(1,14)','corr(5,7)','corr(6,10)') ) +
  scale_x_continuous(
    breaks = c(0.01, 0.25, 0.5, 0.75, 1),  # Set the x-axis breaks
    labels = expression(bolditalic(U)[1], bolditalic(U)[25], bolditalic(U)[50], bolditalic(U)[75], bolditalic(U)[100])  # Subscript formatting
  ) +
  geom_hline(yintercept = 0, color = "#707070", linetype = "dashed")


library(ggcorrplot)
# Figure3 in Appendix F
std_dev <- sqrt(diag(out_array_B_new[[1]][[2]]))
u <- out_array_B_new[[1]][[2]]/ (std_dev %*% t(std_dev))
u<-u[1:50,1:50]
data<-as.data.frame(u)
ggcorrplot(data,tl.cex = 2)

std_dev <- sqrt(diag(out_array_B_new[[50]][[2]]))
u <- out_array_B_new[[50]][[2]]/ (std_dev %*% t(std_dev))
u<-u[1:50,1:50]
data<-as.data.frame(u)
ggcorrplot(data)

std_dev <- sqrt(diag(out_array_B_new[[100]][[2]]))
u <- out_array_B_new[[100]][[2]]/ (std_dev %*% t(std_dev))
u<-u[1:50,1:50]
data<-as.data.frame(u)
ggcorrplot(data,tl.cex = 5)

std_dev <- sqrt(diag(out_array_B_new[[1]][[4]]))
u <- out_array_B_new[[1]][[4]]/ (std_dev %*% t(std_dev))
u<-u[1:50,1:50]
data<-as.data.frame(u)
ggcorrplot(data,tl.cex = 5)



