library(readxl)
library(ggplot2)
library(tidyr)
# Read data
data <- read_excel("real_data_2.xlsx",col_names = FALSE)

matrix_data <- as.matrix(data)
X = matrix_data[, (ncol(matrix_data) - 4):ncol(matrix_data)]
train_X = X[1:100,]

mydata <- as.data.frame(train_X)
colnames(mydata) <- c("Mkt-RF", "SMB", "HML", "RMW", "CMA")
mydata$day <- 1:100

# Reshape data to long format for ggplot
data_long <- pivot_longer(mydata, cols = -day, names_to = "Factor", values_to = "Value")

# Figure1 in section 5.1
plot <- ggplot(data_long, aes(x = day, y = Value, color = Factor)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +  # Increase line size and add transparency for better visibility
  labs(x = "Time (trading days)", y = "Factor return") +
  scale_x_continuous(
    breaks = c(1, 25, 50, 75, 100),  # Set the x-axis breaks
    labels = c(1, 25, 50, 75, 100)  # Subscript formatting
  ) +
  scale_color_manual(values = c("Mkt-RF" = "blue", "SMB" = "green", "HML" = "red", "RMW" = "purple", "CMA" = "orange")) +
  geom_hline(yintercept = 0, color = "#707070", linetype = "dashed")  





