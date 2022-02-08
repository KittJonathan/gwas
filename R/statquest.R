# StatQuest videos

# Last updated 2022-02-07

# Linear Regression in R ----

# https://www.youtube.com/watch?v=u1cc1r_Y7M0

mouse.data <- data.frame(
  weight = c(0.9, 1.8, 2.4, 3.5, 3.9, 4.4, 5.1, 5.6, 6.3),
  size = c(1.4, 2.6, 1.0, 3.7, 5.5, 3.2, 3.0, 4.9, 6.3))

mouse.data

plot(mouse.data$weight, mouse.data$size)

mouse.regression <- lm(size ~ weight, data = mouse.data)
