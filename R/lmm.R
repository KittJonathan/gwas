# Linear mixed models

# Last updated 2022-01-27

# Links ----

# https://ourcodingclub.github.io/tutorials/mixed-models/
# https://towardsdatascience.com/how-linear-mixed-model-works-350950a82911
# https://poissonisfish.com/2017/12/11/linear-mixed-effect-models-in-r/
# https://mfviz.com/hierarchical-models/
# https://onlinelibrary.wiley.com/doi/full/10.1046/j.1439-037X.2003.00049.x

# Load packages ----

library(lme4)
library(tidyverse)

# Linear models and linear mixed effects models in R: Tutorial 1 ----

# https://bodo-winter.net/tutorials.html

# Create the dataset
pitch <- c(233, 204, 242, 130, 112, 142)
sex <- c(rep("female", 3), rep("male", 3))
my.df <- data.frame(sex, pitch)
my.df

# Linear model
xmdl <- lm(pitch ~ sex, my.df)
summary(xmdl)

mean(my.df[my.df$sex == "female", ]$pitch)
mean(my.df[my.df$sex == "male", ]$pitch)

# Other dataset
age <- c(14, 23, 35, 48, 52, 67)
pitch <- c(252, 244, 240, 233, 212, 204)
my.df <- data.frame(age, pitch)
xmdl <- lm(pitch ~ age, my.df)
summary(xmdl)

# Substract the mean age from each age value
my.df$age.c <- my.df$age - mean(my.df$age)
xmdl <- lm(pitch ~ age.c, my.df)
summary(xmdl)

# Construct the residual plot
plot(fitted(xmdl), residuals(xmdl))
