# Linear mixed models

# Last updated 2022-01-28

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

# Simulating equal variances
plot(rnorm(100), rnorm(100))

# Normality of residuals
hist(residuals(xmdl))
qqnorm(residuals(xmdl))

# Absence of influential data points
dfbeta(xmdl)

# A very basic tutorial for performing linear mixed effects analyses: Tutorial 2 ----

# Load data
politeness <- read.csv("http://www.bodowinter.com/tutorial/politeness_data.csv")

# Look at the data
head(politeness)
tail(politeness)
summary(politeness)
str(politeness)
colnames(politeness)
which(is.na(politeness$frequency))
which(!complete.cases(politeness))

# Relationship between politeness & pitch
boxplot(frequency ~ attitude*gender,
        col = c("white", "lightgray"), data = politeness)

# Construct the model
politeness.model <- lmer(frequency ~ attitude + (1|subject) + (1|scenario),
                         data = politeness)

summary(politeness.model)

# Add gender as an additional fixed effect
politeness.model <- lmer(frequency ~ attitude + gender + (1|subject) + (1|scenario),
                         data = politeness)
summary(politeness.model)

# Construct the null model
politeness.null <- lmer(frequency ~ gender + (1|subject) + (1|scenario),
                        data = politeness, REML = FALSE)

# Re-do the full model, with REML = FALSE
politeness.model <- lmer(frequency ~ attitude + gender + (1|subject) + (1|scenario),
                         data = politeness, REML = FALSE)

# Perform the likelihood ratio test
anova(politeness.null, politeness.model)

# Random slopes versus random intercepts
coef(politeness.model)

# Create a model which allows random slopes
politeness.model <- lmer(frequency ~ attitude + gender + 
                           (1+attitude|subject) + (1+attitude|scenario),
                         data = politeness, REML = FALSE)
coef(politeness.model)

# Create a random slope null model
politeness.null <- lmer(frequency ~ gender + 
                          (1+attitude|subject) + (1+attitude|scenario),
                        data = politeness, REML = FALSE)

# Likelihood ratio test
anova(politeness.null, politeness.model)
