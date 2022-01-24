# Linear mixed models

# Last updated 2022-01-24

# Links ----

# https://ourcodingclub.github.io/tutorials/mixed-models/
# https://towardsdatascience.com/how-linear-mixed-model-works-350950a82911
# https://poissonisfish.com/2017/12/11/linear-mixed-effect-models-in-r/

# Introduction to linear mixed models ----

# https://ourcodingclub.github.io/tutorials/mixed-models/

# Explore the data (3 sites in 8 different mountain ranges)

load("data/CC-Linear-mixed-models-master/dragons.RData")
head(dragons)
hist(dragons$testScore)  # seems close to a normal distribution

dragons$bodyLength2 <- scale(dragons$bodyLength,  # center (mean = 0) and scale (sd = 1) 
                             center = TRUE,
                             scale = TRUE)

# Fit all data in one analysis

