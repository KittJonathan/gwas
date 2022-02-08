# Linear Models with R

# Julian J. Faraway

# http://www.utstat.toronto.edu/~brunner/books/LinearModelsWithR.pdf

# Last updated 2022-02-07

# Chapter 1 - Introduction ----

library(faraway)
data(pima)
pima
summary(pima)
sort(pima$diastolic)

pima$diastolic[pima$diastolic == 0] <- NA
pima$glucose[pima$glucose == 0] <- NA
pima$triceps[pima$triceps == 0] <- NA
pima$insulin[pima$insulin == 0] <- NA
pima$bmi[pima$bmi == 0] <- NA

pima$test <- factor(pima$test)
summary(pima)
