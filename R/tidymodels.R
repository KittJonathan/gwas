# Tidymodels

# Last updated 2022-02-08

# Load packages ----

library(tidymodels)
library(tidyverse)
library(readr)
library(broom.mixed)
library(dotwhisker)

# Get started ----

# https://www.tidymodels.org/start/models/

# Import the sea urchins data

urchins <- read_csv("https://tidymodels.org/start/models/urchins.csv") %>% 
  setNames(c("food_regime", "initial_volume", "width")) %>% 
  mutate(food_regime = factor(food_regime, levels = c("Initial", "Low", "High")))

# Quick look at the data

urchins

# Plot the data

ggplot(data = urchins,
       mapping = aes(x = initial_volume,
                     y = width,
                     group = food_regime,
                     col = food_regime)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_colour_viridis_d(option = "plasma", end = .7) +
  theme_minimal()
