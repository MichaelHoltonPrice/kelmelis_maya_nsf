# Load necessary libraries
library(ggplot2)
library(patchwork)

# Clear the workspace
rm(list=ls())

# Define parameters for population_size data
num_times <- 400
start_date <- 600  # AD 600
tau <- seq(start_date, start_date + num_times - 1)  # Calendar dates

# Mixture of Gaussians for calendar dates
population_size <- dnorm(tau, mean=700, sd=25) + dnorm(tau, mean=850, sd=50)

# Define parameters for age_at_death_density data
age_groups <- 0:80  # Age groups starting at 0

# Decreasing exponential for age at death
lambda_ <- 0.05  # Rate parameter for exponential distribution
age_at_death_density <- exp(-lambda_ * age_groups)

# Creating the population_size plot with only the curve
population_plot <- ggplot() +
  geom_line(aes(x = tau, y = population_size)) +
  theme_void()  # Removes axes, labels, and background

# Creating the age_at_death_density plot rotated by 90 degrees (swapping axes)
age_death_plot <- ggplot() +
  geom_line(aes(x = age_at_death_density, y = -age_groups)) +
  theme_void()  # Removes axes, labels, and background

# Creating spacers and empty plots
empty_plot <- plot_spacer()

# Arrange the plots manually
widths <- c(2, 1, 2, 1)  # Adjust the width of each column
heights <- c(1, 1, 2)  # Adjust the height of each row
layout <- empty_plot      + empty_plot + empty_plot      + empty_plot +
          population_plot + empty_plot + population_plot + empty_plot +
          empty_plot      + empty_plot + empty_plot      + age_death_plot +
          plot_layout(widths = widths, heights = heights)

# Define the file name and path for the output
output_file <- "plot_layout.png"

# Save the plot to a PNG file with high DPI
ggsave(output_file, layout, width = 10, height = 10, dpi = 300)