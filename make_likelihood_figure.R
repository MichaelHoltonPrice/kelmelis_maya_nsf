# Load necessary libraries
library(ggplot2)
library(patchwork)

# Clear the workspace
rm(list=ls())

# Define parameters for your data
num_times <- 400
start_date <- 600  # AD 600
tau <- seq(start_date, start_date + num_times - 1)  # Calendar dates

# Mixture of Gaussians for calendar dates
population_size <- dnorm(tau, mean=700, sd=25) + dnorm(tau, mean=850, sd=50)

# Creating the population_size plot with only the curve
population_plot <- ggplot() +
  geom_line(aes(x = tau, y = population_size)) +
  theme_void()  # Removes axes, labels, and background

# Creating spacers and empty plots
empty_plot <- plot_spacer()

# Arrange the plots manually
layout <- (empty_plot + empty_plot + empty_plot) /
          (population_plot + population_plot + empty_plot) /
          (empty_plot + empty_plot + empty_plot)

# Define the file name and path for the output
output_file <- "plot_layout.png"

# Save the plot to a PNG file with high DPI
ggsave(output_file, layout, width = 10, height = 10, dpi = 300)