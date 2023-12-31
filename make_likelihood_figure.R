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
w1 = 0.5
w2 = 1 - w1
population_size <- w1 * dnorm(tau, mean=700, sd=25) + w2*dnorm(tau, mean=850, sd=50)

# Create discretized probability masses for the population marginal distribution
tau_marginal = seq(start_date, start_date + num_times, by=100)
pop_size_marginal <- diff(w1 * pnorm(tau_marginal, mean=700, sd=25) + w2*pnorm(tau_marginal, mean=850, sd=50))
pop_size_marginal <- pop_size_marginal / sum(pop_size_marginal)

# Define parameters for age_at_death_density data
age_groups <- 0:80  # Age groups starting at 0

# Decreasing exponential for age at death
lambda_ <- 0.05  # Rate parameter for exponential distribution
age_at_death_density <- exp(-lambda_ * age_groups)

# Build the age_marginal from the longer vector with four cells
num_parts <- 4
age_marginal <- rep(NA, num_parts)
indices_per_part <- (length(age_groups) - 1) / num_parts

for (i in 1:num_parts) {
  start_index <- (i - 1) * indices_per_part + 1
  end_index <- i * indices_per_part
  age_marginal[i] <- sum(age_at_death_density[start_index:end_index])
}
age_marginal <- age_marginal / sum(age_marginal)

# Make the joint probability from the two marginals
joint_probability <- outer(pop_size_marginal, age_marginal)

# Create the population_size plot
population_plot <- ggplot() +
  geom_line(aes(x = tau, y = population_size)) +
  theme_void()  # Removes axes, labels, and background

# Create the age_at_death_density plot (rotated 90 degrees)
age_death_plot <- ggplot() +
  geom_line(aes(x = age_at_death_density, y = -age_groups)) +
  theme_void()  # Removes axes, labels, and background

# Create a text grid with shape 1 x 4 for the population size marginal data
# Create a grid plot with squares for the population size marginal data
pop_size_marginal_df <- data.frame(
  x = 1:4,
  y = rep(1, 4),
  label = round(pop_size_marginal, 2)
)

pop_size_marginal_grid_plot <- ggplot(pop_size_marginal_df, aes(x = x, y = y, label = label)) +
  geom_tile(color = "black", fill = "white", size = 0.5) +  # Draw squares with boundaries
  geom_text() +  # Add text labels
  theme_void() +  # Removes axes, labels, and background
  xlim(0.5, 4.5) + ylim(0.5, 1.5) +  # Set limits to enclose the squares properly
  coord_fixed(ratio = 1)  # Set aspect ratio to 1:1

# Creating spacers and empty plots
empty_plot <- plot_spacer()

# Arrange the plots manually
widths <- c(2, 1, 2, 1)  # Adjust the width of each column
heights <- c(1, 1, 2)  # Adjust the height of each row
layout <- empty_plot                  + empty_plot + empty_plot      + empty_plot +
          population_plot             + empty_plot + population_plot + empty_plot +
          pop_size_marginal_grid_plot + empty_plot + empty_plot      + age_death_plot +
          plot_layout(widths = widths, heights = heights)

# Define the file name and path for the output
output_file <- "plot_layout.png"

# Save the plot to a PNG file with high DPI
ggsave(output_file, layout, width = 10, height = 10, dpi = 300)