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

## Create the age_at_death_density plot (rotated 90 degrees and spanning 3.5 to 1.5 on the y-axis)
x_min <- 0
x_max <- 1.5
x <- age_at_death_density
x_vect <- x_min + ((x - min(x)) * (x_max - x_min)) / (max(x) - min(x))
y_vect <- seq(3.35, 1.05, length.out = length(age_at_death_density))

age_death_plot <- ggplot() +
  geom_line(aes(x = x_vect, y = y_vect)) +
  xlim(x_min, x_max) + ylim(0.5, 4.5) + 
  coord_fixed(ratio = 1) +
  theme_void()

# Make the joint probability from the two marginals
joint_probability <- outer(pop_size_marginal, rev(age_marginal))

# Create the population_size plot
population_plot <- ggplot() +
  geom_line(aes(x = tau, y = population_size)) +
  theme_void()

# Create a text grid with shape 1 x 4 for the population size marginal data
pop_size_marginal_df <- data.frame(
  x = 1:4,
  y = rep(4, 4),
  label = round(pop_size_marginal, 2)
)

label_df <- data.frame(
  x = seq(0.5, 4.5, length.out = 5),  # Positions for 5 labels
  y = rep(4.75, 5),  # Constant y position for all labels
  label = seq(600, 1000, by = 100)  # Labels from 600 to 1000
)

axis_df <- data.frame(
  x = c(2.5),
  y = c(5.25),
  label = 'Calendar Date [AD]'
)

pop_size_marginal_grid_plot <- ggplot(pop_size_marginal_df, aes(x = x, y = y, label = label)) +
  geom_tile(color = "black", fill = "white", size = 0.5) +  # Draw squares with boundaries
  geom_text() +  # Add text labels within the squares
  geom_text(data = label_df, size=3, aes(x = x, y = y, label = label), vjust = 0) +  # Add year labels above the grid lines using separate data frame
  geom_text(data = axis_df, size=6, aes(x = x, y = y, label = label), vjust = 0) +  # Add an axis label
  theme_void() +  # Removes axes, labels, and background
  xlim(0.5, 4.5) + ylim(0.5, 5.5) +  # Adjust limits to fit the grid and labels
  coord_fixed(ratio = 1)  # Set aspect ratio to 1:1

# Transform the joint probability matrix into a data frame
joint_prob_df <- expand.grid(x = 1:num_parts, y = 1:num_parts)
joint_prob_df$label <- sprintf("%.2f", as.vector(joint_probability))

# Create a 4x4 grid plot for the joint probability data
age_label_df <- data.frame(
  x = rep(4.5, 5),  # Position to the right of the grid
  y = 0.5:4.5,  # Y positions corresponding to the grid
  label = c("0", "20", "40", "60", "80")  # Age labels
)

age_axis_label_df <- data.frame(
  x = 4.75,  # Position further right from the grid
  y = 2.5,   # Middle position on the y-axis
  label = "Age [years]"
)
joint_probability_grid_plot <- ggplot(joint_prob_df, aes(x = x, y = y, label = label)) +
  geom_tile(color = "black", fill = "white", size = 0.5) +  # Draw squares with boundaries
  geom_text() +  # Add text labels
  theme_void() +  # Removes axes, labels, and background
  geom_text(data = label_df, size=3, aes(x = x, y = y, label = label), vjust = 0) +  # Add year labels above the grid lines using separate data frame
  geom_text(data = axis_df, size=6, aes(x = x, y = y, label = label), vjust = 0) +  # Add an axis label
  xlim(0.5, 4.5) + ylim(0.5, 5.5) +  # Set limits to enclose the squares properly
  coord_fixed(ratio = 1)  # Set aspect ratio to 1:1
joint_probability_grid_plot <- joint_probability_grid_plot +
    geom_text(data = age_label_df, aes(x = x, y = y, label = label), angle=-90, vjust = -0.5, hjust=0.5, size = 3) 
    #geom_text(data = age_label_df, aes(x = x, y = y, label = label), angle = -90, vjust = 0.0, hjust = 0, size = 3) 
#  geom_text(data = age_axis_label_df, aes(x = x, y = y, label = label), vjust = 0.5, size = 6)  # Add age axis label


# Create a text plot for cell (1,1)
text_plot1 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Model 1\nRadiocarbon Samples\n\n", 
           size = 5, hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "pt"))

# Create a text plot for cell (1,3)
text_plot2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Model 2\nRadiocarbon Samples\n+\nSkeletal Age at Death", 
           size = 5, hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "pt"))

# Creating spacers and empty plots
empty_plot <- plot_spacer()

# Arrange the plots manually
widths <- c(2, 1, 2, 2)  # Adjust the width of each column
heights <- c(1, 1, 2)  # Adjust the height of each row
layout <- text_plot1                  + empty_plot + text_plot2                  + empty_plot +
          population_plot             + empty_plot + population_plot             + empty_plot +
          pop_size_marginal_grid_plot + empty_plot + joint_probability_grid_plot + age_death_plot +
          plot_layout(widths = widths, heights = heights)


# Define the file name and path for the output
output_file <- "plot_layout.png"

# Save the plot to a PNG file with high DPI
ggsave(output_file, layout, width = 10, height = 10, dpi = 300)