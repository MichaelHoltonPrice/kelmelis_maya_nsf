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



y_vect <- seq(3.5, 1.5, length.out = length(age_at_death_density))

# Create breaks for the y-axis
y_breaks <- c(0, 20, 40, 60, 80)
# Map these breaks to the corresponding y_vect values
mapped_y_breaks <- approx(age_groups, y_vect, y_breaks)$y

age_death_plot <- ggplot() +
  geom_line(aes(x = age_at_death_density, y = y_vect)) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text.x = element_blank(),  # Remove x-axis text if not needed
    axis.ticks.x = element_blank(),  # Remove x-axis ticks if not needed
    panel.grid = element_blank(),   # Remove grid lines
    axis.text.y = element_text(angle = -90)  # Rotate y-axis labels by 90 degrees
  ) +
  scale_y_continuous(breaks = mapped_y_breaks, labels = as.character(y_breaks))






## Create the age_at_death_density plot (rotated 90 degrees and spanning 3.5 to 1.5 on the y-axis)
#y_vect <- seq(3.5, 1.5, len=length(age_at_death_density)) # This replaces the age vector
#y_breaks <- seq(3.5, 1.5, len=5)
##age_death_plot <- ggplot() +
##  geom_line(aes(x = age_at_death_density, y = y_vect)) +
##  theme_void() +  # Removes axes, labels, and background
##  ylim(0.5, 4.5) +
##  coord_fixed(ratio = 1)
#
#age_death_plot <- ggplot() +
#  geom_line(aes(x = age_at_death_density, y = y_vect)) +
#  #ylim(0.5, 4.5) +
#  coord_fixed(ratio = 1) +
#  theme_minimal() +  # Use a minimal theme
#  theme(
#    axis.title = element_blank(),  # Remove axis titles
#    axis.text.y = element_blank(),  # Remove y-axis text
#    axis.ticks.y = element_blank(),  # Remove y-axis ticks
#    panel.grid = element_blank()  # Remove grid lines
#  ) +
#  scale_y_continuous(breaks = y_breaks,
#                     labels = c("0", "20", "40", "60", "80"))  # Labels for the breaks



# Make the joint probability from the two marginals
joint_probability <- outer(pop_size_marginal, rev(age_marginal))

# Create the population_size plot
population_plot <- ggplot() +
  geom_line(aes(x = tau, y = population_size)) +
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    panel.grid = element_blank()  # Remove grid lines
  ) +
  scale_x_continuous(breaks = seq(600, 1000, by = 100),    # Set custom breaks at 100-year intervals
                     labels = seq(600, 1000, by = 100)) +  # Labels for the breaks
  labs(x = "Calendar Date [AD]") +
  theme(axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = -25, b = 10)))


# Create a text grid with shape 1 x 4 for the population size marginal data
pop_size_marginal_df <- data.frame(
  x = 1:4,
  y = rep(4, 4),
  label = round(pop_size_marginal, 2)
)

pop_size_marginal_grid_plot <- ggplot(pop_size_marginal_df, aes(x = x, y = y, label = label)) +
  geom_tile(color = "black", fill = "white", size = 0.5) +  # Draw squares with boundaries
  geom_text() +  # Add text labels
  theme_void() +  # Removes axes, labels, and background
  xlim(0.5, 4.5) + ylim(0.5, 4.5) +  # Adjust limits to fit the 4x4 grid
  coord_fixed(ratio = 1)  # Set aspect ratio to 1:1

# Transform the joint probability matrix into a data frame
joint_prob_df <- expand.grid(x = 1:num_parts, y = 1:num_parts)
joint_prob_df$label <- sprintf("%.2f", as.vector(joint_probability))

# Create a 4x4 grid plot for the joint probability data
joint_probability_grid_plot <- ggplot(joint_prob_df, aes(x = x, y = y, label = label)) +
  geom_tile(color = "black", fill = "white", size = 0.5) +  # Draw squares with boundaries
  geom_text() +  # Add text labels
  theme_void() +  # Removes axes, labels, and background
  xlim(0.5, num_parts + 0.5) + ylim(0.5, num_parts + 0.5) +  # Set limits to enclose the squares properly
  coord_fixed(ratio = 1)  # Set aspect ratio to 1:1

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

# Create a plot with a vertical line
#vertical_line_plot <- ggplot() + 
#  geom_vline(xintercept = 0.5, color = "black", size = 0.5) + 
#  theme_void() +
#  theme(plot.margin = margin(0, 0, 0, 0)) +
#  coord_fixed(ratio = 1)

# Creating spacers and empty plots
empty_plot <- plot_spacer()

# Arrange the plots manually
widths <- c(2, 1, 2, 1)  # Adjust the width of each column
heights <- c(1, 1, 4)  # Adjust the height of each row
#layout <- text_plot1                  + vertical_line_plot + text_plot2                  + empty_plot +
#          population_plot             + vertical_line_plot + population_plot             + empty_plot +
#          pop_size_marginal_grid_plot + vertical_line_plot + joint_probability_grid_plot + age_death_plot +
#          plot_layout(widths = widths, heights = heights)
layout <- text_plot1                  + empty_plot + text_plot2                  + empty_plot +
          population_plot             + empty_plot + population_plot             + empty_plot +
          pop_size_marginal_grid_plot + empty_plot + joint_probability_grid_plot + age_death_plot +
          plot_layout(widths = widths, heights = heights)


# Define the file name and path for the output
output_file <- "plot_layout.png"

# Save the plot to a PNG file with high DPI
ggsave(output_file, layout, width = 10, height = 10, dpi = 300)