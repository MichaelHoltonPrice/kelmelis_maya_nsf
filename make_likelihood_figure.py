import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Define parameters
num_times = 400
start_date = 600 # AD 600
tau = np.arange(start_date, start_date + num_times) # Calendar dates
age_groups = np.arange(0, 81) # Age groups starting at 0

# Mixture of Gaussians for calendar dates
population_size = norm.pdf(tau, loc=700, scale=25) + norm.pdf(tau, loc=850, scale=50)

# Decreasing exponential for age at death
lambda_ = 0.05 # Rate parameter for exponential distribution
age_at_death_density = np.exp(-lambda_ * age_groups)

# Create figure and axes
fig = plt.figure(figsize=(10, 8))

# Top plot (Mixture of Gaussians for Calendar Dates)
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(tau, population_size)
ax1.set_title("Mixture of Gaussians for Calendar Dates")
ax1.set_xlabel("Calendar Years")
ax1.set_ylabel("Density")

# Right plot (Density for Age at Death)
ax2 = fig.add_subplot(2, 2, 4)
ax2.plot(age_at_death_density, age_groups)
ax2.set_xlabel("Density")
ax2.set_ylabel("Age at Death")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

# Adjust layout
plt.tight_layout()

# Save to PDF
plt.savefig("likelihood_plot.pdf")