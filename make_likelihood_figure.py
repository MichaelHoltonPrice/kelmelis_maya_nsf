import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes

# Define parameters
num_times = 400
start_date = 600  # AD 600
tau = np.arange(start_date, start_date + num_times)  # Calendar dates
age_groups = np.arange(0, 81)  # Age groups starting at 0

# Mixture of Gaussians for calendar dates
population_size = norm.pdf(tau, loc=700, scale=25) + norm.pdf(tau, loc=850, scale=50)

# Decreasing exponential for age at death
lambda_ = 0.05  # Rate parameter for exponential distribution
age_at_death_density = np.exp(-lambda_ * age_groups)

# Set the figure size and layout
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()

# Plot for Mixture of Gaussians for Calendar Dates
ax1 = fig.add_subplot(121)
ax1.plot(tau, population_size)
ax1.set_title("Mixture of Gaussians for Calendar Dates")
ax1.set_xlabel("Calendar Years")
ax1.set_ylabel("Density")

# Prepare to rotate the age-at-death density plot by -90 degrees
# Adjust extremes to match the data range after rotation
extremes = (0, max(age_at_death_density), 0, max(age_groups))

transform = Affine2D().rotate_deg(-90)
helper = floating_axes.GridHelperCurveLinear(transform, extremes=extremes)

# Add floating axes for the rotated plot
ax2 = floating_axes.FloatingSubplot(fig, 122, grid_helper=helper)
fig.add_subplot(ax2)

# Plotting on the transformed axes
ax2 = ax2.get_aux_axes(transform)
ax2.plot(age_at_death_density, age_groups)

# Adjust layout
plt.tight_layout()

# Save to PDF
plt.savefig("likelihood_plot.pdf")