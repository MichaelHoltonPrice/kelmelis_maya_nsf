# 0. Clear the workspace, load packages, source the needed functions, create a
#    data directory (if necessary), and create an analysis name based on the
#    current time.
rm(list=ls())
library(baydem)
library(matrixStats) # for logSumExp
library(HDInterval)
library(doParallel)

source("power_analysis_functions.R")

if (!dir.exists('data')) {
  dir.create('data')
}

analysis_name <- as.character(Sys.time())
analysis_name <- gsub(":", "-", analysis_name)

# 1. Set the true values

# The true total fertiltiy rate is 5 (half of babies are female)
TFR0 <- 5

# The true Siler parameter vector is from Gage and Dyke (1986), table TBD
a0 <- c(0.175,1.40,0.368*0.01,0.075*0.001,0.917*0.1)

# The true boost to mortality is kappa = .2 (a 20% increase)
kappa0 <- .2

# The final model parameter is the probability of transitioning from a
# non-famine to a famine year, which is half the probabilty of transitioning
# from a famine to a famine year. We call this parameter p_gb (or p_gb0 for the
# true, simulated value). The transition matrix is
#
# [1 -   p_gb,   p_gb]
# [1 - 2*p_gb, 2*p_gb]
#
# where the first state is the good (non-famine) state. Put slightly differently,
# the individual transition probabilities are
#
# non-famine -->     famine   p_gb
# non-famine --> non-famine   1 - p_gb
#     famine -->     famine   2*p_gb
#     famine --> non-famine   1 - 2*p_gb
p_gb <- .2 # TODO: this should be called p_gb0 for consistency with the other true values

# 2. Run the experiments/simulations
# Load the intcal20 calibration curve and define the error model
# for radiocarbon samples
calib_df <- load_calib_curve('intcal20')
error_spec=list(type="unif_fm",min=.0021,max=.0028)

# Our study period is 400 years, and each time period is one year,
# so there are 400 time periods
num_times <- 400
# Our simulated study period starts in AD 600
start_date <- 600 # AD 600
# Define the vector of calendar dates for the study period
tau <- seq(start_date, start_date + num_times -1)

# N_vect is the vector of samples sizes. The power calculation involves
# determing the number of times each experiment/simulation succeeds for
# each value of the sample size. See the method exp_wrapper in
# power_analysis_functions.R for what constitutes success. For the
# simulations that include age-at-death observations, we utilize N/10
# skeletons with an age-at-death estimate for every radiocarbon samples.
N_vect <- 100*2^(0:5)

# We run 400 simulations for each value of the sample size so that we
# can estimate the success probability.
exp_per_N <- 400

# Set seeds for reproducibility
set.seed(1011)
seed_matrix <- matrix(sample.int(length(N_vect)*exp_per_N),nrow=length(N_vect))

# Loop over experiments to create, for each, the problem needing to be solved.
# This involves three nested for loops: (a) samples sizes in N_vect, (b) 400
# experiments per sample size, and (c) with and without age-at-death observations.
counter <- 0
prob_list <- list()
for (k1 in 1:length(N_vect)) {
  N <- N_vect[k1]
  for (k2 in 1:exp_per_N) {
    for (k3 in 1:2) {
      use_age <- k3 == 2

      counter <- counter + 1
      prob <- list(p_gb=p_gb,
                   a0=a0,
                   kappa0=kappa0,
                   TFR0=TFR0,
                   N=N,
                   num_times=num_times,
                   use_age=use_age,
                   seed=seed_matrix[k1,k2],
                   counter=counter,
                   num_prob=length(N_vect)*exp_per_N*2,
                   analysis_name=analysis_name)
      prob_list[[counter]] <- prob
    }
  }
}

t0 <- Sys.time()

# Call the wrapper that runs experiments inside a parallel for loop
package_list <- c('baydem', 'matrixStats', 'HDInterval', 'doParallel')
registerDoParallel(detectCores())
t0 <- Sys.time()
success_vect <- foreach(n=1:length(prob_list),.combine=cbind, .packages=package_list) %dopar% {
  success <- exp_wrapper(prob_list[[n]])
}
t1 <- Sys.time()
stopImplicitCluster()

# 3. Plots the results
# Now that all experiments have finished, extract the success metrics. 
success_array <- array(NA,dim=c(length(N_vect), exp_per_N, 2))
counter <- 0
for (k1 in 1:length(N_vect)) {
  N <- N_vect[k1]
  for (k2 in 1:exp_per_N) {
    for (k3 in 1:2) {
      counter <- counter + 1
      success_array[k1, k2, k3] <- success_vect[counter]
    }
  }
}

# Plot the success probabilities (statistical power)
power_vect_no_age <- rowSums(success_array[,,1]) / ncol(success_array[,,1])
power_vect_age    <- rowSums(success_array[,,2]) / ncol(success_array[,,2])

y_max <- max(power_vect_no_age, power_vect_age)
pdf('power_curve.pdf')
  plot(N_vect, power_vect_no_age, xlab='N', ylab='Power', ylim=c(0,y_max), col='black')
  points(N_vect, power_vect_age, col='red')
dev.off()