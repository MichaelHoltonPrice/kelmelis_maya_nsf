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
f0 <- 5

# The true Siler parameter vector is from Gage and Dyke (1986), table TBD
a0 <- c(0.175,1.40,0.368*0.01,0.075*0.001,0.917*0.1)

# The true boost to mortality is kappa = .2 (a 20% increase)
kappa0 <- .2

# The final model parameter is the probability of transitioning from a
# non-famine to a famine year, which is half the probabilty of transitioning
# from a famine to a famine year. This parameter is the effect size used in
# the power analysis, and hence is not set here.

# For reproducibility, set the random number
set.seed(1000)

# 2. Iterate over effect size to calculate the power (0.025 to 0.4 by .025)

# p_gb is the effect size and stands for the probability of good to bad
# (non-famine to famine year). That is, we assume a hidden markov model for the
# un-observed famine state (that is, not directly observed; it is indirectly
# inferred), where the transition matrix is
#
# [1 -   p_gb,   p_gb]
# [1 - 2*p_gb, 2*p_gb]
#
# and the first state is the good (non-famine) state. That is, the individual
# transition probabilities are
#
# non-famine -->     famine   p_gb
# non-famine --> non-famine   1 - p_gb
#     famine -->     famine   2*p_gb
#     famine --> non-famine   1 - 2*p_gb



# First, define a function to calculate the population size given the input
# parameter vectors, including the realized famine years. To do so, utilize a
# discrete, age-structured population with one year age intervals and one year
# time steps.

calib_df <- load_calib_curve('intcal20')
error_spec=list(type="unif_fm",min=.0021,max=.0028)

num_times <- 400
start_date <- 600 # AD 600 is the first date
tau <- seq(start_date, start_date + num_times -1)
# The number of experiments for each value of the effect size


p_gb <- .2
N_vect <- 100*2^(0:5)
#N_vect <- N_vect[3:5]
#N_vect <- c(100,1000)
exp_per_N <- 100

set.seed(1011)
seed_matrix <- matrix(sample.int(length(N_vect)*exp_per_N),nrow=length(N_vect))

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
                   f0=f0,
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

registerDoParallel(detectCores())
t0 <- Sys.time()
success_vect <- foreach(n=1:length(prob_list),.combine=cbind) %dopar% {
  success <- exp_wrapper(prob_list[[n]])
}
t1 <- Sys.time()
stopImplicitCluster()

# For clarity, use another set of for loops to unpack success_vect
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

power_vect_no_age <- rowSums(success_array[,,1]) / ncol(success_array[,,1])
power_vect_age    <- rowSums(success_array[,,2]) / ncol(success_array[,,2])

y_max <- max(power_vect_no_age, power_vect_age)
pdf('power_curve.pdf')
  plot(N_vect, power_vect_no_age, xlab='N', ylab='Power', ylim=c(0,y_max), col='black')
  points(N_vect, power_vect_age, col='red')
dev.off()


#exp_obj_no_age <- run_experiment(p_gb, a0, kappa0, f0, N, num_times, 100, FALSE)
#exp_obj_age <- run_experiment(p_gb, a0, kappa0, f0, N, num_times, 100, TRUE)