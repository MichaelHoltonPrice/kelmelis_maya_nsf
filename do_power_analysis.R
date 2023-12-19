# This script runs the power analysis for the NSF proposal by Kelmelis et al.
# named:
#
# A new end-to-end Bayesian and osteological approach to reconstructing
# demographic changes for the ancient Maya
#
# DEMOGRAPHIC MODEL
#
# The target of the power analysis is a simuluated boost to mortality during
# famine years of kappa = 0.2. We assume a baseline Siler mortality hazard of
#
# haz_good = a_1 * exp(-a_2 * x) + a_3 + a_4 * exp(a_5 * x)
#
# where we use a0 <- c(0.175,1.40,0.368*0.01,0.075*0.001,0.917*0.1) for
# the simulated, baseline hazard. This baseline hazard is for good/non-famine
# years. During a famine, the hazard is increased by (1 + kappa) (20% higher)
#
# haz_bad = (1 + kappa) * haz_good
#
# To fully specify the simulated demography we must also define a fertility
# model. We use a simple model where (a) age-specific, instantaneous fertility
# is constant from ages 15 to 45 and 0 otherwise and (b) the total fertility
# rate (TFR) is 5. For details (e.g., what precisely we mean by instantaneous
# fertility), see the documentation for the function calc_F in the file
# power_analysis_functions.R.
#
# HIDDEN MARKOV MODEL
#
# The other major component of the model is a hidden Markov model that
# determines which years are famine years and which are non-famine years.
# This requires just one parameter, p_gb, which is the probability of
# transitioning from a non-famine to a famine year (_gb in p_gb marks that this
# is the probability of going from a good to bad year). We further assume
# that the probability of going from a bad year to another bad year is 2*p_gb.
# That is, once a bad year is reached one gets "stuck" in it. This adds some
# realism since there will, on average, be clusters of good and bad years. The
# transition matrix, which we usually call H in the code, is thus (applying
# conversation of probabilities):
#
# [1 -   p_gb,   p_gb]
# [1 - 2*p_gb, 2*p_gb]
#
# Broken out and using some English, the trasitions are governed by:
#
# non-famine -->     famine   p_gb
# non-famine --> non-famine   1 - p_gb
#     famine -->     famine   2*p_gb
#     famine --> non-famine   1 - 2*p_gb
#
# Associated with these transition probabilities are actual, realized
# sequences of good and bad years. In the code, we call these realized sequences
# h or famine_mask (a boolean vector). For example, for each
# experiment/simulation that contributes to the power calculation there is a
# true, realized famine_mask. If the mask begins [FALSE, FALSE, TRUE,...], then
# the first two years (AD 800 and 801) are good/non-famine years and the third
# year (AD 802) is a bad/famine year. Given a realized famine mask sequence h
# and the demographic parameter vector theta = [a, kappa, TFR], one can fully
# calculate the realized demography, which is a matrix indexed by year and age,
# where each entry in the matrix gives the probability of dying in that given
# year at that given age (at least, this is the right conceptual mindset; the
# details matter and the code takes advantage of the fact that we actually only
# need the marginal probabilities). For details, see the documentation for the
# function calc_pop_size in the file power_analysis_functions.R.
#
# POWER CALCULATION
#
# The target of the power calculation, as already noted above, is the parameter
# kappa, which is the "boost" to mortality during bad/famine years. For each
# value of the sample size, N, we run 400 simulated experiments. An experiment
# passes if the 95% highest density interval (HDI) from Bayesian samples of the
# estimated value of kappa falls entirely within the region of practical
# equivalence (ROPE). It makes sense that we require the reconstructed kappa
# to be positive (since we want to assess whether we can correctly identify
# boosts to mortality that results from famine years). Hence, we use 0 as the
# lower limit of the ROPE. To make things symmetric, we define the upper limit
# of the ROPE to be 0.4, so that values of the HDI between 0.2 +/- 0.2 qualify
# as a success.
#
# The other major parameter in the power calculation is a flag called use_age.
# If use_age is TRUE, the reconstruction may use simulated skeletal
# age-at-death estimates. We assume fewer estimates of this type are available
# compared to radiocarbon dates. In particular, we use 10% of the radiocarbon
# sample size for the skeletal sample size, rounded to the nearest integer. For
# example, if there are 2000 radiocarbon samples there are 200 skeletal samples.
#
# Ultimately, the power calculation yields two curves: the success probability
# as a function of sample size both with and without skeletal samples.
#
# EXPERIMENTS
# 
# Each experiment consists of the following steps:
#
# (1) First, create simulated data. This uses the true/simulated values for
#     a, kappa, TFR, and p_gb. The simulated dataset includes a true/realized
#     sequence of good and bad years, true/realized sample calendar dates
#     and ages-at-death, and realized radiocarbon measurements that account
#     for the radiocarbon calibration curve and measurement uncertainty.
#
# (2) Second, do posterior Bayesian sampling of the parameters a, kappa, TFR,
#     and p_gb. Crucially, only quantities available in practice are used for
#     this Bayesian inference; for example, the realized sequence of good and
#     bad years is not used and only the radiocarbon measurements are
#     available, not the realized calendar dates. In addition to sampling the
#     preceding parameters, we also sample the estimated famine mask. Below we
#     provide further details on the posterior Bayesian sampling.
#
# (3) Third, aggregate and plot the experimental results. An experiment passes
#     if the HDI for the posterior samples of kappa lies fully within the
#     ROPE, which spans 0 to 0.4.
#
# BAYESIAN SAMPLING
#
# We utilize an iterative sampling with three steps involving the following
# parameters:
#
#  theta = th = [a_1, a_2, a_3, a_4, a_5, kappa, TFR]
#  m          = [p_gb]
#  h          = The estimated famine mask, or sequence of good and bad years
#
# Step 1: Sample h|m (h given m)
#
# The parameter p_gb completely specifies the Hidden markov model for good and
# bad years. We draw a novel sample of the realized sequence, h, each
# iteration. Even though m consists only of p_gb, we utilize a distinct symbol
# for the hidden Markov parameter(s) because future models will almost certainly
# be more sophisticated and involve more parameters. Nevertheless, the sampling
# approach remains the same.
#
# Step 2: Sample th|h (theta given h)
#
# Given a realized famine mask, we can calculate the likelihood, or probability,
# of the parameter vector th. At each step, we propose a new value, th_prop,
# by making a Gaussian draw. We utilize a conventional Metropolis-Hastings
# accept/reject step based on the ratio of the likelihoods. If the proposed
# parameter vector is accepted, we replace th with th_prop; otherwise, we keep
# the original value of th. In this preliminary model, we do not account for
# prior probabilities on the parameter values, which amounts to assigning
# uniform priors.
#
# Step 3: Sample m|h (m given h)
#
# Given a realized famine mask, we calculate the likelihood, or probability,
# of the parameter "vector" m. As in Step 2, we make a Gaussian draw to propose
# a new value of m, m_prop, and utilize an accept/reject step based on the
# ratio of likelihoods.
#
# We repeat the preceding steps 2000 times, and use only the final 1000 values
# for the HDI calculation (i.e., we adopt a "burn-in" period of 1000 samples).

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

# The probability of good to bad transitions is 0.2
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

# N_vect is the vector of sample sizes. The power calculation involves
# determing the number of times each experiment/simulation succeeds for
# each value of the sample size. See the method exp_wrapper in
# power_analysis_functions.R for what constitutes success. For the
# simulations that include age-at-death observations, we utilize N/10
# skeletons with an age-at-death estimate for every radiocarbon samples.
N_vect <- 100*2^(0:6)

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
plot_analysis(analysis_name, N_vect, exp_per_N)