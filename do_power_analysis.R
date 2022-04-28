# 0. Clear the workspace and load packages
rm(list=ls())
library(baydem)
library(matrixStats) # for logSumExp

# 1. Set the true values

# The true total fertiltiy rate is 5 (half of babies are female)
f0 <- 5

# The true Siler parameter vector is from Gage and Dyke (1986), table TBD
a0 <- c(0.175,1.40,0.368*0.01,0.075*0.001,0.917*0.1)

# The true boost to mortality is kappa = 1.2 (a 20% increase)
kappa0 <- 1.2

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

calc_P <- function(a, kappa, is_famine) {
    x <- 0:99
    cum_haz <- a[1]/a[2]*(1 - exp(-a[2]*x)) + a[3]*x - a[4]/a[5]*(1 - exp(a[5]*x))
    if (is_famine) {
        cum_haz <- cum_haz * (1 + kappa)
    }
    # The
    l_of_x <- exp(-cum_haz)

    # Calculate the survival probability, P
    P <- l_of_x[2:100] / l_of_x[1:99]

    return(P)
}

calc_F <- function(f, P) {
    # Use 15 as the age of first reproduction (F1 to F14 are zero)
    # Use 45 as the age at last reproduction and normalize so that the number
    # of offspring sum to f
    # In addition, account for mortality over each time period in weighting the
    # fertility
    F_ <- P[15:45] # Here and elsewhere, use F_ rather than F since F is FALSE
    F_ <- F_ / sum(F_)
    F_ <- F_ * f / 2 # divide by two to track females
    F_ <- c(rep(0,14), F_, rep(0,55))
    return(F_)
}

build_A <- function(P, F_) {
    N <- length(P) + 1
    A <- matrix(0, N, N)
    A[1,] <- F_
    for (n in 1:(N-1)) {
        A[n+1,n] <- P[n]
    }
    return(A)
}

calc_pop_size <- function(A_g, A_b, famine_mask) {
    # For now, use something of a kleuge to get the starting age structure by
    # first getting the age structure after multiplying the entire sequence
    # of matrices out, starting with a uniform age structure (this requires two
    # rather than one set of multiplications but is a perfectly fine, albeit
    # slow way of getting the starting age structure; it could instead be
    # approximated by solving for the stable age structure of A_g or using
    # Tulja's small noise approximation).
    N <- nrow(A_g)

    # To reduce package dependencies, use a for loop for the matrix
    # multiplications rather than taking the power of the matrices for the
    # famine and non-famine cases.
    M <- diag(N)
    for (is_famine in famine_mask) {
        if (is_famine) {
            M <- A_b %*% M
        } else {
            M <- A_g %*% M
        }
    }
    z <- M %*% matrix(1,N,1)
    z <- z / sum(z)

    num_periods <- length(famine_mask)
    pop_size <- rep(NA, num_periods+1)
    pop_size[1] <- 1
    counter <- 0

    for (is_famine in famine_mask) {
        if (is_famine) {
            z <- A_b %*% z
        } else {
            z <- A_g %*% z
        }
        counter <- counter + 1
        pop_size[counter+1] <- sum(z)
    }

    # Use the mean population at the start and end of the period as the
    # population estimate.
    pop_size <- (pop_size[1:num_periods] + pop_size[2:(num_periods+1)])/2
    return(pop_size)
}

sample_famine_mask <- function(p_gb, num_periods) {
  H <- matrix(NA,2,2)
  H[1,1] <- 1 -   p_gb
  H[1,2] <-       p_gb
  H[2,1] <- 1 - 2*p_gb
  H[2,2] <-     2*p_gb

  # The overall probability is the dominant left eigenvector, normalized to
  # sum to 1
  w0 <- eigen(t(H))$vectors[,1]
  w0 <- w0 / sum(w0)
  famine_mask <- sample(c(F,T), 1, replace=FALSE, prob=w0)

  for (n in 2:num_periods) {
      if(famine_mask[n-1]) {
          # Last year was a famine
          w <- H[2,]
      } else {
          # Last year was not a famine
          w <- H[1,]
      }
      famine_mask <- c(famine_mask, sample(c(F,T), 1, replace=FALSE, prob=w))
  }
  return(famine_mask)
}

is_th_valid <- function(th) {
    return(!any(th <= 0))
}

# Calculate the negative log-likelihood
# th has ordering [p_gb, a1, ..., a5, kappa, f]
calc_neg_log_lik <- function(th,
                             rc_meas,
                             calib_df,
                             tau,
                             num_famine_samps,
                             meas_matrix) {
   print('top of calc_neg_log_lik')
   # Extract parameters
   p_gb  <- th[1]
   a1    <- th[2]
   a2    <- th[3]
   a3    <- th[4]
   a4    <- th[5]
   a5    <- th[6]
   kappa <- th[7]
   f     <- th[8]

   # TODO: ensure parameter vector is valid
   if (!is_th_valid(th)) {
       return(-Inf)
   }
   
   a <- th[2:6]

   P_g <- calc_P(a, kappa, FALSE)
   F_g <- calc_F(f, P_g)
   A_g <- build_A(P_g, F_g)
   P_b <- calc_P(a, kappa, TRUE)
   F_b <- calc_F(f, P_b)
   A_b <- build_A(P_b, F_b)

   num_years <- length(tau)
   famine_mask_matrix <- matrix(NA, num_famine_samps, num_years)
   log_lik_for_lse <- rep(NA, num_famine_samps)
   print('before famine sample for loop')
   t_mask <- 0
   t_v <- 0
   t_mult <- 0
   for (n_f in 1:num_famine_samps) {
     t0 <- proc.time()
     famine_mask_matrix[n_f,] <- sample_famine_mask(p_gb, num_years)
     t1 <- proc.time()
     dt <- as.numeric(t1-t0)
     dt <- dt[3]
     t_mask <- t_mask + dt
     t0 <- proc.time()
     pop_size_known <- calc_pop_size(A_g, A_b, famine_mask_matrix[n_f,])
     v <- pop_size_known / sum(pop_size_known)
     t1 <- proc.time()
     dt <- as.numeric(t1-t0)
     dt <- dt[3]
     t_v <- t_v + dt

     t0 <- proc.time()
     log_lik_for_lse[n_f] <- log(1/num_famine_samps) + sum(log(meas_matrix %*% v))
     t1 <- proc.time()
     dt <- as.numeric(t1-t0)
     dt <- dt[3]
     t_mult <- t_mult + dt
   }
   print('after famine sample for loop')
   print(t_mask)
   print(t_v)
   print(t_mult)
   log_lik <- logSumExp(log_lik_for_lse)
   return(-log_lik)
}

# Sample from the posterior distribution of the parameter vector
sample_param_vector <- function(th0,
                                rc_meas,
                                calib_df,
                                tau,
                                num_famine_samps,
                                meas_matrix,
                                th_scale) {
  # Calculate the initial negative log-likelihood, and ensure it is finite
  eta0 <- calc_neg_log_lik(th0,
                           rc_meas,
                           calib_df,
                           tau,
                           num_famine_samps,
                           meas_matrix)

  if (!is.finite(eta0)) {
      stop('eta0 is not finite')
  }

  num_mcmc_samp <- 1000
  th <- th0
  eta <- eta0

  # TODO: add prior to calculation
  TH <- matrix(NA,length(th), num_mcmc_samp)
  for (n_mcmc in 1:num_mcmc_samp) {
      print('----')
      print(n_mcmc)
      print(eta)
      th_prop <- th + rnorm(length(th)) * th_scale
      eta_prop <- calc_neg_log_lik(th_prop,
                                   rc_meas,
                                   calib_df,
                                   tau,
                                   num_famine_samps,
                                   meas_matrix)
      if (!is.finite(eta_prop)) {
          accept <- FALSE
          print('Not finite')
      } else {
          # Calculate the acceptance parameter
          a <- min(1, exp(-(eta_prop - eta)))
          print(a)
          accept <- runif(1) < a
      }
      print(accept)

      if (accept) {
          th <- th_prop
          eta <- eta_prop
      }
      TH[,num_mcmc_samp] <- th
  }

  # At least for now, do not thin
  return(TH)
}

# [1 -   p_gb,   p_gb]
# [1 - 2*p_gb, 2*p_gb]
P_g <- calc_P(a0, kappa0, FALSE)
F_g <- calc_F(f0, P_g)
A_g <- build_A(P_g, F_g)
P_b <- calc_P(a0, kappa0, TRUE)
F_b <- calc_F(f0, P_b)
A_b <- build_A(P_b, F_b)

#famine_mask <- c(F, T, F, F, F, T, T, F, F, F)
p_gb <- .1
famine_mask <- sample_famine_mask(p_gb, 100)
pop_size <- calc_pop_size(A_g, A_b, famine_mask)

# First, define a function to calculate the population size given the input
# parameter vectors, including the realized famine years. To do so, utilize a
# discrete, age-structured population with one year age intervals and one year
# time steps.

calib_df <- load_calib_curve('intcal20')
error_spec=list(type="unif_fm",min=.0021,max=.0028)

start_date <- 600 # AD 600 is the first date
tau <- seq(start_date, start_date+ 100-1)
# The number of experiments for each value of the effect size

# To fully take advantange of multiple cores, run each experiment inside a
# standalone function, run_experiment.
run_experiment <- function(p_gb, a0, kappa0, f0, N, rand_seed) {
  set.seed(rand_seed)
  num_famine_samps <- 100

  # Create simulated data
  P_g <- calc_P(a0, kappa0, FALSE)
  F_g <- calc_F(f0, P_g)
  A_g <- build_A(P_g, F_g)
  P_b <- calc_P(a0, kappa0, TRUE)
  F_b <- calc_F(f0, P_b)
  A_b <- build_A(P_b, F_b)

  # For now, set th0 equal to the true parameter vector
  # th has ordering [p_gb, a1, ..., a5, kappa, f]
  th0 <- c(p_gb, a0, kappa0, f0)
  # Use th0 / 20 as the scale for the proposal distribution in MCMC sampling
  th_scale <- th0 / 5

  famine_mask_known <- sample_famine_mask(p_gb, 100)
  pop_size_known <- calc_pop_size(A_g, A_b, famine_mask)
  # Normalize the population size to sum to 1
  # dtau is 1, so pop_size_known  and M have the correct normalizations.
  # pop_size_known is the same as v in the JAS article.
  pop_size_known <- pop_size_known / sum(pop_size_known)
  true_calendar_dates <- sample(length(pop_size_known),
                                       N,
                                       replace=T,
                                       prob=pop_size_known)
  # Make the dates calendar dates, AD
  true_calendar_dates <- true_calendar_dates + start_date - 1
  rc_meas <-
    draw_rc_meas_using_date(true_calendar_dates,
                            calib_df,
                            error_spec,
                            is_AD=TRUE)
  meas_matrix <- calc_meas_matrix(tau, rc_meas$phi_m, rc_meas$sig_m, calib_df)
  TH <- sample_param_vector(th0,
                            rc_meas,
                            calib_df,
                            tau,
                            num_famine_samps,
                            meas_matrix,
                            th_scale)

  return(TH)
}

p_gb <- .025
N <- 100

run_experiment(p_gb, a0, kappa0, f0, N, 100)
