calc_P <- function(a, kappa, is_famine) {
    # Calculate the age-specifc survival probabilities (the off-diagonal of the
    # population projection matrix)
    # a         -- The Siler parameter vector
    # kappa     -- The additional boost to mortality in famine years (likely 0.2)
    # is_famine -- Is this a famine year?

    # Create a vector of ages from 0 to 79 years
    x <- 0:79
    # Calculate the cumulative hazard
    cum_haz <- a[1]/a[2]*(1 - exp(-a[2]*x)) + a[3]*x - a[4]/a[5]*(1 - exp(a[5]*x))

    if (is_famine) {
        # During a famine year, both the hazard and cumulative hazard are
        # (1 + kappa) larger
        cum_haz <- cum_haz * (1 + kappa)
    }
    # Calculate survival from age 0 to age x from the cumulative hazard
    l_of_x <- exp(-cum_haz)

    # The survival probability from the beginning to end of a time period is
    # calculated from the ratio of the l_of_x values.
    P <- l_of_x[2:80] / l_of_x[1:79]

    return(P)
}

calc_F <- function(TFR, P) {
    # Calculate the age-specific fertilities (the top row of the population
    # projection matrix)
    # TFR -- The total fertility rate (TFR)
    # P -- The vector of age-specific survival probabilities
    #
    # The total fertility rate (TFR) is the mean number of children a woman
    # would have in the absence of mortality. This is not just the sum of
    # the age-specific fertilities, F, since these values account for
    # age-specific mortality. Equation 3.3.13 in Keyfitz and Caswell (2005)
    # gives the following approximate formula that relates instantaneous,
    # age-specific fertility (m_i) to F_i [we use the symbol m differently
    # below in the rest of the code base]:
    #
    # F_i = l(0.5) * (m_i + P_i * m_{i+1}) / 2,
    #
    # where l(0.5) is instant survival half way through the first age period. The
    # TFR is the sum of these m_i values. We make some further approximations to
    # keep the model simple. In particular, we assume that:
    #
    # (a) l(0.5) is close to 1
    # (b) m_i is constant across age classes
    # (c) (1 + P_i)/2 is approximately P_i
    #
    # Given these assumptions, the relationship simplifies to
    #
    # F_i = m P_i
    #
    # Aside from the preceding considerations, the population projection matrix
    # only tracks female births, so we must divide by 2. We also assume that the
    # reproductive span is 15 to 45 years.

    # Initialize F_ with the correct relative values (proportional to P_i)
    F_ <- P[15:45] # Here and elsewhere, use F_ rather than F since F is FALSE

    # Normalize F_ to yield the desired TFR
    F_ <- F_ / sum(F_)
    F_ <- F_ * TFR / 2 # divide by two to track just females

    # Build out the full vector, with 0 values included
    F_ <- c(rep(0,14), F_, rep(0,35))
    return(F_)
}

project_pop <- function(z, P, F_) {
  # Project the population one timestep into the future
  # z  -- The vector of population sizes by age
  # P  -- The age-specific survival probabilities
  # F_ -- The age-specific fertilities
  z_0 <- sum(z[15:45] * F_[15:45])
  return(c(z_0, z[1:79]*P))
}

calc_pop_size <- function(P_g, F_g, P_b, F_b,
                          famine_mask, return_aad_prob=FALSE) {
    # Calculate the population size given the population model for both
    # good and bad years (P_g, F_g, P_b, and F_b) and the mask of actual
    # famine years.  If return_aad_prob is TRUE, also calculate and
    # return the age at death values across time periods (aad stands
    # for age at death)
    #
    # P_g -- The age-specific survival probabilities for good/non-famine years
    # F_g -- The age-specific fertilities for good/non-famine years
    # P_b -- The age-specific survival probabilities for bad/famine years
    # F_b -- The age-specific fertilities for bad/famine years
    # famine_mask -- True if a year is bad/famine and False otherwise

    # For now, use something of a kleuge to get the starting age structure by
    # first getting the age structure after multiplying the entire sequence
    # of matrices out, starting with a uniform age structure (this requires two
    # rather than one set of multiplications but is a perfectly fine, albeit
    # slow way of getting the starting age structure; it could instead be
    # approximated by solving for the stable age structure of the Leslie matrix
    # or using Tulja's small noise approximation).
    # N_a The number of age groups
    N_a <- length(F_g)

    if (any(is.na(P_g))) {
      print(P_g)
      stop("P_g contains one or more NAs")
    }
    if (any(is.na(F_g))) {
      stop("F_g contains one or more NAs")
    }
    if (any(is.na(P_b))) {
      stop("P_b contains one or more NAs")
    }
    if (any(is.na(F_b))) {
      stop("F_b contains one or more NAs")
    }

    # To reduce package dependencies, use a for loop for the matrix
    # multiplications rather than taking the power of the matrices for the
    # famine and non-famine cases.
    z <- rep(1,N_a)
    counter <- 0
    for (is_famine in famine_mask) {
      counter <- counter + 1
      if (is_famine) {
        z <- project_pop(z, P_b, F_b)
      } else {
        z <- project_pop(z, P_g, F_g)
      }
    }
    z <- z / sum(z)

    num_periods <- length(famine_mask)
    pop_size <- rep(NA, num_periods+1)
    pop_size[1] <- 1
    counter <- 0

    if (return_aad_prob) {
      Z <- matrix(NA, length(z), num_periods+1)
      Z[,1] <- z
    }
    for (is_famine in famine_mask) {
        if (is_famine) {
            z <- project_pop(z, P_b, F_b)
        } else {
            z <- project_pop(z, P_g, F_g)
        }
        counter <- counter + 1
        pop_size[counter+1] <- sum(z)
        if (return_aad_prob) {
          Z[,counter+1] <- z
        }
    }

    # Use the mean population at the start and end of the period as the
    # population estimate.
    pop_size <- (pop_size[1:num_periods] + pop_size[2:(num_periods+1)])/2
    if (!return_aad_prob) {
      return(pop_size)
    } else {
      if (any(is.na(Z))) {
        stop('Z no entries of Z should be negative')
      }
      if (any(Z < 0 )) {
        stop('There should be no negative entries in Z')
      }
      # Create the age-at-death probability vector, aad_prob
      # TODO: the following for loop is likely why this function takes a good
      #       deal longer to execute when return_aad_prob is TRUE. It could be
      #       vectorized.
      aad_prob <- rep(0, N_a)
      for (n_a in 1:(N_a-1)) {
        for (n_t in 1:num_periods) {
          # For this age group, add the number of  deaths between time periods
          aad_prob[n_a] <- aad_prob[n_a] - Z[n_a+1,n_t+1] + Z[n_a,n_t]
        }
        # Everybody in the last age group dies
        aad_prob[N_a] <- aad_prob[N_a] + Z[N_a,n_t]
      }
      aad_prob <- aad_prob / sum(aad_prob)
      return(list(pop_size=pop_size, aad_prob=aad_prob))
    }
}

sample_famine_mask <- function(p_gb, num_periods) {
  # Create a sampled famine mask
  # p_gb        -- The probability of transitioning from a good to bad year
  # num_periods -- The numper of periods (years) to sample
  H <- matrix(NA,2,2)
  H[1,1] <- 1 -   p_gb
  H[1,2] <-       p_gb
  H[2,1] <- 1 - 2*p_gb
  H[2,2] <-     2*p_gb

  # First, sample from the overall probabilities of good and bad years to set
  # the mask for the first year. This overall probability is the dominant left
  # eigenvector of the matrix H, normalized to sum to 1.
  w0 <- eigen(t(H))$vectors[,1]
  w0 <- w0 / sum(w0)
  famine_mask <- sample(c(F,T), 1, replace=FALSE, prob=w0)

  # Now that we've initialized the mask for the first year, use the transition
  # probabilities to sample the mask for the remaining years.
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
    # Is the input the vector valid? Each element of th must be positive
    # th -- The parameter vector
    return(!any(th <= 0))
}

calc_log_prob_th_given_h <- function(th,
                                     h, # the famine mask
                                     meas_matrix,
                                     aad_vect=c()) {
    # Calculate the probability of th given a known famine mask, h
    # th          -- The parameter vector, [a1, ..., a5, kappa, TFR]
    # h           -- The famine mask (TRUE if a year is bad and FALSE otherwise)
    # meas_matrix -- The measurement matrix for the radiocarbon samples
    use_age <- length(aad_vect) != 0
    # Extract parameters
    a1    <- th[1]
    a2    <- th[2]
    a3    <- th[3]
    a4    <- th[4]
    a5    <- th[5]
    kappa <- th[6]
    TFR   <- th[7]
    a <- th[1:5]

    # Check that the parameter vector is valid. If not, return -Inf, which leads
    # to the corresponding sample being rejected.
    if (!is_th_valid(th)) {
        return(-Inf)
    }
  
    # For the mortality, some parameterizations with very small values of the
    # mortality hazard at high ages evaluate to NaN. Return Inf if this happens.
    # Do the same for F_g in case I am forgetting an edge case.
    P_g <- calc_P(a, kappa, FALSE) # age-specific survival probabilities good years
    if (any(is.na(P_g))) {
      return(Inf)
    }
    F_g <- calc_F(TFR, P_g) # age-specific fertilities good years
    if (any(is.na(F_g))) {
      return(Inf)
    }
    P_b <- calc_P(a, kappa, TRUE) # age-specific survival probabilities bad years
    if (any(is.na(P_b))) {
      return(Inf)
    }
    F_b <- calc_F(TFR, P_b) # age-specific fertilities bad years
    if (any(is.na(F_g))) {
      return(Inf)
    }
  
    # Calculate the population size (including the age-at-death vectors, if
    # necessary) given the demographic rates and mask of famine years.
    pop_size <- calc_pop_size(P_g, F_g, P_b, F_b, h, return_aad_prob=use_age)
    if (use_age) {
      aad_prob <- pop_size$aad_prob
      pop_size <- pop_size$pop_size
    }
    
    # Apply the likelihood function for the radiocarbon dates (v is the same as
    # in Price et al. 2021 -- End-to-end Bayesian analysis for summarizing sets
    # of radiocarbon dates).
    v <- pop_size/ sum(pop_size)
    lik_vect <- meas_matrix %*% v
    log_lik <- sum(log(lik_vect))

    # If necessary, use ages-at-death in the likelihood calculation
    if (use_age) {
      for (k in 1:length(aad_vect)) {
        ind_age <- aad_vect[k]
        log_lik <- log_lik + log(aad_prob[ind_age])
      }
    }
    return(log_lik)
}

calc_log_prob_m_given_h <- function(m, h) {
    # Calculate the log likelihood of the parameter m given actual/sampled hidden
    # states h. m has only one entry, p_gb = m[1], but for more complicated hidden
    # markov models would include more terms.
    # m -- The parameter vector for the hidden markov component of the model
    # h -- The vector mask of famine years (TRUE for a famine/bad year, FALSE otherwise)

    # p_gb cannot be negative
    p_gb <- m[1]
    if (p_gb <= 0) {
        return(Inf)
    }
  
    # p_gb cannot be greater than 0.5
    if (p_gb >= 0.5) {
        return(Inf)
    }
  
    # Create a matrix of transition probabilities, H (unrelated to little h),
    # so we can calculate the stable probabilities of good and bad years
    H <- matrix(NA,2,2)
    H[1,1] <- 1 -   p_gb
    H[1,2] <-       p_gb
    H[2,1] <- 1 - 2*p_gb
    H[2,2] <-     2*p_gb
  
    # The overall probability of good versus bad years is calculated from the dominant
    # left eigenvector of H (normalized to sum to 1)
    w0 <- eigen(t(H))$vectors[,1]
    w0 <- w0 / sum(w0)

    # Calculate the contribution of the first entry to the log likelihood

    if (h[1]) {
        # Initial period is famine
        log_lik <- log(w0[2])
    } else {
        # Initial period is not famine
        log_lik <- log(w0[1])
    }

    # Iterate over the remaning time periods to calculate the log likelihoods
    num_periods <- length(h)
    for (n in 2:num_periods) {
        if(h[n-1]) {
            if (h[n]) {
                # Famine then Famine
                log_lik <- log_lik + log(2*p_gb)
            } else {
                # Famine then not Famine
                log_lik <- log_lik + log(1 - 2*p_gb)
            }
        } else {
            if (h[n]) {
                # Not Famine then Famine
                log_lik <- log_lik + log(p_gb)
            } else {
                # Not Famine then Not Famine
                log_lik <- log_lik + log(1 - p_gb)
            }
        }
    }
    return(log_lik)
}

# Sample from the posterior distribution of the parameter vector
sample_param_vector <- function(th0,
                                m0,
                                rc_meas,
                                calib_df,
                                tau,
                                num_famine_samps,
                                meas_matrix,
                                th_scale,
                                m_scale,
                                aad_vect,
                                verbose=FALSE) {
    # Do the Bayesian sampling for the full model. This occurs in three,
    # repeating accept/reject steps (| stands for given):
    #
    # (a) Sample  h|m, where h is the famine mask and m = c(p_gb)
    # (b) Sample th|h, where th is the parameter vector for the demographic model
    # (c) Sample  m|h
    # (d) repeat (a) through (c)
    # First, make a draw for the latent states given the Markov transition
    # probabilities, m0 (in this case, m0 has only one entry, p_gb, that
    # determines all the other entries)
    # Create an initial sample of the famine mask using the known value of m0
    # (by ergodicity, how one initializes h0 does not change the sampling
    # stastitics if enough samples are used, and using m0 to initialize h0
    # reduces how long the "burn-in" needs to be)
    num_times <- length(tau) # number of time periods
    h0 <- sample_famine_mask(m0, num_times)

    # Do we have any age-at-death estimates?
    use_age <- length(aad_vect) != 0

    # Define our sampling parameters (these could be input control parameters)
    burn_in <- 1000
    thinning <- 10 # This is here for the future, but we are not actually using it yet
    total_samples <- 100

    num_mcmc_steps <- burn_in + thinning * total_samples

    # Initialize all the variables we need to sample
    th <- th0
    m <- m0
    h <- h0

    # TODO: Consider adding prior to calculation so we do not implicitly assume uniform priors
    # Define matrices to store our samples (one for each of th, m, and h)
    TH <- matrix(NA,length(th), num_mcmc_steps)
    H <- matrix(NA,length(h), num_mcmc_steps) # not the same as the transition matrix used elsewhere
    M <- matrix(NA,length(m), num_mcmc_steps)


    # A dev switch to assess the run times of different parts of the sampling
    assess_run_times <- FALSE
    if (assess_run_times) {
      t_mask <- 0
      t_th_given_h <- 0
      t_m_given_h <- 0
    }

    # Enter the core sampling loop
    for (n_mcmc in 1:num_mcmc_steps) {
        # (a) draw h|m
        if (verbose) {
            print('**********')
            print(n_mcmc)
            print('(a) draw h|m')
        }
        # Store h before we sample it
        h_before <- h
        if (assess_run_times) {
            t0 <- proc.time()
        }
        # Actually sample h
        h <- sample_famine_mask(m, num_times)
        if (assess_run_times) {
            t1 <- proc.time()
            dt <- as.numeric(t1-t0)
            dt <- dt[3]
            t_mask <- t_mask + dt
        }
 

        # Calculate the log-likelihood for this aspect of the sampling
        log_lik_curr <- calc_log_prob_th_given_h(th,
                                                 h, # the famine mask
                                                 meas_matrix,
                                                 aad_vect)
 
        # Though rare, it's possible that this new sample yields an untenable
        # probability for the next sampling step. If so, reject the h (which
        # means we need to recalculate the log_lik_curr value)
        if (!is.finite(log_lik_curr)) {
            h <- h_before
            log_lik_curr <- calc_log_prob_th_given_h(th,
                                                     h, # the famine mask
                                                     meas_matrix,
                                                     aad_vect)
        }

        # (a) draw th|h
        if (verbose) {
            print('(b) sample th|h')
        }
        # Use a scaled normal draw to propose a new value of th
        th_prop <- th + rnorm(length(th))*th_scale
  
        if (assess_run_times) {
            t0 <- proc.time()
        }
        # Calculate the log-likelihood for this aspect of the sampling
        log_lik_prop <- calc_log_prob_th_given_h(th_prop,
                                                 h, # the famine mask
                                                 meas_matrix,
                                                 aad_vect)
        if (assess_run_times) {
            t1 <- proc.time()
            dt <- as.numeric(t1-t0)
            dt <- dt[3]
            t_th_given_h <- t_th_given_h + dt
        }

        # Do the accept/reject logic
        if (!is.finite(log_lik_prop)) {
            # Always reject the new sample if the associated log likelihood is not finite 
            accept <- FALSE
            if (verbose) {
                print('Not finite')
            }
        } else {
            # The log likelihood is finite. Calculate the acceptance parameter
            # and randomly choose whether to accept the sample.
            a <- min(1, exp(log_lik_prop - log_lik_curr))
            if (verbose) {
                print(a)
            }
            accept <- runif(1) < a
        }
        if (verbose) {
            print(accept)
        }

        # Set th to equal th_prop if accept is TRUE
        if (accept) {
            th <- th_prop
        }

        # (c) draw m|h
        if (verbose) {
            print('(c) sample m|h')
        }
        # Use a scaled normal draw to propose a new value of m
        m_prop <- m + rnorm(length(m))*m_scale
        if (assess_run_times) {
            t0 <- proc.time()
        }
        log_lik_curr <- calc_log_prob_m_given_h(m, h)
        if (assess_run_times) {
            t1 <- proc.time()
            dt <- as.numeric(t1-t0)
            dt <- dt[3]
            t_m_given_h <- t_m_given_h + dt
        }
        if (!is.finite(log_lik_curr)) {
            stop('log_lik_curr is not finite. must handle this')
        }
        # Really prob_h_given_m in the sampling formula, but they are the same
        log_lik_prop <- calc_log_prob_m_given_h(m_prop, h)

        # Do the accept/reject logic
        if (!is.finite(log_lik_prop)) {
            # Always reject the new sample if the associated log likelihood is not finite 
            if (verbose) {
                print('Not finite')
            }
            accept <- FALSE
        } else {
            # The log likelihood is finite. Calculate the acceptance parameter
            # and randomly choose whether to accept the sample.
            a <- min(1, exp(log_lik_prop - log_lik_curr))
            if (verbose) {
                print(a)
            }
            accept <- runif(1) < a
        }
        if (verbose) {
            print(accept)
        }

        # Set th to equal th_prop if accept is TRUE
        if (accept) {
            m <- m_prop
        }

        # Store the final samples for this iteration
        TH[,n_mcmc] <- th
        H [,n_mcmc]  <- h
        M [,n_mcmc]  <- m
        if (verbose) {
            print(paste0('p_gb = ', m[1]))
        }
    }

    if (assess_run_times) {
        print('Run times of mask, th_given_h, and m_given_h')
        print(t_mask)
        print(t_th_given_h)
        print(t_m_given_h)
    }

    # At least for now, do not apply the thinning
    return(list(TH=TH,H=H,M=M))
}

run_experiment <- function(p_gb,
                           a0,
                           kappa0,
                           TFR0,
                           N,
                           num_times,
                           rand_seed,
                           use_age,
                           verbose=FALSE) {
    # Run one experiment. There are
    # two fundamental steps:
    #
    # (a) create the known/true simulated data
    # (b) run the Bayesian sampling algorithm.
    #
    # p_gb      -- The true probaiblity of transitioning from good to bad years  # TODO: this should be p_gb0
    # a0        -- The true Siler parameter vector
    # kappa0    -- The true mortality boost in famine/bad years
    # TFR0      -- The true total fertility rate
    # num_times -- The number of time periods to simulate
    # rand_seed -- A random number seed
    # use_age   -- Whether to use age-at-death observations in the likelihood calculation

    # Set the random number seed for reproducibility
    set.seed(rand_seed)

    # Create simulated data
    P_g <- calc_P(a0, kappa0, FALSE)
    F_g <- calc_F(TFR0, P_g)
    P_b <- calc_P(a0, kappa0, TRUE)
    F_b <- calc_F(TFR0, P_b)

    # Build the full, simulation parameter vectors
    # th has ordering [a1, ..., a5, kappa, TFR]
    th0 <- c(a0, kappa0, TFR0)
    m0 <- p_gb
    # Use th0 / 10 as the scale for the proposal distribution in MCMC sampling
    # dividing by 5 seems about right for small sample sizes
    # dividing by 10 seems about right for N=100. Using th0 to set the scale
    # should improve sampling for now (which is just a conveneince), but will not
    # be possible with an actual dataset.
    th_scale <- th0 / 10
    m_scale <- m0 / 10

    # Create a simulated, true (known) famine mask
    famine_mask_known <- sample_famine_mask(p_gb, num_times)
    # Create the simulated, true (known) population objects
    pop_size_known <- calc_pop_size(P_g, F_g, P_b, F_b,
                                  famine_mask_known, return_aad_prob=use_age)
    if (use_age) {
        aad_prob_known <- pop_size_known$aad_prob
        pop_size_known <- pop_size_known$pop_size
        Nage <- round(N*.1) # Assume 1 skeleton for 10 radiocarbon dates
    }
    # Normalize the population size to sum to 1
    # dtau is 1, so pop_size_known  and M have the correct normalizations.
    # pop_size_known is the same as v in the JAS article.
    pop_size_known <- pop_size_known / sum(pop_size_known)

    # Sample the simulated, true (known) calendar dates for the radiocarbon samples
    true_calendar_dates <- sample(length(pop_size_known),
                                         N,
                                         replace=T,
                                         prob=pop_size_known)
    # Make the dates calendar dates, AD
    start_date <- 600
    true_calendar_dates <- true_calendar_dates + start_date - 1
    # Sample the radiocarbon measurements, which accounts for measurement
    # uncertainty and the radiocarbon calibration curve.
    rc_meas <- draw_rc_meas_using_date(true_calendar_dates,
                                       calib_df,
                                       error_spec,
                                       is_AD=TRUE)
    # Create the time grid vector, tau
    tau <- seq(start_date, start_date + num_times -1)
    # Calculate the measurement matrix, M
    meas_matrix <- calc_meas_matrix(tau, rc_meas$phi_m, rc_meas$sig_m, calib_df)
    if (use_age) {
        # If necessary, sample the simulated, true (known) ages-at-death
        aad_vect_known <- sample(length(aad_prob_known),
                                 Nage,
                                 replace=T,
                                 prob=aad_prob_known)
    } else {
        aad_vect_known <- c()
    }
    
    # Do the actual sampling
    samp_obj <- sample_param_vector(th0,
                                    m0,
                                    rc_meas,
                                    calib_df,
                                    tau,
                                    num_famine_samps,
                                    meas_matrix,
                                    th_scale,
                                    m_scale,
                                    aad_vect_known,
                                    verbose=verbose)
  
    # Create the return object list and return it
    sim_info <- list(th_known=th0,
                     m_known=m0,
                     h_known=famine_mask_known,
                     true_calendar_dates=true_calendar_dates,
                     pop_size_known=pop_size_known,
                     rc_meas=rc_meas,
                     tau=tau,
                     aad_vect_known=aad_vect_known)
    return(list(sim_info=sim_info, samp_obj=samp_obj))
}

exp_wrapper <- function(prob) {
  # A wrapper to run the "experiment" / simulation and calculate success
  # for this experiment
  #
  # prob -- The problem, which is a list of values specifying this experiment
  print(paste0('Starting problem ', prob$counter, ' of ', prob$num_prob))
  # Run the actual experiment
  exp_obj <- run_experiment(prob$p_gb,
                            prob$a0,
                            prob$kappa0,
                            prob$TFR0,
                            prob$N,
                            prob$num_times,
                            prob$seed,
                            prob$use_age)
  # For the Bayesian power calculation, the region of practical equivalence
  # (ROPE) is 0 to 0.4 (the true value is 0.2)
  kappa_lo <- 0
  kappa_hi <- 0.4

  # Extract the vector of kappa values from the posterior sampling. The
  # first 1000 samples are "burn-in", so we discard them.
  kappa_vect <- exp_obj$samp_obj$TH[6,1000:2000]
  
  # Use the HDInterval package to determine the highest density interval (HDI)
  hdi_obj <- hdi(kappa_vect)
  hdi_interval <- as.numeric(hdi_obj)

  # The "experiment" is a success if the HDI lies entirely inside the ROPE
  success <- (kappa_lo <= hdi_interval[1]) && (hdi_interval[2] <= kappa_hi)

  # Save the data to file
  file_path <- file.path('data',
                         paste0(prob$analysis_name, '_', prob$counter, '.rds'))
  save_obj <- list(exp_obj=exp_obj,
                   kappa_lo=kappa_lo,
                   kappa_hi=kappa_hi,
                   hdi_obj=hdi_obj,
                   hdi_interval=hdi_interval,
                   success=success)
  saveRDS(save_obj, file_path)
  return(success)
}

plot_analysis <- function(analysis_name, N_vect, exp_per_N) {
  # Get all files in the data directory
  files <- list.files('data')

  # Subset to files beginning with the input analysis_name
  ind <- which(unlist(lapply(files, function(f){startsWith(f,analysis_name)})))
  files <- files[ind]

  # Extract the counter (i.e., experiment number within this analysis)
  counter_vect <-
    unlist(lapply(files,
                  function(f){as.numeric(strsplit(strsplit(f,'_')[[1]][2],
                                                           '\\.')[[1]][1])}))
  success_array <- array(NA,dim=c(length(N_vect), exp_per_N, 2))
  counter <- 0
  for (k1 in 1:length(N_vect)) {
    N <- N_vect[k1]
    for (k2 in 1:exp_per_N) {
      for (k3 in 1:2) {
        counter <- counter + 1
        if (counter %in% counter_vect) {
          file_path <- file.path('data',
                                 paste0(analysis_name,'_',counter,'.rds'))
          result <- readRDS(file_path)
          print('----')
          print(counter)
          success_array[k1, k2, k3] <- result$success
          if (k3 == 2) {
            print(result$success)
          }
        }
      }
    }
  }
  
  power_vect_no_age <- rep(NA, length(N_vect))
  power_vect_age    <- rep(NA, length(N_vect))
  for (k1 in 1:length(N_vect)) {
    outcomes_no_age <- success_array[k1,,1]
    outcomes_no_age <- outcomes_no_age[!is.na(outcomes_no_age)]
    outcomes_age    <- success_array[k1,,2]
    outcomes_age    <- outcomes_age[!is.na(outcomes_age)]
    power_vect_no_age[k1] <- sum(outcomes_no_age) / length(outcomes_no_age)
    power_vect_age   [k1] <- sum(outcomes_age) / length(outcomes_age)
  }
   
  y_max <- max(power_vect_no_age, power_vect_age, na.rm=TRUE)
  pdf(paste0('power_curve_',analysis_name,'.pdf'))
    plot(N_vect,
         power_vect_no_age,
         xlab='',
         ylab='Power',
         ylim=c(0,y_max),
         col='black',
         xaxt='n',
         type='o')
    points(N_vect, power_vect_age, col='red', type='o')
    axis(1, N_vect, N_vect)
    axis(1, N_vect, round(N_vect/10), line=2.5)
    x_offset <- 2400
    y_offset <- 1.00
    mtext("Num Rad Carb Samples", side=1, line=y_offset,       at=2400)
    mtext("Num Skeletal Samples", side=1, line=y_offset + 2.5, at=2400)
  dev.off()
}