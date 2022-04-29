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

#build_A <- function(P, F_) {
#    N <- length(P) + 1
#    A <- matrix(0, N, N)
#    A[1,] <- F_
#    for (n in 1:(N-1)) {
#        A[n+1,n] <- P[n]
#    }
#    return(A)
#}

project_pop <- function(z, P, F_) {
  z_0 <- sum(z[15:45] * F_[15:45])
  return(c(z_0, z[1:99]*P))
}

calc_pop_size <- function(P_g, F_g, P_b, F_b,
                          famine_mask, return_age_vect=FALSE) {
    # For now, use something of a kleuge to get the starting age structure by
    # first getting the age structure after multiplying the entire sequence
    # of matrices out, starting with a uniform age structure (this requires two
    # rather than one set of multiplications but is a perfectly fine, albeit
    # slow way of getting the starting age structure; it could instead be
    # approximated by solving for the stable age structure of the Leslie matrix
    # or using Tulja's small noise approximation).
    N <- length(F_g)

    # To reduce package dependencies, use a for loop for the matrix
    # multiplications rather than taking the power of the matrices for the
    # famine and non-famine cases.
    z <- rep(1,N)
    for (is_famine in famine_mask) {
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

    if (return_age_vect) {
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
        if (return_age_vect) {
          Z[,counter+1] <- z
        }
    }

    # Use the mean population at the start and end of the period as the
    # population estimate.
    pop_size <- (pop_size[1:num_periods] + pop_size[2:(num_periods+1)])/2
    if (!return_age_vect) {
      return(pop_size)
    } else {
      age_vect <- rowSums(Z)
      age_vect <- age_vect / sum(age_vect)
      return(list(pop_size=pop_size, age_vect=age_vect))
    }
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
calc_log_prob_th_given_h <- function(th,
                                     h, # the famine mask
                                     meas_matrix,
                                     known_ages=c()) {
  use_age <- length(known_ages) != 0
  # Extract parameters
  a1    <- th[1]
  a2    <- th[2]
  a3    <- th[3]
  a4    <- th[4]
  a5    <- th[5]
  kappa <- th[6]
  f     <- th[7]

  # TODO: ensure parameter vector is valid
  if (!is_th_valid(th)) {
      return(-Inf)
  }
   
  a <- th[1:5]

  P_g <- calc_P(a, kappa, FALSE)
  F_g <- calc_F(f, P_g)
  P_b <- calc_P(a, kappa, TRUE)
  F_b <- calc_F(f, P_b)

  pop_size<- calc_pop_size(P_g, F_g, P_b, F_b, h, return_age_vect=use_age)
  if (use_age) {
    age_vect <- pop_size$age_vect
    pop_size <- pop_size$pop_size
  }
  v <- pop_size/ sum(pop_size)
  lik_vect <- meas_matrix %*% v
  log_lik <- sum(log(lik_vect))
  if (use_age) {
    for (n_a in 1:length(known_ages)) {
      ind_age <- known_ages[n_a]
      log_lik <- log_lik + log(age_vect[ind_age])
    }
  }
  return(log_lik)
}

calc_markov_probs <- function(p_gb) {
  H <- matrix(NA,2,2)
  H[1,1] <- 1 -   p_gb
  H[1,2] <-       p_gb
  H[2,1] <- 1 - 2*p_gb
  H[2,2] <-     2*p_gb

  # The overall probability is the dominant left eigenvector, normalized to
  # sum to 1
  w0 <- eigen(t(H))$vectors[,1]
  w0 <- w0 / sum(w0)
  return(w0)
}

calc_log_prob_m_given_h <- function(m, h) {
  p_gb <- m[1]
  if (p_gb <= 0) {
      return(Inf)
  }

  if (p_gb >= 0.5) {
      return(Inf)
  }

  H <- matrix(NA,2,2)
  H[1,1] <- 1 -   p_gb
  H[1,2] <-       p_gb
  H[2,1] <- 1 - 2*p_gb
  H[2,2] <-     2*p_gb

  # The overall probability is the dominant left eigenvector, normalized to
  # sum to 1
  w0 <- eigen(t(H))$vectors[,1]
  w0 <- w0 / sum(w0)

  if (h[1]) {
    # Initial period is famine
    log_lik <- log(w0[2])
  } else {
    # Initial period is not famine
    log_lik <- log(w0[1])
  }

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
                                known_ages,
                                verbose=FALSE) {
  # First, make a draw for the latent states given the Markov transition
  # probabilities, m0 (in this case, m0 has only one entry, p_gb, that
  # determines all the other entries)
  num_times <- length(tau)
  h0 <- sample_famine_mask(m0, num_times) # the number of time periods,
                                          # 100, need not be hard coded

  use_age <- length(known_ages) != 0
  burn_in <- 1000
  thinning <- 10
  total_samples <- 100

  num_mcmc_steps <- burn_in + thinning * total_samples
  th <- th0
  m <- m0
  h <- h0
  #eta <- eta0

  # TODO: add prior to calculation
  TH <- matrix(NA,length(th), num_mcmc_steps)
  H <- matrix(NA,length(h), num_mcmc_steps)
  M <- matrix(NA,length(m), num_mcmc_steps)
  for (n_mcmc in 1:num_mcmc_steps) {
      if (verbose) {
        print('**********')
        print(n_mcmc)
        print('(a) draw h|m')
      }
      h_before <- h
      h <- sample_famine_mask(m, num_times)
      # Though rare, it's possible that this new sample yields an untenable
      # probability for the next sampling step. If so, reject the h
      log_lik_curr <- calc_log_prob_th_given_h(th,
                                               h, # the famine mask
                                               meas_matrix,
                                               known_ages)
 
      if (!is.finite(log_lik_curr)) {
          h <- h_before
          log_lik_curr <- calc_log_prob_th_given_h(th,
                                                   h, # the famine mask
                                                   meas_matrix,
                                                   known_ages)
      }

      if (verbose) {
        print('(b) sample th|h')
      }
      th_prop <- th + rnorm(length(th))*th_scale
      log_lik_prop <- calc_log_prob_th_given_h(th_prop,
                                               h, # the famine mask
                                               meas_matrix,
                                               known_ages)

      if (!is.finite(log_lik_prop)) {
          accept <- FALSE
          if (verbose) {
            print('Not finite')
          }
      } else {
          # Calculate the acceptance parameter
          a <- min(1, exp(log_lik_prop - log_lik_curr))
          if (verbose) {
            print(a)
          }
          accept <- runif(1) < a
      }
      if (verbose) {
        print(accept)
      }

      if (accept) {
          th <- th_prop
      }

      if (verbose) {
        print('(c) sample m|h')
      }
      m_prop <- m + rnorm(length(m))*m_scale
      log_lik_curr <- calc_log_prob_m_given_h(m, h)
      if (!is.finite(log_lik_curr)) {
          stop('log_lik_curr is not finite. must handle this')
      }
      # Really prob_h_given_m
      log_lik_prop <- calc_log_prob_m_given_h(m_prop, h)

      if (!is.finite(log_lik_prop)) {
          if (verbose) {
            print('Not finite')
          }
          accept <- FALSE
      } else {
          # Calculate the acceptance parameter
          a <- min(1, exp(log_lik_prop - log_lik_curr))
          if (verbose) {
            print(a)
          }
          accept <- runif(1) < a
      }
      if (verbose) {
        print(accept)
      }

      if (accept) {
          m <- m_prop
      }

      TH[,n_mcmc] <- th
      H [,n_mcmc]  <- h
      M [,n_mcmc]  <- m
      if (verbose) {
        print(paste0('p_gb = ', m[1]))
      }
  }

  # At least for now, do not thin
  return(list(TH=TH,H=H,M=M))
}

# To fully take advantange of multiple cores, run each experiment inside a
# standalone function, run_experiment.
run_experiment <- function(p_gb,
                           a0,
                           kappa0,
                           f0,
                           N,
                           num_times,
                           rand_seed,
                           use_age,
                           verbose=FALSE) {
  set.seed(rand_seed)

  # Create simulated data
  P_g <- calc_P(a0, kappa0, FALSE)
  F_g <- calc_F(f0, P_g)
  P_b <- calc_P(a0, kappa0, TRUE)
  F_b <- calc_F(f0, P_b)

  # For now, set th0 equal to the true parameter vector
  # th has ordering [p_gb, a1, ..., a5, kappa, f]
  th0 <- c(a0, kappa0, f0)
  m0 <- p_gb
  # Use th0 / 20 as the scale for the proposal distribution in MCMC sampling
  # dividing by 5 seems about right for small sample sizes
  # dividing by 10 seems about right for N=100
  th_scale <- th0 / 10
  m_scale <- m0 / 10

  famine_mask_known <- sample_famine_mask(p_gb, num_times)
  pop_size_known <- calc_pop_size(P_g, F_g, P_b, F_b,
                                  famine_mask_known, return_age_vect=use_age)
  if (use_age) {
      age_vect <- pop_size_known$age_vect
      pop_size_known <- pop_size_known$pop_size
      Nage <- round(N*.2)
  }
  # Normalize the population size to sum to 1
  # dtau is 1, so pop_size_known  and M have the correct normalizations.
  # pop_size_known is the same as v in the JAS article.
  pop_size_known <- pop_size_known / sum(pop_size_known)
  true_calendar_dates <- sample(length(pop_size_known),
                                       N,
                                       replace=T,
                                       prob=pop_size_known)
  # Make the dates calendar dates, AD
  start_date <- 600
  true_calendar_dates <- true_calendar_dates + start_date - 1
  rc_meas <-
    draw_rc_meas_using_date(true_calendar_dates,
                            calib_df,
                            error_spec,
                            is_AD=TRUE)
  tau <- seq(start_date, start_date + num_times -1)
  meas_matrix <- calc_meas_matrix(tau, rc_meas$phi_m, rc_meas$sig_m, calib_df)
  if (use_age) {
    known_ages <- sample(length(age_vect),
                         Nage,
                         replace=T,
                         prob=age_vect)
  } else {
      known_ages <- c()
  }
  samp_obj <- sample_param_vector(th0,
                                  m0,
                                  rc_meas,
                                  calib_df,
                                  tau,
                                  num_famine_samps,
                                  meas_matrix,
                                  th_scale,
                                  m_scale,
                                  known_ages,
                                  verbose=verbose)

  sim_info <- list(th_known=th0,
                   m_known=m0,
                   h_known=famine_mask_known,
                   true_calendar_dates=true_calendar_dates,
                   pop_size_known=pop_size_known,
                   rc_meas=rc_meas,
                   tau=tau,
                   known_ages=known_ages)
  return(list(sim_info=sim_info, samp_obj=samp_obj))
}

exp_wrapper <- function(prob) {
  print(paste0('Starting problem ', prob$counter, ' of ', prob$num_prob))
  exp_obj <- run_experiment(prob$p_gb,
                            prob$a0,
                            prob$kappa0,
                            prob$f0,
                            prob$N,
                            prob$num_times,
                            prob$seed,
                            prob$use_age)
  f_lo <- 4
  f_hi <- 6
  hdi_obj <- hdi(exp_obj$samp_obj$TH[7,1000:2000])
  hdi_interval <- as.numeric(hdi_obj)
  #success <- (hdi_interval[1] <= f_lo) && (f_hi <= hdi_interval[2])
  success <- (f_lo <= hdi_interval[1]) && (hdi_interval[2] <= f_hi)
  file_path <- file.path('data',
                         paste0(prob$analysis_name, '_', prob$counter, '.rds'))
  save_obj <- list(exp_obj=exp_obj,
                   f_lo=f_lo,
                   f_hi=f_hi,
                   hdi_obj=hdi_obj,
                   hdi_interval=hdi_interval,
                   success=success)
  saveRDS(save_obj, file_path)
  return(success)
}

#plot_sample_n <- function(exp_obj, n) {
#  th <- exp_obj$samp_obj$TH[,1]
#  a <- th[1:5]
#  kappa <- th[6]
#  f     <- th[7]
#
#  start_date <- 600 # AD 600 is the first date
#  tau <- exp_obj$sim_info$tau
#
#  P_g <- calc_P(a, kappa, FALSE)
#  F_g <- calc_F(f, P_g)
#  P_b <- calc_P(a, kappa, TRUE)
#  F_b <- calc_F(f, P_b)
#
#  h <- exp_obj$samp_obj$H[,n]
#  v<- calc_pop_size(P_g, F_g, P_b, F_b, h)
#  v <- v / sum(v)
#
#  v_known <- exp_obj$sim_info$pop_size_known
#  ymax <- max(v, v_known)
#  pdf(paste0('plot_', n, '.pdf'))
#    plot(tau, v_known, col='green', type='l', lwd=3, ylim=c(0,ymax))
#    lines(tau, v, col='black', lwd=3)
#  dev.off()
#} 

plot_analysis <- function(analysis_name, N_vect, exp_per_N) {
  # Get all files in the data directory
  files <- list.files('data')

  # Subset to files beginning with the input analysis_name
  ind <- unlist(lapply(files, function(f){startsWith(f,analysis_name)}))
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
          success_array[k1, k2, k3] <- result$success
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
    plot(N_vect, power_vect_no_age, xlab='N', ylab='Power', ylim=c(0,y_max), col='black')
    points(N_vect, power_vect_age, col='red')
  dev.off()
}