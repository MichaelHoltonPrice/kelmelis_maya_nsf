# TODO: consider making the number of age classes an input rather than hard
#       coding it
calc_P <- function(a, kappa, is_famine) {
    #x <- 0:99
    x <- 0:79
    cum_haz <- a[1]/a[2]*(1 - exp(-a[2]*x)) + a[3]*x - a[4]/a[5]*(1 - exp(a[5]*x))
    if (is_famine) {
        cum_haz <- cum_haz * (1 + kappa)
    }
    # The
    l_of_x <- exp(-cum_haz)

    # Calculate the survival probability, P
    #P <- l_of_x[2:100] / l_of_x[1:99]
    P <- l_of_x[2:80] / l_of_x[1:79]

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
    #F_ <- c(rep(0,14), F_, rep(0,55))
    F_ <- c(rep(0,14), F_, rep(0,35))
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
  #return(c(z_0, z[1:99]*P))
  return(c(z_0, z[1:79]*P))
}

calc_pop_size <- function(P_g, F_g, P_b, F_b,
                          famine_mask, return_aad_prob=FALSE) {
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
                                     aad_vect=c()) {
  use_age <- length(aad_vect) != 0
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

  # For the mortality, some parameterizations with very small values of the
  # mortality hazard at high ages evaluate to NaN. Return Inf if this happens.
  # Do the same for F_g in case I am forgetting an edge case.
  P_g <- calc_P(a, kappa, FALSE)
  if (any(is.na(P_g))) {
    return(Inf)
  }
  F_g <- calc_F(f, P_g)
  # TODO: consider whether this should throw an error
  if (any(is.na(F_g))) {
    return(Inf)
  }
  P_b <- calc_P(a, kappa, TRUE)
  if (any(is.na(P_b))) {
    return(Inf)
  }
  F_b <- calc_F(f, P_b)
  # TODO: consider whether this should throw an error
  if (any(is.na(F_g))) {
    return(Inf)
  }

  pop_size <- calc_pop_size(P_g, F_g, P_b, F_b, h, return_aad_prob=use_age)
  if (use_age) {
    aad_prob <- pop_size$aad_prob
    pop_size <- pop_size$pop_size
  }
  v <- pop_size/ sum(pop_size)
  lik_vect <- meas_matrix %*% v
  log_lik <- sum(log(lik_vect))
  if (use_age) {
    for (k in 1:length(aad_vect)) {
      ind_age <- aad_vect[k]
      log_lik <- log_lik + log(aad_prob[ind_age])
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
                                aad_vect,
                                verbose=FALSE) {
  # First, make a draw for the latent states given the Markov transition
  # probabilities, m0 (in this case, m0 has only one entry, p_gb, that
  # determines all the other entries)
  num_times <- length(tau)
  h0 <- sample_famine_mask(m0, num_times)

  use_age <- length(aad_vect) != 0
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
  assess_run_times <- FALSE
  if (assess_run_times) {
    t_mask <- 0
    t_th_given_h <- 0
    t_m_given_h <- 0
  }

  for (n_mcmc in 1:num_mcmc_steps) {
      if (verbose) {
        print('**********')
        print(n_mcmc)
        print('(a) draw h|m')
      }
      h_before <- h
      # probability for the next sampling step. If so, reject the h
      if (assess_run_times) {
        t0 <- proc.time()
      }
      h <- sample_famine_mask(m, num_times)
      if (assess_run_times) {
        t1 <- proc.time()
        dt <- as.numeric(t1-t0)
        dt <- dt[3]
        t_mask <- t_mask + dt
      }
 

      # Though rare, it's possible that this new sample yields an untenable
      # probability for the next sampling step. If so, reject the h
      log_lik_curr <- calc_log_prob_th_given_h(th,
                                               h, # the famine mask
                                               meas_matrix,
                                               aad_vect)
 
      if (!is.finite(log_lik_curr)) {
          h <- h_before
          log_lik_curr <- calc_log_prob_th_given_h(th,
                                                   h, # the famine mask
                                                   meas_matrix,
                                                   aad_vect)
      }

      if (verbose) {
        print('(b) sample th|h')
      }
      th_prop <- th + rnorm(length(th))*th_scale

      if (assess_run_times) {
        t0 <- proc.time()
      }
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

  if (assess_run_times) {
      print('Run times of mask, th_given_h, and m_given_h')
      print(t_mask)
      print(t_th_given_h)
      print(t_m_given_h)
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
    aad_vect_known <- sample(length(aad_prob_known),
                         Nage,
                         replace=T,
                         prob=aad_prob_known)
  } else {
      aad_vect_known <- c()
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
                                  aad_vect_known,
                                  verbose=verbose)

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
  print(paste0('Starting problem ', prob$counter, ' of ', prob$num_prob))
  exp_obj <- run_experiment(prob$p_gb,
                            prob$a0,
                            prob$kappa0,
                            prob$f0,
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
          #print(counter)
          #print(result)
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