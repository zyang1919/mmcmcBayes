asgn_func <- function(data, 
                      priors = NULL, 
                      mcmc = list(nburn = 5000, niter = 10000, thin = 1), 
                      seed = 2021) {
  
  # Ensure MCMCpack is installed
  if (!requireNamespace("MCMCpack", quietly = TRUE)) {
    stop("MCMCpack package is required but not installed. Install it using install.packages('MCMCpack')")
  }
  
  ybar <- data
  
  # If this is Stage 1 (no priors provided), use empirical estimates
  if (is.null(priors)) {
    prior_mean_mu <- mean(ybar, na.rm = TRUE)  
    prior_sd_mu <- sd(ybar, na.rm = TRUE)  
    
    prior_mean_sigma2 <- var(ybar, na.rm = TRUE)  
    prior_sd_sigma2 <- prior_mean_sigma2 / 2  # Weakly informative
    
    priors <- list(
      alpha = 1,  # Default prior for alpha
      mu = prior_mean_mu,
      sigma2 = prior_mean_sigma2
    )
  }
  
  # Set default MCMC values if not provided
  nburn <- ifelse(is.null(mcmc$nburn), 5000, mcmc$nburn)
  niter <- ifelse(is.null(mcmc$niter), 10000, mcmc$niter)
  thin  <- ifelse(is.null(mcmc$thin), 1, mcmc$thin)
  
  # Extract priors for this stage
  alpha_pri <- as.numeric(priors$alpha);siga <- 1
  mu_pri <- as.numeric(priors$mu);sigm <- (prior_sd_mu)^2
  sigma2_pri <- as.numeric(priors$sigma2)
  
  # Regularized inverse gamma prior
  As <- (sigma2_pri^2) / (sigma2_pri^2 / 2) + 2  
  Bs <- sigma2_pri * ((sigma2_pri^2) / (sigma2_pri^2 / 2) + 1)
  
  # Initialize MCMC values
  alpha_t <- 1
  mu_t <- mu_pri
  sig_t <- sigma2_pri
  
  ## MCMC Info
  total_iter <- nburn + niter * thin  # Total iterations including thinning

  n <- nrow(ybar)
  
  ## Containers for posterior samples
  alpha_samples <- numeric(niter)
  mu_samples <- numeric(niter)
  sig_samples <- numeric(niter)
  
  ## MCMC Sampling
  save_index <- 1  # Tracks where to save samples after burn-in
  for (k in 1:total_iter) {
    
    ## --- Update Alpha (α) using Metropolis-Hastings ---
    alpha_c <- rnorm(1, alpha_t, sd=1)
    
    log_prior_alpha_c <- - (alpha_c - alpha_pri)^2 / (2 * siga)
    log_prior_alpha_t <- - (alpha_t - alpha_pri)^2 / (2 * siga)
    
    ## Include the missing term from the full conditional
    log_gamma_term_c <- -n * log(4 * gamma(3/2) * (alpha_c^2) + 4 * gamma(1/2))
    log_gamma_term_t <- -n * log(4 * gamma(3/2) * (alpha_t^2) + 4 * gamma(1/2))
    
    log_likelihood_c <- sum(log(((1 - alpha_c * ybar)^2 + 1)) - ((ybar - mu_t)^2) / (2 * sig_t))
    log_likelihood_t <- sum(log(((1 - alpha_t * ybar)^2 + 1)) - ((ybar - mu_t)^2) / (2 * sig_t))
    
    ## Compute the Metropolis-Hastings acceptance ratio
    log_alpha_ratio <- (log_prior_alpha_c + log_likelihood_c + log_gamma_term_c) - 
      (log_prior_alpha_t + log_likelihood_t + log_gamma_term_t)
    
    if (log(runif(1,0,1)) < log_alpha_ratio) {
      alpha_t <- alpha_c
    }
    
    ## --- Update Mu (μ) using Metropolis-Hastings ---
    mu_c <- rnorm(1, mu_t, sd=1)
    
    log_prior_mu_c <- - (mu_c - mu_pri)^2 / (2 * sigm)
    log_prior_mu_t <- - (mu_t - mu_pri)^2 / (2 * sigm)
    
    log_likelihood_mu_c <- sum(log(((1 - alpha_t * ybar)^2 + 1)) - ((ybar - mu_c)^2) / (2 * sig_t))
    log_likelihood_mu_t <- sum(log(((1 - alpha_t * ybar)^2 + 1)) - ((ybar - mu_t)^2) / (2 * sig_t))
    
    log_mu_ratio <- (log_prior_mu_c + log_likelihood_mu_c) - (log_prior_mu_t + log_likelihood_mu_t)
    
    if (log(runif(1,0,1)) < log_mu_ratio) {
      mu_t <- mu_c
    }
    
    ## --- Update Sigma² (σ²) using Inverse Gamma Sampling with Stability Check ---
    sig_c <- MCMCpack::rinvgamma(1, shape = (sig_t^2 + 2), scale = ((sig_t^2 + 1) * sig_t))
    
    ## Prevent numerical instability (avoid very small values leading to large variances)
    if (sig_c < 1e-300) sig_c <- 1e-300  
    
    log_prior_sig_c <- -As * log(sig_c) - (Bs / sig_c)
    log_prior_sig_t <- -As * log(sig_t) - (Bs / sig_t)
    
    log_likelihood_sig_c <- sum(log(((1 - alpha_t * ybar)^2 + 1)) - ((ybar - mu_t)^2) / (2 * sig_c))
    log_likelihood_sig_t <- sum(log(((1 - alpha_t * ybar)^2 + 1)) - ((ybar - mu_t)^2) / (2 * sig_t))
    
    log_sigma_ratio <- (log_prior_sig_c + log_likelihood_sig_c) - (log_prior_sig_t + log_likelihood_sig_t)
    
    if (log(runif(1,0,1)) < log_sigma_ratio) {
      sig_t <- sig_c
    }
    
    ## Store posterior samples after burn-in with thinning
    if (k > nburn && (k - nburn) %% thin == 0) {
      alpha_samples[save_index] <- alpha_t
      mu_samples[save_index] <- mu_t
      sig_samples[save_index] <- sig_t
      save_index <- save_index + 1
    }
  }
  
  ## Compute posterior means for next stage priors
  alpha_post <- mean(alpha_samples)
  mu_post <- mean(mu_samples)
  sig_post <- mean(sig_samples)
  
  ## Return results
  return(c(alpha_post, mu_post, sig_post))
}
