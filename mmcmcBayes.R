source("asgn_func.R")

mmcmcBayes <- function(cancer_data, normal_data,
                       stage = 1, max_stages = 3,
                       num_splits = 10, 
                       test = "BF", 
                       mcmc = NULL, 
                       priors = NULL, 
                       bf_thresholds = NULL, 
                       pvalue = NULL) {
  
  ####### Compute Mean Methylation #######
  calMean <- function(data) {
    data <- as.data.frame(data)
    required_cols <- c("CpG_ID", "Chromosome")
    
    if (!all(required_cols %in% colnames(data))) {
      cat("Warning: Missing metadata columns in methylation data.\n")
      return(NULL)
    }
    
    data_numeric <- data[, !(colnames(data) %in% required_cols)]
    if (ncol(data_numeric) == 0) return(NULL)
    
    mean_meth <- colMeans(data_numeric, na.rm = TRUE)
    return(matrix(mean_meth, ncol = 1))
  }
  
  ####### Compute ASGN Density #######
  asgn_density <- function(x, alpha, mu, sigma2) {
    numer <- sqrt(2) * ((1 - alpha * x)^2 + 1)
    denom <- 4 * gamma(3/2) * (alpha^2) + 4 * gamma(1/2)
    e <- exp(-((x - mu)^2) / (2 * sigma2))
    return((numer / denom) * e)
  }
  
  ####### Compute Bayes Factor #######
  calBF <- function(ybar_cancer, ybar_normal, posterior_cancer, posterior_normal) {
    likelihood_cancer <- sapply(ybar_cancer, asgn_density, 
                                alpha = posterior_cancer[1], 
                                mu = posterior_cancer[2], 
                                sigma2 = posterior_cancer[3])
    
    likelihood_normal <- sapply(ybar_normal, asgn_density, 
                                alpha = posterior_normal[1], 
                                mu = posterior_normal[2], 
                                sigma2 = posterior_normal[3])
    
    BF <- sum(likelihood_cancer) / sum(likelihood_normal)
    return(BF)
  }
  
  ####### Compute Anderson-Darling Test #######
  calAD <- function(ybar_cancer, ybar_normal, posterior_cancer, posterior_normal) {
    if (!requireNamespace("kSamples", quietly = TRUE)) {
      stop("The 'kSamples' package is required but not installed. Install it using install.packages('kSamples')")
    }
    
    density_cancer <- sapply(ybar_cancer, asgn_density, 
                             alpha = posterior_cancer[1], 
                             mu = posterior_cancer[2], 
                             sigma2 = posterior_cancer[3])
    
    density_normal <- sapply(ybar_normal, asgn_density, 
                             alpha = posterior_normal[1], 
                             mu = posterior_normal[2], 
                             sigma2 = posterior_normal[3])
    
    ad_test_result <- kSamples::ad.test(density_cancer, density_normal)
    p_value <- ad_test_result$ad[2, 3]
    
    return(p_value)
  }
  
  ####### Process the Current Stage #######
  if (is.null(mcmc)) {
    # Default MCMC settings
    mcmc <- list(nburn = 5000, niter = 10000, thin = 1)  
  }
  if (is.null(bf_thresholds)) {
    # Default BF thresholds
    bf_thresholds <- list(stage1 = 10, stage2 = 15, stage3 = 20)  
  }
  if (is.null(pvalue)) {
    # Default p-value thresholds
    pvalue <- list(stage1 = 10^(-4), stage2 = 10^(-6), stage3 = 10^(-8))  
  }
  
  total_cpgs <- nrow(cancer_data)
  if (total_cpgs == 0) {
    cat("Empty segment at Stage", stage, ". Skipping.\n")
    return(NULL)
  }
  
  ybar_cancer <- calMean(cancer_data)
  ybar_normal <- calMean(normal_data)
  
  if (is.null(ybar_cancer) || is.null(ybar_normal)) {
    cat("Invalid methylation data at Stage", stage, ". Skipping.\n")
    return(NULL)
  }
  
  posterior_cancer <- asgn_func(ybar_cancer, priors, mcmc)
  posterior_normal <- asgn_func(ybar_normal, priors, mcmc)
  
  ####### Select Statistical Test #######
  if (test == "BF") {
    BF <- calBF(ybar_cancer, ybar_normal, posterior_cancer, posterior_normal)
    threshold <- bf_thresholds[[paste0("stage", stage)]]
    decision_value <- BF
    decision_criteria <- (BF >= threshold)
    cat("Bayes Factor at Stage", stage, ":", BF, "\n")
  } else if (test == "AD") {
    p_value <- calAD(ybar_cancer, ybar_normal, posterior_cancer, posterior_normal)
    threshold <- pvalue[[paste0("stage", stage)]]
    decision_value <- p_value
    decision_criteria <- (p_value < threshold)
    cat("AD Test p-value at Stage", stage, ":", p_value, "\n")
  } else {
    stop("Invalid test. Choose either 'BF' or 'AD'.")
  }
  
  ####### Stopping Condition #######
  if (stage == max_stages || !decision_criteria) {
    if (decision_criteria) {
      detected_DMR <- data.frame(
        Chromosome = cancer_data$Chromosome[1],
        Start_CpG = cancer_data$CpG_ID[1],
        End_CpG = cancer_data$CpG_ID[nrow(cancer_data)],
        CpG_Count = total_cpgs,
        Decision_Value = decision_value
      )
      cat("Significant region detected at Stage", stage, "\n")
      return(detected_DMR)
    } else {
      cat("No significant difference at Stage", stage, "\n")
      return(NULL)
    }
  }
  
  ####### Split Data Into Sub-Segments #######
  cat("Stage", stage, ": Splitting into", num_splits, "sub-segments...\n")
  segment_size <- ceiling(total_cpgs / num_splits)
  split_indices <- split(seq_len(nrow(cancer_data)), cut(seq_len(nrow(cancer_data)), 
                                                         num_splits, labels = FALSE))
  
  cancer_sub_segments <- lapply(split_indices, function(idx) cancer_data[idx, , drop = FALSE])
  normal_sub_segments <- lapply(split_indices, function(idx) normal_data[idx, , drop = FALSE])
  
  ####### Process Each Sub-Segment #######
  results <- lapply(seq_along(cancer_sub_segments), function(i) {
    mmcmcBayes(cancer_sub_segments[[i]], normal_sub_segments[[i]], stage + 1, max_stages, num_splits, 
               test, mcmc, priors, bf_thresholds, pvalue)
  })
  
  results <- Filter(Negate(is.null), results)
  
  return(do.call(rbind, results))
}
