traceplot_asgn <- function(result) {
  if (!is.list(result) || is.null(result$samples)) {
    stop("Invalid input. Provide the result from `asgn_func()`.")
  }
  
  par(mfrow = c(3, 1))  # Arrange plots for alpha, mu, sigma2
  
  param_names <- c("alpha", "mu", "sigma2")
  
  for (param in param_names) {
    trace_vals <- result$samples[[param]]
    plot(trace_vals, type = "l",
         main = paste("Trace Plot:", param),
         xlab = "Iteration", ylab = param)
    abline(h = mean(trace_vals), col = "red", lwd = 2, lty = 2)  # Posterior mean
  }
  
  par(mfrow = c(1, 1))  # Reset plotting layout
}
