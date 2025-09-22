# Summary function
gerar_resumo <- function(obj, vtheta) {
  vies <- c(obj$mediaMC) - vtheta
  dp <- apply(obj$nlminb, 2, sd)
  eqm <- sqrt(vies^2 + dp^2)
  
  data.frame(
    mean = round(apply(obj$nlminb, 2, mean), 4), 
    biasR = round(c(obj$ViesR), 4), #Relative bias
    sd = round(dp, 4), # Standard deviation
    MSE = round(eqm, 4), # mean square error
    kurtosis = round(kurtosis(obj$nlminb), 4),
    skewness = round(skewness(obj$nlminb), 4),
    round(obj$TC, 4) # Coverage probability
  )
}


parameter_boxplots <- function(object_list, true_values, parameter_labels) {
  
  sample_sizes <- names(object_list)  # e.g., "20", "30", ..., "500"
  plots <- list()
  n_params <- ncol(object_list[[1]]$nlminb)  # number of parameters estimated
  
  for (param_index in seq_len(n_params)) {
    
    # Extract the current parameter across all simulations
    param_data <- map_dfc(object_list, ~ .x$nlminb[, param_index])
    colnames(param_data) <- sample_sizes
    
    # Reshape to long format
    long_data <- pivot_longer(param_data, cols = everything(), names_to = "n", values_to = "Value")
    
    # Order x-axis by sample size
    long_data$n <- factor(long_data$n, levels = sample_sizes)
    
    # Create boxplot
    p <- ggplot(long_data, aes(x = n, y = Value)) +
      geom_boxplot() +
      geom_hline(yintercept = true_values[param_index], color = "red", lwd = 0.4) +
      labs(x = "Sample size (n)", y = parameter_labels[param_index]) +
      theme_bw()
    
    plots[[param_index]] <- p
  }
  
  names(plots) <- parameter_labels
  return(plots)
}