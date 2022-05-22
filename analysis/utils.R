library(rcompanion)

test_assumptions <- function(model, transform, plot_hist = TRUE, hist_xlim = c(-6,6)) {
  transform <- paste("(", transform, ")")
  
  # Normality
  if (plot_hist){
    model_residuals = residuals(model)
    
    plotNormalHistogram(
      model_residuals, 
      main = paste("Linear Model Residuals", transform), 
      breaks = 100, 
      xlim = hist_xlim,
      xlab = "Residuals"
    )
  }
  
  # Q-Q Plot
  plot(model, which = c(2), main = paste("Residuals Q-Q Plot", transform))
  
  # Homocedasticity
  plot(model, which = c(1), main = paste("Variance of residuals", transform))
}