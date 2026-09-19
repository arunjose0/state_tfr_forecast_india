library(rstan)
library(dplyr)

cat("Extracting posterior estimates of the global floor (F)...\n\n")

models <- c(
  "Main Model" = "models/fit_main.rds",
  "Floor Low"  = "models/fit_floor_low.rds",
  "Floor High" = "models/fit_floor_high.rds"
)

for (m_name in names(models)) {
  file_path <- models[[m_name]]
  
  if (file.exists(file_path)) {
    fit <- readRDS(file_path)
    
    # Extract the posterior draws for the global floor 'F'
    F_draws <- rstan::extract(fit, "F")$F
    
    F_median <- median(F_draws)
    F_lower  <- quantile(F_draws, 0.025)
    F_upper  <- quantile(F_draws, 0.975)
    
    cat(sprintf("%-12s: Median F = %.2f (95%% CI: %.2f - %.2f)\n", 
                m_name, F_median, F_lower, F_upper))
  } else {
    cat(sprintf("%-12s: File not found (%s)\n", m_name, file_path))
  }
}
cat("\nExtraction complete.\n")
