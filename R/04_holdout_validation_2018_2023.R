# ==============================================================================
# 04_holdout_validation_2018_2023.R
#
# Out-of-sample validation: Fits (or loads) a model trained strictly on data
# <= 2017. Evaluates forecast accuracy against observed SRS/NFHS data from
# 2018-2023 using Mean Absolute Error (MAE), Root Mean Squared Error (RMSE),
# and Continuous Ranked Probability Score (CRPS), reported for both the
# training-period fit and the out-of-sample test period. Generates faceted
# calibration plots.
# ==============================================================================

library(tidyverse)
library(here)
library(rstan)
library(scoringRules)

options(warn = 1, mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

cat("--- HOLDOUT VALIDATION (2018-2023) ---\n")

# ── 1. Load Data & Filter to Pre-2018 (<= 2017) ──────────────────────────────
prep <- readRDS(here("results", "prep_objects.rds"))
df_srs <- prep$df_srs
df_nfhs <- prep$df_nfhs
transition_map <- prep$transition_map

# Exclude Ladakh since it has no training data before 2018
df_srs_all  <- df_srs |> filter(state_clean != "Ladakh")
df_nfhs_all <- df_nfhs |> filter(state_clean != "Ladakh")

# Rebuild state mappings and groups to keep indices sequential
all_states <- sort(unique(c(df_srs_all$state_clean, df_nfhs_all$state_clean)))
S <- length(all_states)
state_map <- tibble(state = all_states, state_id = seq_len(S))

df_srs_all  <- df_srs_all  |> select(-state_id) |> left_join(state_map, by = c("state_clean" = "state"))
df_nfhs_all <- df_nfhs_all |> select(-state_id) |> left_join(state_map, by = c("state_clean" = "state"))

state_map_with_groups <- state_map |> left_join(transition_map, by = c("state" = "state_clean"))
group_of_state <- state_map_with_groups$group_id

df_srs_train  <- df_srs_all  |> filter(year <= 2017)
df_nfhs_train <- df_nfhs_all |> filter(year <= 2017)
df_srs_test <- df_srs_all |> filter(year > 2017)
df_nfhs_test <- df_nfhs_all |> filter(year > 2017)

# ── 2. Rebuild state/time grids for training ──────────────────────────────────
state_year_grids <- list()
T_s <- numeric(S)
tfr_init <- numeric(S)
tfr_last <- numeric(S)
last_obs_year <- numeric(S)

for (s in 1:S) {
  yrs_srs  <- df_srs_train  |> filter(state_id == s) |> pull(year)
  yrs_nfhs <- df_nfhs_train |> filter(state_id == s) |> pull(year)
  comb_yrs <- sort(unique(c(yrs_srs, yrs_nfhs)))
  
  state_year_grids[[s]] <- comb_yrs
  T_s[s] <- length(comb_yrs)
  
  sub_srs  <- df_srs_train  |> filter(state_id == s) |> arrange(year)
  sub_nfhs <- df_nfhs_train |> filter(state_id == s) |> arrange(year)
  tfr_init[s] <- if (nrow(sub_srs) > 0) sub_srs$tfr_total[1] else sub_nfhs$tfr_total[1]
  
  last_yr <- max(comb_yrs)
  last_obs_year[s] <- last_yr
  last_obs_srs <- df_srs_train |> filter(state_id == s, year == last_yr)
  if (nrow(last_obs_srs) > 0) {
    tfr_last[s] <- last_obs_srs$tfr_total[1]
  } else {
    last_obs_nfhs <- df_nfhs_train |> filter(state_id == s, year == last_yr)
    tfr_last[s] <- last_obs_nfhs$tfr_total[1]
  }
}

max_T <- max(T_s)
dt <- matrix(1.0, nrow = S, ncol = max_T)
for (s in 1:S) {
  comb_yrs <- state_year_grids[[s]]
  if (length(comb_yrs) > 1) {
    for (t in 2:length(comb_yrs)) dt[s, t] <- comb_yrs[t] - comb_yrs[t - 1]
  }
}

df_srs_train <- df_srs_train |> rowwise() |> mutate(time_id = which(state_year_grids[[state_id]] == year)) |> ungroup()
df_nfhs_train <- df_nfhs_train |> rowwise() |> mutate(time_id = which(state_year_grids[[state_id]] == year)) |> ungroup()

stan_data_holdout <- list(
  N_srs = nrow(df_srs_train), N_nfhs = nrow(df_nfhs_train), S = S,
  y_srs = log(df_srs_train$tfr_total), state_srs = df_srs_train$state_id, time_id_srs = df_srs_train$time_id,
  y_nfhs = log(df_nfhs_train$tfr_total), state_nfhs = df_nfhs_train$state_id, time_id_nfhs = df_nfhs_train$time_id,
  max_T = max_T, dt = dt, T_s = T_s, tfr_init = tfr_init, tfr_last = tfr_last,
  n_fore = 2023 - 2017, group_of_state = group_of_state,
  prior_timescale_mean = log(50), prior_timescale_sd = 0.5,
  log_F_mean = log(1.2), log_F_sd = 0.30,
  sigma_prior_regime = 1, tight_sigma_stat_sd = 0.15,
  model_type = 2
)

# ── 3. Fit/Load Holdout Model ────────────────────────────────────────────────
fit_file <- here("results", "fit_holdout.rds")

if (file.exists(fit_file)) {
  cat("Loading existing holdout fit from: models/fit_holdout.rds\n")
  fit_holdout <- readRDS(fit_file)
} else {
  cat("Fitting holdout model (data <= 2017)...\n")
  stan_model_file <- here("stan", "tfr_state_model.stan")
  compiled_model <- stan_model(file = stan_model_file)
  
  set.seed(42)
  cat("Running holdout estimation...\n")
  fit_holdout <- sampling(compiled_model, data = stan_data_holdout,
                          iter = 2000, warmup = 1000, chains = 4,
                          control = list(adapt_delta = 0.90, max_treedepth = 10), seed = 42)
# saveRDS(fit_holdout, fit_file)  # NOT SAVED: fit object is ~1-2GB, exceeds GitHub's 100MB file limit. Uncomment this line and run locally to save for downstream scripts.
}

# ── 4. Compile Validation Data Frame (SRS only, >= 2018) ──────────────────────
df_validation <- df_srs_all  |>
  filter(year >= 2018) |>
  select(state_clean, state_id, year, tfr_obs = tfr_total) |>
  mutate(source = "SRS") |>
  left_join(transition_map, by = c("state_clean" = "state_clean")) |>
  arrange(state_clean, year)

# ── 5. Extract Holdout Forecasts, Evaluate Metrics, and Compute TEST CRPS ─────
cat("\nExtracting holdout forecasts...\n")
tfr_fore_holdout <- rstan::extract(fit_holdout, "tfr_fore")$tfr_fore # [draws, S, n_fore]

df_val_results_list <- list()
for (i in 1:nrow(df_validation)) {
  row_val <- df_validation[i, ]
  s       <- row_val$state_id
  yr      <- row_val$year
  tfr_obs <- row_val$tfr_obs
  
  # Forecast step corresponding to the target year
  h <- yr - last_obs_year[s]
  
  if (h >= 1 && h <= dim(tfr_fore_holdout)[3]) {
    draws  <- tfr_fore_holdout[, s, h]
    q_vals <- quantile(draws, probs = c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975))
    
    # TEST CRPS for this held-out state-year, from the full forecast draws
    crps_val <- scoringRules::crps_sample(y = tfr_obs, dat = draws)
    
    df_val_results_list[[i]] <- row_val |>
      mutate(
        model_median   = q_vals["50%"],
        model_mean     = mean(draws),
        model_q2.5     = q_vals["2.5%"],
        model_q5       = q_vals["5%"],
        model_q25      = q_vals["25%"],
        model_q75      = q_vals["75%"],
        model_q95      = q_vals["95%"],
        model_q97.5    = q_vals["97.5%"],
        inside_50      = (tfr_obs >= q_vals["25%"]) & (tfr_obs <= q_vals["75%"]),
        inside_90      = (tfr_obs >= q_vals["5%"]) & (tfr_obs <= q_vals["95%"]),
        model_error    = tfr_obs - model_median,
        crps           = crps_val
      )
  }
}
df_val_results <- bind_rows(df_val_results_list)

# ── 6. Extract Training Fitted Values, Evaluate Metrics, and Compute TRAIN CRPS
cat("Extracting training fits...\n")
tfr_hist_fit_draws <- rstan::extract(fit_holdout, "tfr_hist_fit")$tfr_hist_fit # [draws, S, max_T]

df_srs_train_fitted_list <- list()
for (i in 1:nrow(df_srs_train)) {
  row_train <- df_srs_train[i, ]
  s         <- row_train$state_id
  t         <- row_train$time_id
  tfr_obs   <- row_train$tfr_total
  
  draws <- tfr_hist_fit_draws[, s, t]
  fitted_med <- median(draws)
  
  # TRAIN CRPS: same scoring rule, but against the in-sample fitted draw
  # distribution rather than the forecast distribution. This compares the
  # observed training value against the model's fitted uncertainty at that
  # historical time point.
  crps_train_val <- scoringRules::crps_sample(y = tfr_obs, dat = draws)
  
  df_srs_train_fitted_list[[i]] <- row_train |>
    mutate(
      fitted_median = fitted_med,
      train_error   = tfr_obs - fitted_med,
      crps_train    = crps_train_val
    )
}
df_srs_train_fitted <- bind_rows(df_srs_train_fitted_list) |>
  left_join(transition_map, by = c("state_clean" = "state_clean"))

# ── 7. Generate Holdout Comparison and Summary Tables ─────────────────────────
state_table_csv <- here("results", "holdout_state_comparison.csv")
write_csv(df_val_results |> select(
  state = state_clean, group = group_name, year, source, actual_tfr = tfr_obs,
  forecast_median = model_median, ci_lower_95 = model_q2.5, ci_upper_95 = model_q97.5,
  absolute_error = model_error, crps
), state_table_csv)
cat(sprintf("Saved state-level holdout comparison to: %s\n", state_table_csv))

groups_vec <- c("Pre2001", "By2011", "By2016", "By2021", "NotYet", "Total")
summary_rows <- list()
idx <- 1

for (g in groups_vec) {
  sub_df_test <- if (g == "Total") df_val_results else df_val_results |> filter(group_name == g)
  n_states_sub <- length(unique(sub_df_test$state_id))
  
  sub_df_train <- if (g == "Total") df_srs_train_fitted else df_srs_train_fitted |> filter(group_name == g)
  
  rmse_train <- sqrt(mean(sub_df_train$train_error^2, na.rm = TRUE))
  mae_train  <- mean(abs(sub_df_train$train_error), na.rm = TRUE)
  crps_train <- mean(sub_df_train$crps_train, na.rm = TRUE)
  
  rmse_test <- sqrt(mean(sub_df_test$model_error^2, na.rm = TRUE))
  mae_test  <- mean(abs(sub_df_test$model_error), na.rm = TRUE)
  crps_test <- mean(sub_df_test$crps, na.rm = TRUE)
  
  summary_rows[[idx]] <- tibble(
    Group = g, N_states = n_states_sub,
    Train_RMSE = round(rmse_train, 3), Train_MAE = round(mae_train, 3), Train_CRPS = round(crps_train, 3),
    Test_RMSE  = round(rmse_test, 3),  Test_MAE  = round(mae_test, 3),  Test_CRPS  = round(crps_test, 3)
  )
  idx <- idx + 1
}

df_summary_table <- bind_rows(summary_rows)
summary_csv <- here("results", "holdout_summary_table.csv")
write_csv(df_summary_table, summary_csv)
cat("\n=== OUT-OF-SAMPLE HOLDOUT PERFORMANCE SUMMARY (Train & Test, incl. CRPS) ===\n")
print(df_summary_table)

# ── 8. Observed vs. Predicted Scatter Plot ────────────────────────────────────
cat("\nGenerating observed vs. predicted validation scatter plot...\n")

df_val_results <- df_val_results |>
  mutate(group_name = factor(group_name, levels = c("Pre2001", "By2011", "By2016", "By2021", "NotYet")))

p_holdout <- ggplot() +
  geom_errorbar(data = df_val_results, aes(x = tfr_obs, ymin = model_q2.5, ymax = model_q97.5, color = group_name),
                alpha = 0.35, width = 0.02) +
  geom_point(data = df_val_results, aes(x = tfr_obs, y = model_median, color = group_name),
             size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  scale_color_brewer(palette = "Set1", name = "Transition Group") +
  scale_x_continuous(limits = c(0.8, 3.5), breaks = seq(1.0, 3.5, 0.5)) +
  scale_y_continuous(limits = c(0.8, 3.5), breaks = seq(1.0, 3.5, 0.5)) +
  labs(
    x = "Actual Held-Out TFR (SRS)",
    y = "Predicted TFR (Model Forecast Median)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_holdout

ggsave(here("figures", "04_holdout", "holdout_observed_vs_predicted.png"), plot = p_holdout, width = 8, height = 7, dpi = 300, bg = "white")

cat("\n[SUCCESS] Holdout validation complete.\n")
