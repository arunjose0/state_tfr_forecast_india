# ==============================================================================
# 03_sensitivity_analysis.R
#
# Sensitivity Analysis:
# Evaluates the robustness of the TFR forecast to two distinct sources of
# uncertainty:
#   (A) Prior uncertainty (Wide vs Tight sigma priors)
#   (B) Floor location (Low vs High theoretical limits)
# ==============================================================================

library(tidyverse)
library(here)
library(rstan)

options(warn = 1, mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

cat("--- SENSITIVITY ANALYSIS (PRIORS & FLOOR) ---\n")

# ── 1. Setup ─────────────────────────────────────────────────────────────────
prep <- readRDS(here("results", "prep_objects.rds"))
stan_model_file <- here("stan", "tfr_state_model.stan")
compiled_model <- stan_model(file = stan_model_file)

states <- prep$state_map$state
last_obs_year <- prep$last_obs_year
n_fore <- 2050 - min(last_obs_year)

stan_data_main <- list(
  N_srs = nrow(prep$df_srs), N_nfhs = nrow(prep$df_nfhs), S = prep$S,
  y_srs = log(prep$df_srs$tfr_total), state_srs = prep$df_srs$state_id, time_id_srs = prep$df_srs$time_id,
  y_nfhs = log(prep$df_nfhs$tfr_total), state_nfhs = prep$df_nfhs$state_id, time_id_nfhs = prep$df_nfhs$time_id,
  max_T = prep$max_T, dt = prep$dt, T_s = prep$T_s, tfr_init = prep$tfr_init, tfr_last = prep$tfr_last, 
  n_fore = n_fore, group_of_state = prep$group_of_state,
  prior_timescale_mean = log(50), prior_timescale_sd = 0.5,
  log_F_mean = log(1.2), log_F_sd = 0.30,
  sigma_prior_regime = 1, tight_sigma_stat_sd = 0.15,
  model_type = 2
)

# ── 2. Define Regimes ────────────────────────────────────────────────────────
stan_data_wide <- stan_data_main
stan_data_wide$sigma_prior_regime <- 2
stan_data_wide$log_F_mean <- log(1.2)
stan_data_wide$log_F_sd <- 0.60

stan_data_tight <- stan_data_main
stan_data_tight$sigma_prior_regime <- 3
stan_data_tight$log_F_mean <- log(1.2)
stan_data_tight$log_F_sd <- 0.15

stan_data_low <- stan_data_main
stan_data_low$sigma_prior_regime <- 1
stan_data_low$log_F_mean <- log(0.6)
stan_data_low$log_F_sd <- 0.30

stan_data_high <- stan_data_main
stan_data_high$sigma_prior_regime <- 1
stan_data_high$log_F_mean <- log(1.8)
stan_data_high$log_F_sd <- 0.30

regimes <- list(
  wide_combined  = stan_data_wide,
  tight_combined = stan_data_tight,
  floor_low      = stan_data_low,
  floor_high     = stan_data_high
)

# ── 3. Fit / Load Models ─────────────────────────────────────────────────────
fits <- list()
fits[["main"]] <- readRDS(here("results", "fit_main.rds"))

for (r_name in names(regimes)) {
  out_path <- here("results", sprintf("fit_%s.rds", r_name))
  if (file.exists(out_path)) {
    cat(sprintf("Loading existing fit: %s\n", r_name))
    fits[[r_name]] <- readRDS(out_path)
  } else {
    cat(sprintf("Fitting: %s...\n", r_name))
    fits[[r_name]] <- sampling(compiled_model, data = regimes[[r_name]], 
                               iter = 2000, warmup = 1000, chains = 4, 
                               control = list(adapt_delta = 0.90, max_treedepth = 10), seed = 42)
# saveRDS(fits[[r_name]], out_path)  # NOT SAVED: fit object is ~1-2GB, exceeds GitHub's 100MB file limit. Uncomment this line and run locally to save for downstream scripts.
  }
}

# ── 4. Forecasts & Intervals ──────────────────────────────────────────────────
cat("Extracting forecasts...\n")
grp_names <- c("Pre2001", "By2011", "By2016", "By2021", "NotYet")
group_of_state <- prep$group_of_state
target_years <- c(2030, 2040, 2050)

extract_forecasts <- function(fit, fit_name) {
  tfr_draws <- rstan::extract(fit, "tfr_fore")$tfr_fore
  df_list <- list()
  for (s_idx in seq_along(states)) {
    s_proj_years <- seq(last_obs_year[s_idx] + 1, length.out = n_fore)
    for (t_yr in target_years) {
      if (t_yr %in% s_proj_years) {
        p_idx <- which(s_proj_years == t_yr)
        draws_yr <- tfr_draws[, s_idx, p_idx]
        q <- quantile(draws_yr, probs = c(0.05, 0.50, 0.95))
        df_list[[length(df_list) + 1]] <- tibble(
          regime = fit_name, state = states[s_idx], year = t_yr,
          median_tfr = q["50%"], ci90_width = q["95%"] - q["5%"],
          is_above_replacement = median_tfr > 2.1
        )
      }
    }
  }
  bind_rows(df_list)
}

forecast_summaries <- map2_df(fits, names(fits), ~extract_forecasts(.x, .y))
main_forecasts <- forecast_summaries |> filter(regime == "main") |>
  select(state, year, main_median = median_tfr, main_width = ci90_width, main_status = is_above_replacement)
group_lookup <- tibble(state = states, group = grp_names[group_of_state])

forecast_comparisons <- forecast_summaries |>
  left_join(main_forecasts, by = c("state", "year")) |>
  left_join(group_lookup, by = "state") |>
  mutate(
    signed_diff_vs_main = median_tfr - main_median,
    abs_diff_vs_main = abs(median_tfr - main_median),
    width_pct_change_vs_main = (ci90_width - main_width) / main_width * 100,
    status_flipped = ifelse(regime != "main" & (is_above_replacement != main_status), TRUE, FALSE)
  )

write_csv(forecast_comparisons |> select(regime, state, group, year, median_tfr, signed_diff_vs_main, abs_diff_vs_main), 
          here("results", "forecast_comparison_5regime.csv"))
write_csv(forecast_comparisons |> select(regime, state, group, year, ci90_width, width_pct_change_vs_main), 
          here("results", "interval_width_comparison_5regime.csv"))

# ── 5. Summary ────────────────────────────────────────────────────────────────
cat("\nSummary of 2050 Projections:\n")
df_2050 <- forecast_comparisons |> filter(year == 2050, regime != "main")

for (r in c("wide_combined", "tight_combined", "floor_low", "floor_high")) {
  sub_df <- df_2050 |> filter(regime == r)
  mean_shift <- mean(sub_df$abs_diff_vs_main)
  max_shift <- max(sub_df$abs_diff_vs_main)
  avg_width_pct <- mean(sub_df$width_pct_change_vs_main)
  flipped <- sub_df |> filter(status_flipped == TRUE) |> pull(state)
  
  cat(sprintf("\nRegime: %s\n", toupper(r)))
  cat(sprintf("  Mean Abs Shift: %.3f | Max Abs Shift: %.3f\n", mean_shift, max_shift))
  cat(sprintf("  Avg CI Width Change: %+.1f%%\n", avg_width_pct))
  cat(sprintf("  Status Flipped: %s\n", if(length(flipped) > 0) paste(flipped, collapse=", ") else "None"))
}

# ── 6. Figures ────────────────────────────────────────────────────────────────
cat("\nGenerating overlay figures...\n")

extract_full_trajectories <- function(fit, fit_name) {
  tfr_draws <- rstan::extract(fit, "tfr_fore")$tfr_fore
  df_list <- list()
  for (s_idx in seq_along(states)) {
    s_proj_years <- seq(last_obs_year[s_idx] + 1, length.out = n_fore)
    for (p_idx in seq_along(s_proj_years)) {
      yr <- s_proj_years[p_idx]
      q <- quantile(tfr_draws[, s_idx, p_idx], probs = c(0.05, 0.50, 0.95))
      df_list[[length(df_list) + 1]] <- tibble(
        regime = fit_name, state = states[s_idx], year = yr,
        median = q["50%"], q5 = q["5%"], q95 = q["95%"]
      )
    }
  }
  bind_rows(df_list)
}

full_traj <- map2_df(fits, names(fits), ~extract_full_trajectories(.x, .y)) |> 
  mutate(regime = factor(regime, levels = c("main", "wide_combined", "tight_combined", "floor_low", "floor_high"))) |> 
  filter(year >= 2024, year <= 2050)

pal_A <- c(main = "#000000", wide_combined = "#1B9E77", tight_combined = "#D95F02")
pal_B <- c(main = "#000000", floor_low = "#7570B3", floor_high = "#E7298A")

plot_theme <- theme_minimal(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(color = "gray30", size = 15),
    plot.caption = element_text(size = 12, color = "gray40", hjust = 0),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray92", color = NA),
    strip.text = element_text(face = "bold", size = 12)
  )

# Figure A
p_compare_A <- full_traj |> filter(regime %in% c("main", "wide_combined", "tight_combined")) |>
  ggplot(aes(x = year, y = median, color = regime, fill = regime)) +
  geom_hline(yintercept = 2.1, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.08, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~state, ncol = 6) +
  expand_limits(y = 0) +
  scale_color_manual(values = pal_A) + scale_fill_manual(values = pal_A) +
  scale_x_continuous(breaks = seq(2025, 2050, 5)) +
  labs(
    title = "Sigma/Floor Uncertainty Sensitivity: All 36 States",
    subtitle = "MAIN vs WIDE-COMBINED vs TIGHT-COMBINED (2024-2050, 90% CI)",
    x = "Year", y = "Total Fertility Rate (TFR)", color = "Prior Specification", fill = "Prior Specification",
    caption = "MAIN: log_F_mean=log(1.2), log_F_sd=0.30, sigma~exponential(2).\nWIDE-COMBINED: log_F_sd=0.60, sigma~student_t(3,0,2.5).\nTIGHT-COMBINED: log_F_sd=0.15, sigma~normal priors."
  ) + plot_theme

ggsave(here("figures", "03_sensitivity", "sigma_floor_uncertainty_all36.png"), plot = p_compare_A, width = 28, height = 32, dpi = 300, bg = "white", limitsize = FALSE)

# Figure B
p_compare_B <- full_traj |> filter(regime %in% c("main", "floor_low", "floor_high")) |>
  ggplot(aes(x = year, y = median, color = regime, fill = regime)) +
  geom_hline(yintercept = 2.1, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.08, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~state, ncol = 6) +
  expand_limits(y = 0) +
  scale_color_manual(values = pal_B) + scale_fill_manual(values = pal_B) +
  scale_x_continuous(breaks = seq(2025, 2050, 5)) +
  labs(
    title = "Floor Location Sensitivity: All 36 States",
    subtitle = "MAIN vs FLOOR-LOW (log_F_mean=log(0.6)) vs FLOOR-HIGH (log_F_mean=log(1.8)) (2024-2050, 90% CI)",
    x = "Year", y = "Total Fertility Rate (TFR)", color = "Prior Specification", fill = "Prior Specification",
    caption = "MAIN: log_F_mean=log(1.2). FLOOR-LOW: log_F_mean=log(0.6). FLOOR-HIGH: log_F_mean=log(1.8)."
  ) + plot_theme

ggsave(here("figures", "03_sensitivity", "floor_location_all36.png"), plot = p_compare_B, width = 28, height = 32, dpi = 300, bg = "white", limitsize = FALSE)

cat("\n[SUCCESS] Sensitivity analysis complete.\n")
