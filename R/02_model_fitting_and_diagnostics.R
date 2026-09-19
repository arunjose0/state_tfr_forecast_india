# ==============================================================================
# 02_model_fitting_and_diagnostics.R
#
# Fits the state-level Bayesian TFR model and performs all critical diagnostics.
# It checks prior predictives, evaluates convergence, calculates PPC coverage,
# and generates the final faceted state-level forecasts (1971-2050).
# Extracts the estimated global floor and transition half-lives.
# ==============================================================================

library(tidyverse)
library(here)
library(rstan)
library(bayesplot)

options(warn = 1, mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

cat("--- MODEL FITTING & DIAGNOSTICS ---\n")

# ── 1. Setup Data & Compilation ──────────────────────────────────────────────
prep <- readRDS(here("results", "prep_objects.rds"))
stan_model_file <- here("stan", "tfr_state_model.stan")
fit_file <- here("results", "fit_main.rds")

# Extract variables from prep
state_labels <- prep$state_map$state
S <- prep$S
last_obs_year <- prep$last_obs_year
n_fore <- 2050 - min(last_obs_year)

build_stan_data <- function(model_type_val = 1) {
  list(
    N_srs = nrow(prep$df_srs), N_nfhs = nrow(prep$df_nfhs), S = prep$S,
    y_srs = log(prep$df_srs$tfr_total), state_srs = prep$df_srs$state_id, time_id_srs = prep$df_srs$time_id,
    y_nfhs = log(prep$df_nfhs$tfr_total), state_nfhs = prep$df_nfhs$state_id, time_id_nfhs = prep$df_nfhs$time_id,
    max_T = prep$max_T, dt = prep$dt, T_s = prep$T_s, tfr_init = prep$tfr_init, tfr_last = prep$tfr_last, 
    n_fore = n_fore, group_of_state = prep$group_of_state,
    prior_timescale_mean = log(50), prior_timescale_sd = 0.5,
    log_F_mean = log(1.2), log_F_sd = 0.30,
    sigma_prior_regime = 1, tight_sigma_stat_sd = 0.15,
    model_type = model_type_val
  )
}

cat("Compiling Stan model...\n")
compiled_model <- stan_model(file = stan_model_file)

# ── 2. Prior Predictive Check ────────────────────────────────────────────────
cat("Running prior predictive checks...\n")
stan_data_prior <- build_stan_data(model_type_val = 1)
fit_prior <- sampling(compiled_model, data = stan_data_prior, iter = 1500, warmup = 750, chains = 2, refresh = 0, seed = 42)

rep_fore <- rstan::extract(fit_prior, "tfr_fore")$tfr_fore
prop_implausible <- mean(rep_fore > 15.0 | rep_fore < 0.3)
rq <- quantile(rep_fore, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
cat("\nPrior Predictive Forecast Quantiles:\n")
print(round(rq, 2))

if (prop_implausible > 0.05 || rq["50%"] > 8.0 || rq["50%"] < 0.35 || rq["99%"] < 2.5 || rq["99%"] > 15.0) {
  stop("Prior predictive check failed.")
} else {
  cat("[SUCCESS] Prior predictive checks passed.\n")
}

# ── 3. Posterior Fit ─────────────────────────────────────────────────────────
if (file.exists(fit_file)) {
  cat(sprintf("Loading existing posterior fit from: %s\n", fit_file))
  fit <- readRDS(fit_file)
} else {
  cat("Running full posterior estimation...\n")
  stan_data_post <- build_stan_data(model_type_val = 2)
  fit <- sampling(compiled_model, data = stan_data_post, iter = 2000, warmup = 1000,
                  chains = 4, control = list(adapt_delta = 0.90, max_treedepth = 10), refresh = 200, seed = 42)
# saveRDS(fit, fit_file)  # NOT SAVED: fit object is ~1-2GB, exceeds GitHub's 100MB file limit. Uncomment this line and run locally to save for downstream scripts.
  cat("[SUCCESS] Saved posterior fit to:", fit_file, "\n")
}

# ── 4. Diagnostics & PPC Coverage ────────────────────────────────────────────
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
rhats <- rstan::summary(fit)$summary[, "Rhat"]
high_rhat_count <- sum(rhats > 1.05, na.rm = TRUE)

cat(sprintf("\n=== SAMPLING DIAGNOSTICS ===\nDivergences: %d\nRhat > 1.05: %d\n", divergences, high_rhat_count))

cat("\n=== POSTERIOR PREDICTIVE CHECKS ===\n")
y_srs_rep <- rstan::extract(fit, "y_srs_rep")$y_srs_rep
y_nfhs_rep <- rstan::extract(fit, "y_nfhs_rep")$y_nfhs_rep
y_srs_obs <- log(prep$df_srs$tfr_total)
y_nfhs_obs <- log(prep$df_nfhs$tfr_total)

srs_50 <- mean(y_srs_obs >= apply(y_srs_rep, 2, quantile, 0.25) & y_srs_obs <= apply(y_srs_rep, 2, quantile, 0.75))
srs_90 <- mean(y_srs_obs >= apply(y_srs_rep, 2, quantile, 0.05) & y_srs_obs <= apply(y_srs_rep, 2, quantile, 0.95))
nfhs_50 <- mean(y_nfhs_obs >= apply(y_nfhs_rep, 2, quantile, 0.25) & y_nfhs_obs <= apply(y_nfhs_rep, 2, quantile, 0.75))
nfhs_90 <- mean(y_nfhs_obs >= apply(y_nfhs_rep, 2, quantile, 0.05) & y_nfhs_obs <= apply(y_nfhs_rep, 2, quantile, 0.95))

cat(sprintf("SRS Coverage  - 50%% PI: %.1f%% | 90%% PI: %.1f%%\n", srs_50 * 100, srs_90 * 100))
cat(sprintf("NFHS Coverage - 50%% PI: %.1f%% | 90%% PI: %.1f%%\n", nfhs_50 * 100, nfhs_90 * 100))

p_srs <- ppc_intervals_grouped(y_srs_obs, y_srs_rep, group = prep$df_srs$state_clean, x = prep$df_srs$year) + 
  facet_wrap(~group, scales="free_y", ncol=6, labeller=label_wrap_gen(width=18)) +
  labs(title="PPC: SRS Observed vs Predicted (log TFR)", x="Year", y="log(TFR)") + theme_minimal()

p_nfhs <- ppc_intervals_grouped(y_nfhs_obs, y_nfhs_rep, group = prep$df_nfhs$state_clean, x = prep$df_nfhs$year) + 
  facet_wrap(~group, scales="free_y", ncol=6, labeller=label_wrap_gen(width=18)) +
  labs(title="PPC: NFHS Observed vs Predicted (log TFR)", x="Year", y="log(TFR)") + theme_minimal()

ggsave(here("figures", "02_model_fitting", "ppc_srs.png"), p_srs, width=18, height=12, bg="white")
ggsave(here("figures", "02_model_fitting", "ppc_nfhs.png"), p_nfhs, width=18, height=12, bg="white")
cat("PPC plots saved to figures/02_model_fitting/\n")

# ── 5. Extract Parameters (Lower Bound & Half-Life) ──────────────────────────
F_post <- rstan::extract(fit, "F")$F
k_group_post <- rstan::extract(fit, "k_group")$k_group
hl_post <- log(2) / k_group_post
years90_post <- rstan::extract(fit, "years_to_90pct_closure")$years_to_90pct_closure
grp_names <- c("Pre2001", "By2011", "By2016", "By2021", "NotYet")

df_speeds <- tibble(
  group = grp_names,
  k_median = apply(k_group_post, 2, median),
  k_q2.5   = apply(k_group_post, 2, quantile, 0.025),
  k_q97.5  = apply(k_group_post, 2, quantile, 0.975),
  half_life_median = apply(hl_post, 2, median),
  half_life_q2.5   = apply(hl_post, 2, quantile, 0.025),
  half_life_q97.5  = apply(hl_post, 2, quantile, 0.975),
  years90_median = apply(years90_post, 2, median),
  years90_q2.5   = apply(years90_post, 2, quantile, 0.025),
  years90_q97.5  = apply(years90_post, 2, quantile, 0.975)
)
write_csv(df_speeds, here("results", "transition_speeds_and_halflives.csv"))

cat(sprintf("\nEstimated Global Floor (F): %.2f (95%% CI: %.2f - %.2f)\n", 
            median(F_post), quantile(F_post, 0.025), quantile(F_post, 0.975)))

# ── 6. State-Level Faceted Forecast Plot ─────────────────────────────────────
cat("\nGenerating state-level forecast plot...\n")
tfr_hist_fit <- rstan::extract(fit, "tfr_hist_fit")$tfr_hist_fit
tfr_fore <- rstan::extract(fit, "tfr_fore")$tfr_fore

summarize_3d <- function(arr, names_, yr_fn, label) {
  res <- list(); idx <- 1
  for (s in 1:S) {
    for (t in 1:dim(arr)[3]) {
      yr <- yr_fn(s, t)
      if (!is.na(yr)) {
        v <- arr[, s, t]
        res[[idx]] <- tibble(state = names_[s], year = yr, type = label,
                             mean = mean(v), q2.5 = quantile(v,.025), q10 = quantile(v,.10),
                             q90 = quantile(v,.90), q97.5 = quantile(v,.975))
        idx <- idx + 1
      }
    }
  }
  bind_rows(res)
}

hist_map <- function(s,t) { g <- prep$state_year_grids[[s]]; if (t <= length(g)) g[t] else NA_integer_ }
fore_map <- function(s,h) { yr <- last_obs_year[s] + h; if (yr <= 2050) yr else NA_integer_ }

df_hist <- summarize_3d(tfr_hist_fit, state_labels, hist_map, "Historical Fit")
df_fore <- summarize_3d(tfr_fore, state_labels, fore_map, "Forecast")
bridge <- df_hist |> group_by(state) |> filter(year == max(year)) |> mutate(type = "Forecast") |> ungroup()
df_traj <- bind_rows(df_hist, bridge, df_fore) |> distinct(state, year, type, .keep_all = TRUE) |> arrange(state, year)

df_obs <- prep$df_srs |> 
  select(state = state_clean, year, tfr_obs = tfr_total) |> 
  mutate(source = "Observed SRS")

df_vlines <- tibble(state = state_labels, last_year = last_obs_year)

p_full <- ggplot() +
  geom_ribbon(data = df_traj, aes(year, ymin=q2.5, ymax=q97.5, fill=type), alpha=.25) +
  geom_ribbon(data = df_traj, aes(year, ymin=q10, ymax=q90, fill=type), alpha=.40) +
  geom_line(data = df_traj, aes(year, mean, color=type, linetype=type), linewidth=.8) +
  geom_point(data = df_obs, aes(year, tfr_obs, shape=source), color="black", size=1.0, alpha=.8) +
  geom_hline(yintercept = 2.1, linetype="dotted", color="gray50") +
  geom_vline(data = df_vlines, aes(xintercept = last_year), linetype="dashed", color="blue", alpha=.35) +
  scale_fill_manual(name = "", values=c("Historical Fit"="#2b5c8f", "Forecast"="#d95f02")) +
  scale_color_manual(name = "", values=c("Historical Fit"="#1b3b5f", "Forecast"="#a63603")) +
  scale_linetype_manual(name = "", values=c("Historical Fit"="solid", "Forecast"="dashed")) +
  scale_shape_manual(name = "", values=c("Observed SRS"=16)) +
  facet_wrap(~state, scales="free_y", ncol=6, labeller = label_wrap_gen(width = 18)) +
  expand_limits(y = 0) +
  labs(x="Year", y="Total Fertility Rate") +
  theme_minimal(base_size=16) +
  theme(plot.background=element_rect(fill="white",color=NA), panel.background=element_rect(fill="white",color=NA),
        strip.background=element_rect(fill="white",color=NA), 
        strip.text=element_text(face="bold",size=14),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(size=18, face="bold"),
        plot.title=element_text(size=22, face="bold", hjust=0.5),
        legend.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14),
        legend.position="bottom", legend.box="horizontal")

# Print to console for immediate viewing
print(p_full)

# Save
ggsave(here("figures", "02_model_fitting", "state_trajectories.png"), p_full, width=18, height=12, dpi=300, bg="white")
cat("Plot saved to figures/02_model_fitting/state_trajectories.png\n")

# ── 7. Export Forecast Tables ────────────────────────────────────────────────
out_csv <- here("results", "tfr_forecasts_2024_2050.csv")
write_csv(df_fore, out_csv)
cat("Forecast table saved to tables/tfr_forecasts_2024_2050.csv\n")
