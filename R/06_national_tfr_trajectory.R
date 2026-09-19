# ==============================================================================
# 06_national_tfr_trajectory.R
#
# Computes the National (India) TFR trajectory by taking a population-weighted
# average of the 36 modelled states.
# Weights: 2011 Census proportion of women aged 15-49.
# The aggregation is done at the posterior draw level to correctly propagate
# uncertainty to the national level.
# ==============================================================================

library(tidyverse)
library(here)
library(rstan)

cat("--- NATIONAL TFR TRAJECTORY AGGREGATION ---\n")

# ── 1. Load Data & Fits ───────────────────────────────────────────────────────
cat("Loading model fit, prep objects, and population weights...\n")

prep <- readRDS(here("results", "prep_objects.rds"))
fit_file <- here("results", "fit_main.rds")
if (!file.exists(fit_file)) stop("fit_main.rds not found. Ensure model is fitted.")
fit_main <- readRDS(fit_file)

weights_path <- here("data", "state_women_15_49_proportions.csv")
if(!file.exists(weights_path)) stop("Weights CSV not found!")
df_weights <- read_csv(weights_path, show_col_types = FALSE)

state_labels <- prep$state_map$state

# Align weights to the exact state order used in the Stan model
df_weights_aligned <- tibble(state_clean = state_labels) |>
  left_join(df_weights, by = "state_clean")

w <- df_weights_aligned$proportion_of_india
if(any(is.na(w))) stop("Missing weights for some modelled states!")

# ── 2. Extract Posterior Draws ────────────────────────────────────────────────
tfr_hist <- rstan::extract(fit_main, "tfr_hist_fit")$tfr_hist_fit  # [D, S, max_T]
tfr_fore <- rstan::extract(fit_main, "tfr_fore")$tfr_fore          # [D, S, n_fore]

D <- dim(tfr_hist)[1]
S <- length(state_labels)
n_fore_val <- dim(tfr_fore)[3]

# ── 3. Draw-Level Aggregation ─────────────────────────────────────────────────
cat("Aggregating 36 state trajectories to national level...\n")

min_yr <- 1971
max_yr <- 2050
common_years <- min_yr:max_yr
n_years <- length(common_years)

india_draws <- matrix(0, nrow = D, ncol = n_years)

parent_map <- list(
  "Telangana" = "Andhra Pradesh",
  "Ladakh" = "Jammu & Kashmir",
  "Uttarakhand" = "Uttar Pradesh",
  "Jharkhand" = "Bihar",
  "Chhattisgarh" = "Madhya Pradesh"
)

dense_state_draws <- array(0, dim = c(D, S, n_years))

for (s in 1:S) {
  s_hist_years <- prep$state_year_grids[[s]]
  s_last_obs <- prep$last_obs_year[s]
  s_fore_years <- seq(s_last_obs + 1, length.out = n_fore_val)
  
  x_all <- c(s_hist_years, s_fore_years)
  
  for (d in 1:D) {
    y_hist <- tfr_hist[d, s, 1:length(s_hist_years)]
    y_fore <- tfr_fore[d, s, ]
    y_all <- c(y_hist, y_fore)
    
    # rule = 2 extends the earliest value backward
    dense_state_draws[d, s, ] <- approx(x = x_all, y = y_all, xout = common_years, rule = 2)$y
  }
}

# Apply parent imputation for years BEFORE the child state's first observation
for (child in names(parent_map)) {
  parent <- parent_map[[child]]
  c_idx <- which(state_labels == child)
  p_idx <- which(state_labels == parent)
  
  if (length(c_idx) == 1 && length(p_idx) == 1) {
    c_min_yr <- min(prep$state_year_grids[[c_idx]])
    impute_indices <- which(common_years < c_min_yr)
    
    if (length(impute_indices) > 0) {
      dense_state_draws[, c_idx, impute_indices] <- dense_state_draws[, p_idx, impute_indices]
    }
  }
}

# Accumulate the weighted sum for India
for (s in 1:S) {
  india_draws <- india_draws + (dense_state_draws[, s, ] * w[s])
}

india_quantiles <- tibble(
  year = common_years,
  median = apply(india_draws, 2, median),
  q5 = apply(india_draws, 2, quantile, probs = 0.05),
  q95 = apply(india_draws, 2, quantile, probs = 0.95),
  q25 = apply(india_draws, 2, quantile, probs = 0.25),
  q75 = apply(india_draws, 2, quantile, probs = 0.75),
  q2.5 = apply(india_draws, 2, quantile, probs = 0.025),
  q97.5 = apply(india_draws, 2, quantile, probs = 0.975)
)

# ── 4. Load Raw SRS India Data for Validation ─────────────────────────────────
srs_path <- here("data", "cleaned_srs.csv")
df_srs <- read_csv(srs_path, show_col_types = FALSE)

# Extract India aggregate data directly from the cleaned dataset 
# (It was included during cleaning but filtered out in 01_prepare_stan_data for modeling)
df_srs_india <- df_srs |>
  filter(state_clean == "India" & year >= min_yr) |>
  select(year, tfr_total)

# ── 5. Plot ───────────────────────────────────────────────────────────────────
cat("Generating plot...\n")

india_quantiles_hist <- india_quantiles |> filter(year <= 2023)
india_quantiles_fore <- india_quantiles |> filter(year >= 2023)

p_nat <- ggplot() +
  geom_hline(yintercept = 2.1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_ribbon(data = india_quantiles, aes(x = year, ymin = q2.5, ymax = q97.5, fill = "95% CI"), alpha = 0.1) +
  geom_ribbon(data = india_quantiles, aes(x = year, ymin = q5, ymax = q95, fill = "90% CI"), alpha = 0.2) +
  geom_ribbon(data = india_quantiles, aes(x = year, ymin = q25, ymax = q75, fill = "50% CI"), alpha = 0.4) +
  geom_line(data = india_quantiles_hist, aes(x = year, y = median, color = "Model Estimate"), linewidth = 0.8) +
  geom_line(data = india_quantiles_fore, aes(x = year, y = median, color = "Model Estimate"), linewidth = 0.8, linetype = "dashed") +
  geom_point(data = df_srs_india, aes(x = year, y = tfr_total, color = "Raw SRS (India)"), size = 2) +
  geom_line(data = df_srs_india, aes(x = year, y = tfr_total, color = "Raw SRS (India)"), linewidth = 0.6, alpha = 0.5) +
  scale_fill_manual(name = "Model Credible Intervals", values = c("95% CI" = "#2C7BB6", "90% CI" = "#2C7BB6", "50% CI" = "#2C7BB6")) +
  scale_color_manual(name = "TFR Estimate", values = c("Model Estimate" = "#000000", "Raw SRS (India)" = "#D7191C")) +
  scale_x_continuous(breaks = seq(1970, 2050, 10)) +
  scale_y_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) +
  labs(
    x = "Year", y = "Total Fertility Rate (TFR)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 22),
    axis.title = element_text(face = "bold")
  )

print(p_nat)
ggsave(here("figures", "06_national_trajectory", "national_tfr_trajectory.png"), plot = p_nat, width = 12, height = 8, dpi = 300, bg = "white")
write_csv(india_quantiles, here("results", "national_tfr_quantiles.csv"))

