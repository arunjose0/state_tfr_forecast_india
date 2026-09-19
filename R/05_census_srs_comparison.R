# ==============================================================================
# 05_census_srs_comparison.R
#
# Validates model estimates and raw SRS period averages against Census
# indirect TFR estimates (2001 and 2011) at the state level.
# 
# Generates a 4-panel scatter plot (Census 2001/2011 vs Model/Raw SRS) and
# outputs the regression summaries and comparison tables.
# ==============================================================================

library(tidyverse)
library(here)
library(rstan)

cat("--- CENSUS VS SRS/MODEL VALIDATION ---\n")

# ── 1. Load Setup & Model ────────────────────────────────────────────────────
prep <- readRDS(here("results", "prep_objects.rds"))
df_srs <- prep$df_srs
state_year_grids <- prep$state_year_grids
all_states <- prep$state_map$state

# Load main fit to get posterior TFR draws
fit_file <- here("results", "fit_main.rds")
if (!file.exists(fit_file)) stop("fit_main.rds not found. Ensure model is fitted.")
fit_main <- readRDS(fit_file)
log_tfr_draws <- rstan::extract(fit_main, "x")$x
n_draws <- dim(log_tfr_draws)[1]

# ── 2. Load and Clean Census Data ────────────────────────────────────────────
census_csv_path <- here("data", "census_tfr_2001_2011.csv")
if (!file.exists(census_csv_path)) stop("Census data not found.")
census_tfr <- read_csv(census_csv_path, show_col_types = FALSE)

clean_state_names_census <- function(st) {
  st |>
    str_trim() |>
    recode(
      "Jammu and Kashmir"                  = "Jammu & Kashmir",
      "Andaman and Nicobar"                = "Andaman & Nicobar Islands",
      "Andaman and Nicobar "               = "Andaman & Nicobar Islands",
      "Andaman and Nicobar Islands"        = "Andaman & Nicobar Islands",
      "Chhatisgarh"                        = "Chhattisgarh",
      "Pondicherry"                        = "Puducherry",
      "Puducherry "                        = "Puducherry",
      "Chandigarh "                        = "Chandigarh",
      "Lakshadweep "                       = "Lakshadweep",
      "Ladakh "                            = "Ladakh",
      "Himachal"                           = "Himachal Pradesh",
      "Telengana"                          = "Telangana",
      "Uttaranchal"                        = "Uttarakhand",
      "Orissa"                             = "Odisha"
    )
}

census_long <- census_tfr |>
  pivot_longer(cols = c(tfr_2001, tfr_2011),
               names_to = "census_year_col", values_to = "census_tfr") |>
  mutate(
    census_year = as.integer(gsub("tfr_", "", census_year_col)),
    state_clean = clean_state_names_census(state),
    log_census_tfr = log(census_tfr)
  ) |>
  filter(!is.na(census_tfr)) |>
  filter(!state_clean %in% c("Dadra and Nagar Haveli", "Daman and Diu", "India"))

# ── 3. Define Period Windows ─────────────────────────────────────────────────
window_defs <- list(
  "2001" = list(census_year = 2001L, window_start = 1994L, window_end = 2000L),
  "2011" = list(census_year = 2011L, window_start = 2004L, window_end = 2010L)
)

# ── 4. Extract Comparisons ───────────────────────────────────────────────────
cat("Extracting period-averages for comparison...\n")
comparison_rows <- list()
idx <- 0

for (wname in names(window_defs)) {
  wdef <- window_defs[[wname]]
  census_yr <- wdef$census_year
  w_start   <- wdef$window_start
  w_end     <- wdef$window_end
  
  census_sub <- census_long |> filter(census_year == census_yr)
  
  for (i in 1:nrow(census_sub)) {
    row <- census_sub[i, ]
    grp_info <- prep$transition_map |> filter(state_clean == row$state_clean)
    grp_name <- if (nrow(grp_info) > 0) grp_info$group_name[1] else "Unknown"
    
    # Model: period-average
    s <- which(all_states == row$state_clean)
    if (length(s) > 0) {
      comb_yrs <- state_year_grids[[s]]
      window_indices <- which(comb_yrs >= w_start & comb_yrs <= w_end)
      
      if (length(window_indices) > 0) {
        draw_matrix <- log_tfr_draws[, s, window_indices, drop = FALSE]
        dim(draw_matrix) <- c(n_draws, length(window_indices))
        period_avg_draws <- rowMeans(draw_matrix)
        
        q <- quantile(period_avg_draws, probs = c(0.025, 0.50, 0.975))
        
        idx <- idx + 1
        comparison_rows[[idx]] <- tibble(
          state = row$state_clean, group = grp_name, census_year = census_yr,
          log_census_tfr = row$log_census_tfr, comparison_type = "Model",
          fitted_median = q["50%"], ci_lo = q["2.5%"], ci_hi = q["97.5%"],
          abs_diff = abs(q["50%"] - row$log_census_tfr)
        )
      }
    }
    
    # Raw SRS: period-average
    srs_window <- df_srs |>
      filter(state_clean == row$state_clean, year >= w_start, year <= w_end)
    
    n_srs <- nrow(srs_window)
    if (n_srs > 0) {
      log_srs_vals <- log(srs_window$tfr_total)
      raw_mean <- mean(log_srs_vals)
      
      if (n_srs >= 2) {
        raw_se <- sd(log_srs_vals) / sqrt(n_srs)
        raw_ci_lo <- raw_mean - 1.96 * raw_se
        raw_ci_hi <- raw_mean + 1.96 * raw_se
      } else {
        raw_ci_lo <- NA_real_
        raw_ci_hi <- NA_real_
      }
      
      idx <- idx + 1
      comparison_rows[[idx]] <- tibble(
        state = row$state_clean, group = grp_name, census_year = census_yr,
        log_census_tfr = row$log_census_tfr, comparison_type = "Raw SRS",
        fitted_median = raw_mean, ci_lo = raw_ci_lo, ci_hi = raw_ci_hi,
        abs_diff = abs(raw_mean - row$log_census_tfr)
      )
    }
  }
}

df_comparison <- bind_rows(comparison_rows) |>
  mutate(
    census_year_label = paste("Census", census_year),
    comparison_type = factor(comparison_type, levels = c("Model", "Raw SRS")),
    group = factor(group, levels = c("Pre2001", "By2011", "By2016", "By2021", "NotYet"))
  )

write_csv(df_comparison, here("results", "census_vs_model_srs.csv"))

# ── 5. Regression Summaries ──────────────────────────────────────────────────
cat("\n=== REGRESSION SUMMARY ===\n")
reg_rows <- list()
reg_idx <- 0
ann_data <- tibble()

for (cy in c(2001, 2011)) {
  for (ctype in c("Model", "Raw SRS")) {
    df_sub <- df_comparison |> filter(census_year == cy, comparison_type == ctype)
    if (nrow(df_sub) < 3) next
    
    fit_lm <- lm(fitted_median ~ log_census_tfr, data = df_sub)
    slope  <- coef(fit_lm)[2]
    intcpt <- coef(fit_lm)[1]
    r_sq   <- summary(fit_lm)$r.squared
    mae    <- mean(df_sub$abs_diff)
    n_val  <- nrow(df_sub)
    
    reg_idx <- reg_idx + 1
    reg_rows[[reg_idx]] <- tibble(
      census_year = cy, comparison_type = ctype, N = n_val,
      MAE = round(mae, 4), OLS_slope = round(slope, 4),
      OLS_intercept = round(intcpt, 4), R_squared = round(r_sq, 4)
    )
    
    cat(sprintf("[%d %s] N=%d | MAE=%.4f | y = %.4fx + %.4f | R^2 = %.4f\n",
                cy, ctype, n_val, mae, slope, intcpt, r_sq))
    
    eq_label <- sprintf("y = %.3fx + %.4f\nR^2 = %.4f", slope, intcpt, r_sq)
    ann_data <- bind_rows(ann_data, tibble(
      census_year_label = paste("Census", cy),
      comparison_type = factor(ctype, levels = c("Model", "Raw SRS")),
      label = eq_label
    ))
  }
}

df_reg_summary <- bind_rows(reg_rows)
write_csv(df_reg_summary, here("results", "census_regression_summary.csv"))

# ── 6. Faceted Scatter Plot ──────────────────────────────────────────────────
cat("\nGenerating 4-panel comparison plot...\n")

all_vals <- c(df_comparison$log_census_tfr, df_comparison$fitted_median, df_comparison$ci_lo, df_comparison$ci_hi)
all_vals <- all_vals[!is.na(all_vals)]
ax_lo <- floor(min(all_vals) * 10) / 10
ax_hi <- ceiling(max(all_vals) * 10) / 10

ann_data <- ann_data |> mutate(x_pos = ax_lo + (ax_hi - ax_lo) * 0.05, y_pos = ax_hi - (ax_hi - ax_lo) * 0.05)

p <- ggplot(df_comparison, aes(x = log_census_tfr, y = fitted_median)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, color = group),
                alpha = 0.35, width = 0.015, linewidth = 0.8) +
  geom_point(aes(color = group), size = 2, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray30", linewidth = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              color = "black", linewidth = 0.7, linetype = "solid") +
  geom_text(data = ann_data, aes(x = x_pos, y = y_pos, label = label),
            hjust = 0, vjust = 1, size = 3.2, fontface = "italic", color = "gray20") +
  scale_color_brewer(palette = "Set1", name = "Transition Group") +
  scale_x_continuous(limits = c(ax_lo, ax_hi)) +
  scale_y_continuous(limits = c(ax_lo, ax_hi)) +
  coord_equal() +
  facet_grid(rows = vars(census_year_label), cols = vars(comparison_type)) +
  labs(
    x = "log(Census Indirect TFR Estimate)",
    y = "log(TFR) Estimate (Model or Raw SRS)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",

    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 11),
    strip.background  = element_rect(fill = "gray90", color = NA)
  )

print(p)
ggsave(here("figures", "05_census_validation", "census_srs_comparison.png"), p, width = 8, height = 7, dpi = 300, bg = "white")

cat("\n[SUCCESS] Census validation script completed.\n")
