# ==============================================================================
# 01_prepare_stan_data.R
#
# Reads the cleaned CSV datasets (SRS & NFHS), aligns them with the original
# hardcoded NFHS-round-based transition groupings, and constructs the 
# state-year time grids. Finally, it exports prep_objects.rds for modeling.
# ==============================================================================

library(tidyverse)
library(here)

cat("--- STAN DATA PREPARATION ---\n")

# ── 1. Load Cleaned Data ───────────────────────────────────────────────────────
srs_path <- here("data", "cleaned_srs.csv")
nfhs_path <- here("data", "cleaned_nfhs.csv")

if (!file.exists(srs_path)) stop("Cleaned SRS data not found. Run 00 first.")
if (!file.exists(nfhs_path)) stop("Cleaned NFHS data not found. Run 00 first.")

df_srs <- read_csv(srs_path, show_col_types = FALSE) |> filter(state_clean != "India")
df_nfhs <- read_csv(nfhs_path, show_col_types = FALSE)

srs_states  <- unique(df_srs$state_clean)
nfhs_states <- unique(df_nfhs$state_clean)
all_states  <- sort(unique(c(srs_states, nfhs_states)))
S <- length(all_states)

cat(sprintf("Loaded successfully: %d States | %d SRS obs | %d NFHS obs\n", S, nrow(df_srs), nrow(df_nfhs)))

if (min(df_srs$tfr_total) <= 0 || min(df_nfhs$tfr_total) <= 0) {
  stop("CRITICAL ERROR: non-positive TFR detected.")
}

# ── 2. Transition Group Mapping ────────────────────────────────────────────────
cat("Applying transition group mapping (NFHS-round-based)...\n")
transition_map <- tribble(
  ~state_clean,                   ~group_id, ~group_name,
  "Goa",                          1L, "Pre2001",
  "Kerala",                       1L, "Pre2001",
  "Puducherry",                   1L, "Pre2001",
  "Tamil Nadu",                   1L, "Pre2001",
  "Andaman & Nicobar Islands",    2L, "By2011",
  "Andhra Pradesh",               2L, "By2011",
  "Chandigarh",                   2L, "By2011",
  "Himachal Pradesh",             2L, "By2011",
  "Karnataka",                    2L, "By2011",
  "Lakshadweep",                  2L, "By2011",
  "Punjab",                       2L, "By2011",
  "Sikkim",                       2L, "By2011",
  "Telangana",                    2L, "By2011",
  "West Bengal",                  2L, "By2011",
  "Arunachal Pradesh",            3L, "By2016",
  "Delhi",                        3L, "By2016",
  "Gujarat",                      3L, "By2016",
  "Haryana",                      3L, "By2016",
  "Jammu & Kashmir",              3L, "By2016",
  "Maharashtra",                  3L, "By2016",
  "Odisha",                       3L, "By2016",
  "Tripura",                      3L, "By2016",
  "Uttarakhand",                  3L, "By2016",
  "Ladakh",                       3L, "By2016",
  "Assam",                        4L, "By2021",
  "Chhattisgarh",                 4L, "By2021",
  "Dadra & Nagar Haveli and Daman & Diu", 4L, "By2021",
  "Madhya Pradesh",               4L, "By2021",
  "Mizoram",                      4L, "By2021",
  "Nagaland",                     4L, "By2021",
  "Rajasthan",                    4L, "By2021",
  "Bihar",                        5L, "NotYet",
  "Jharkhand",                    5L, "NotYet",
  "Manipur",                      5L, "NotYet",
  "Meghalaya",                    5L, "NotYet",
  "Uttar Pradesh",                5L, "NotYet"
)

unmatched_states <- setdiff(all_states, transition_map$state_clean)
if (length(unmatched_states) > 0) {
  stop(sprintf("CRITICAL ERROR: unmatched states: %s", paste(unmatched_states, collapse = ", ")))
}

# ── 3. State/Time Grid Construction ────────────────────────────────────────────
cat("Constructing state-year time grids...\n")
state_map <- tibble(state = all_states, state_id = seq_len(S))
df_srs  <- df_srs  |> left_join(state_map, by = c("state_clean" = "state"))
df_nfhs <- df_nfhs |> left_join(state_map, by = c("state_clean" = "state"))

state_map_with_groups <- state_map |> left_join(transition_map, by = c("state" = "state_clean"))
group_of_state <- state_map_with_groups$group_id

state_year_grids <- list()
T_s <- numeric(S); tfr_init <- numeric(S); tfr_last <- numeric(S); last_obs_year <- numeric(S)

for (s in 1:S) {
  yrs_srs  <- df_srs  |> filter(state_id == s) |> pull(year)
  yrs_nfhs <- df_nfhs |> filter(state_id == s) |> pull(year)
  comb_yrs <- sort(unique(c(yrs_srs, yrs_nfhs)))
  state_year_grids[[s]] <- comb_yrs
  T_s[s] <- length(comb_yrs)

  sub_srs  <- df_srs  |> filter(state_id == s) |> arrange(year)
  sub_nfhs <- df_nfhs |> filter(state_id == s) |> arrange(year)
  tfr_init[s] <- if (nrow(sub_srs) > 0) sub_srs$tfr_total[1] else sub_nfhs$tfr_total[1]

  last_yr <- max(comb_yrs)
  last_obs_year[s] <- last_yr
  last_obs_srs <- df_srs |> filter(state_id == s, year == last_yr)
  if (nrow(last_obs_srs) > 0) {
    tfr_last[s] <- last_obs_srs$tfr_total[1]
  } else {
    last_obs_nfhs <- df_nfhs |> filter(state_id == s, year == last_yr)
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

df_srs <- df_srs |> rowwise() |> mutate(time_id = which(state_year_grids[[state_id]] == year)) |> ungroup()
df_nfhs <- df_nfhs |> rowwise() |> mutate(time_id = which(state_year_grids[[state_id]] == year)) |> ungroup()

# ── 4. Export Preparations ─────────────────────────────────────────────────────
prep_objects <- list(
  S = S,
  max_T = max_T,
  T_s = T_s,
  dt = dt,
  group_of_state = group_of_state,
  tfr_init = tfr_init,
  tfr_last = tfr_last,
  last_obs_year = last_obs_year,
  state_year_grids = state_year_grids,
  df_srs = df_srs,
  df_nfhs = df_nfhs,
  state_map = state_map,
  transition_map = transition_map
)

out_file <- here("results", "prep_objects.rds")
saveRDS(prep_objects, out_file)
cat("[SUCCESS] Saved data structure to:", out_file, "\n")
