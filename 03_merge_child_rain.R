# =======================================================
# AQRD Final Project
# 03_merge_child_rain.R
#
# Purpose (paper-ready, minimal):
# - Merge NIDS Wave 5 child-level microdata with pre-constructed
#   district × month CHIRPS rainfall data.
# - Construct prenatal rainfall exposure during pregnancy, focusing on
#   trimester 2 as the main exposure window.
# - Define rainfall shocks as deviations from district-specific
#   calendar-month climatological means.
# - Restrict the analysis sample to rural households
#   (GeoType 2011: Traditional areas or Farms).
# - Produce only the minimal set of outputs required for the paper:
#   (i) analysis dataset,
#   (ii) two diagnostic figures,
#   (iii) first-stage regression output.
#
# Outputs:
#   output/dat_child_trimester_rural.rds
#   output/fig_2.png   (Histogram of trimester-2 rainfall shock)
#   output/fig_6.png   (Scatter: trimester-2 shock vs birth weight)
# =======================================================

#### 0) Packages ---------------------------------------------------------------
# Core data manipulation, estimation, and date-handling packages
library(tidyverse)
library(haven)
library(fixest)
library(lubridate)

#### 1) Load household & child data -------------------------------------------
# Household-derived file:
#   provides geographic identifiers (district code, GeoType 2011).
hh_org <- read_dta("data/hhderived_w5_anon_v1.0.0-stata11.dta")

# Child-level file:
#   provides child characteristics and birth information.
child_org <- read_dta("data/child_w5_anon_v1.0.0-stata11.dta")

# Retain only variables required for rainfall exposure construction
# and the first-stage regression.
vars_child <- c(
  "w5_hhid", "pid",
  "w5_c_gen", "w5_c_dob_y", "w5_c_dob_m",
  "w5_c_brwght", "w5_c_brwght_u"
)

dat_child <- child_org |>
  select(all_of(vars_child))

#### 2) Attach geography (district 2011, GeoType 2011) -------------------------
# w5_mdbdc2011: District Municipality code (2011 boundary).
# w5_geo2011   : GeoType (2011 Census):
#                1 = Traditional areas, 2 = Urban areas, 3 = Farms.
hh_geo <- hh_org |>
  select(w5_hhid, w5_mdbdc2011, w5_geo2011) |>
  distinct() |>
  rename(w5_dc2011 = w5_mdbdc2011)

# Merge household geography into child-level data and construct
# basic identifiers used throughout the script.
child <- dat_child |>
  left_join(hh_geo, by = "w5_hhid") |>
  mutate(
    id_child = row_number(),
    byear    = as.integer(w5_c_dob_y),
    bmonth   = as.integer(w5_c_dob_m),
    # Indicator for rural residence (used for sample restriction).
    rural    = as.integer(w5_geo2011 %in% c(1, 3))
  )

#### 3) Load rainfall (district × month, millimeters) --------------------------
# Pre-constructed district-level monthly rainfall panel derived from CHIRPS.
# Expected variables:
#   w5_dc2011, year, month, rain_mm
rain_dc <- read_rds("output/rain_district_monthly_mm.rds")

#### 4) Monthly climatology (district × calendar-month mean) ------------------
# Define "normal" rainfall as the long-run mean rainfall for each
# district and calendar month.
# Rainfall shocks are constructed as deviations from this baseline.
rain_clim_monthly <- rain_dc |>
  group_by(w5_dc2011, month) |>
  summarise(
    clim_mm = mean(rain_mm, na.rm = TRUE),
    .groups = "drop"
  )

#### 5) Birth date and pregnancy-month expansion ------------------------------
# Approximate the child's date of birth as the first day of the reported
# birth month (day information is not available in the survey).
child <- child |>
  mutate(
    birth_date = as.Date(paste0(byear, "-", bmonth, "-01"))
  )

# Expand each child into a sequence of 12 months ending at the birth month.
# Offset indexing:
#   offset = 0  corresponds to birth_date - 11 months,
#   offset = 11 corresponds to the birth month.
#
# Prenatal exposure is defined over months -9 to -2 relative to birth,
# which maps to offsets 2–9 in this construction.
preg_long <- child |>
  rowwise() |>
  mutate(
    preg_months = list(
      tibble(
        offset = 0:11,
        date   = birth_date %m-% months(11 - offset),
        year   = year(date),
        month  = month(date)
      )
    )
  ) |>
  unnest(preg_months) |>
  ungroup()

# Merge rainfall data and climatology; construct rainfall shocks in millimeters.
preg_rain <- preg_long |>
  left_join(rain_dc, by = c("w5_dc2011", "year", "month")) |>
  left_join(rain_clim_monthly, by = c("w5_dc2011", "month")) |>
  mutate(
    shock_mm = rain_mm - clim_mm
  )

#### 6) Trimester construction (months -9 to -2 only) --------------------------
# Assign pregnancy months to trimesters within the exposure window:
#   offsets 2–4 -> trimester 1,
#   offsets 5–7 -> trimester 2 (main exposure used in the paper),
#   offsets 8–9 -> trimester 3.
trimester_data <- preg_rain |>
  filter(offset %in% 2:9) |>
  mutate(
    trimester = case_when(
      offset %in% 2:4 ~ 1L,
      offset %in% 5:7 ~ 2L,
      offset %in% 8:9 ~ 3L
    )
  ) |>
  group_by(id_child, trimester) |>
  summarise(
    rain_level  = mean(rain_mm,  na.rm = TRUE),
    shock_level = mean(shock_mm, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    id_cols     = id_child,
    names_from  = trimester,
    values_from = c(rain_level, shock_level),
    names_glue  = "{.value}_t{trimester}"
  )

# Retain only trimester-2 exposure variables used in the analysis.
child_trim <- child |>
  left_join(trimester_data, by = "id_child") |>
  rename(
    rain_trim2  = rain_level_t2,
    shock_trim2 = shock_level_t2
  )

#### 7) Birth weight in grams --------------------------------------------------
# Convert reported birth weight into grams using unit codes:
#   1 = kilograms, 2 = grams.
child_trim <- child_trim |>
  mutate(
    bw_g = case_when(
      w5_c_brwght_u == 1 ~ as.numeric(w5_c_brwght) * 1000,
      w5_c_brwght_u == 2 ~ as.numeric(w5_c_brwght),
      TRUE               ~ NA_real_
    )
  )

#### 8) Analysis sample (rural only) -------------------------------------------
# Final estimation sample:
#   - Rural households only,
#   - Non-missing birth weight and trimester-2 rainfall shock,
#   - Valid district and birth-year identifiers.
analysis_dat <- child_trim |>
  filter(
    rural == 1,
    !is.na(bw_g),
    !is.na(shock_trim2),
    !is.na(w5_c_gen),
    !is.na(w5_dc2011),
    !is.na(w5_c_dob_y)
  )

dir.create("output", showWarnings = FALSE)
saveRDS(analysis_dat, "output/dat_child_trimester_rural.rds")

#### 9) Required figures -------------------------------------------------------
# Figure 2: Histogram of trimester-2 rainfall shocks.
p_hist <- ggplot(analysis_dat, aes(x = shock_trim2)) +
  geom_histogram()

ggsave("output/fig_2.png", p_hist)

# Figure 6: Scatter plot of trimester-2 rainfall shock vs birth weight.
p_scatter <- ggplot(analysis_dat, aes(x = shock_trim2, y = bw_g)) +
  geom_point()

ggsave("output/fig_6.png", p_scatter)

#### 11) First-stage regression -----------------------------------------------
# First-stage specification:
#   BirthWeight = beta * RainShock_T2 + gamma * Female
#                 + district FE + birth-year FE + error.
mod_t2_shock <- feols(
  bw_g ~ shock_trim2 + w5_c_gen |
    w5_dc2011 + w5_c_dob_y,
  data    = analysis_dat,
  cluster = ~ w5_dc2011
)

summary(mod_t2_shock)
