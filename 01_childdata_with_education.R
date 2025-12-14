# =======================================
# AQRD Final Project 
# 01_childdata_with_education.R
# Purpose (paper-ready, minimal):
# - Read NIDS Wave 5 child-level microdata and household-derived data.
# - Attach district identifier (2011 boundary) to children via household ID.
# - Construct core education variables used in the paper (age, grade gap,
#   pre-primary attendance, and age at Grade 1 entry).
# - Restrict to the analysis sample (school-age children with valid birth weight,
#   observed schooling outcomes, and required fixed-effect identifiers).
# - Estimate OLS regressions with district and birth-year fixed effects and
#   print regression summaries (Figure 3, 4).
# =======================================

#### 0) Packages ---------------------------------------------------------------
library(tidyverse)
library(haven)
library(stringr)
library(fixest)

#### 1) Read data --------------------------------------------------------------
# NIDS Wave 5 child-level data
child_org <- read_dta("data/child_w5_anon_v1.0.0-stata11.dta")

# NIDS Wave 5 household-derived data (for geography)
hh_org <- read_dta("data/hhderived_w5_anon_v1.0.0-stata11.dta")

#### 2) Attach district identifiers -------------------------------------------
# Attach district (2011 boundary) at the household level
hh_geo <- hh_org |>
  select(w5_hhid, w5_mdbdc2011) |>
  distinct() |>
  rename(w5_dc2011 = w5_mdbdc2011)

# Ensure one-to-one merge at the household level
stopifnot(!any(duplicated(hh_geo$w5_hhid)))

child_org <- child_org |>
  left_join(hh_geo, by = "w5_hhid")

#### 3) Variable construction --------------------------------------------------
dat <- child_org |>
  mutate(
    # Age in 2017
    age_2017 = 2017 - w5_c_dob_y,
    
    # Birth weight (grams and kilograms)
    bw_g = case_when(
      w5_c_brwght_u == 1 ~ as.numeric(w5_c_brwght) * 1000,  # kg → g
      w5_c_brwght_u == 2 ~ as.numeric(w5_c_brwght),         # already in g
      TRUE ~ NA_real_
    ),
    bw_kg = bw_g / 1000,
    
    # Female dummy
    female = case_when(
      w5_c_gen == 2 ~ 1,
      w5_c_gen == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Pre-primary attendance
    preprimary = case_when(
      w5_c_edpre_pp == 1 ~ 1,
      w5_c_edpre_pp == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Age at Grade 1 entry (restricted to plausible range)
    # Note: relies on w5_c_edgrd1yr being present in the raw data
    age_started_raw = w5_c_edgrd1yr - w5_c_dob_y,
    age_started = if_else(
      !is.na(age_started_raw) & age_started_raw >= 5 & age_started_raw <= 10,
      age_started_raw,
      NA_real_
    ),
    
    # Current grade in 2017 (robust parsing)
    # Accepts strings like "Grade 1", "Grade 1 ", "Grade 1 (completed)", etc.
    grade_lev = as_factor(w5_c_ed17curlev),
    grade_chr = str_squish(as.character(grade_lev)),
    grade_current = case_when(
      str_detect(grade_chr, "^Grade\\s+[0-9]{1,2}\\b") ~
        as.numeric(str_match(grade_chr, "^Grade\\s+([0-9]{1,2})\\b")[, 2]),
      TRUE ~ NA_real_
    ),
    
    # Age-for-grade gap
    grade_gap = if_else(
      !is.na(grade_current),
      age_2017 - (grade_current + 6),
      NA_real_
    ),
    
    # Birth year fixed effect
    birth_year = w5_c_dob_y
  )

#### 4) Sample selection -------------------------------------------------------
# Analysis sample:
# - School-age children (6–15)
# - Valid birth weight
# - Observed grade gap (uses grade_current internally)
# - District + birth_year observed (needed for FE regression)
dat_ps3 <- dat |>
  filter(
    age_2017 >= 6, age_2017 <= 15,
    !is.na(bw_g), bw_g >= 500, bw_g <= 6000,
    !is.na(grade_gap),
    !is.na(w5_dc2011),
    !is.na(birth_year),
    grade_current >= 1, grade_current <= 12
  )
# NOTE: household-based restriction removed (no n()>=2 filter)

#### 5) OLS with district and birth-year fixed effects --------------------------
# Figure3 Baseline OLS
ols_gap_base <- feols(
  grade_gap ~ bw_kg + female + preprimary |
    w5_dc2011 + birth_year,
  data = dat_ps3,
  cluster = ~ w5_dc2011
)

# Figure3 Controlling for age at school entry
ols_gap_age <- feols(
  grade_gap ~ bw_kg + female + preprimary + age_started |
    w5_dc2011 + birth_year,
  data = dat_ps3,
  cluster = ~ w5_dc2011
)

# Figure4: School entry timing (mechanism)
ols_age_start <- feols(
  age_started ~ bw_kg + female |
    w5_dc2011 + birth_year,
  data = dat_ps3,
  cluster = ~ w5_dc2011
)

summary(ols_gap_base)
summary(ols_gap_age)
summary(ols_age_start)
