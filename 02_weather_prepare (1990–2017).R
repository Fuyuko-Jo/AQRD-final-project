# =======================================================
# AQRD Final Project
# 02_weather_prepare_r.R (1990–2017)
#
# Purpose (paper-ready, minimal):
# - Prepare monthly CHIRPS rainfall (mm) for South Africa over 1990–2017.
# - Aggregate gridded rainfall to administrative units in the NIDS geography:
#   (i) District Municipalities (2011 boundary) for the main analysis, and
#   (ii) Wards (2011 boundary) for Figure 7 (discussion/illustration only).
# - Produce district × month and ward × month rainfall panels to support
#   downstream construction of pregnancy exposure measures (e.g., trimester shocks).
#
# Notes:
# - The raster layers are assumed to be monthly GeoTIFFs whose filenames contain
#   a "YYYY.MM" pattern used to name layers and construct a date index.
# =======================================================

#### 0) Packages ---------------------------------------------------------------
library(tidyverse)
library(terra)
library(lubridate)

#### 1) Spatial extent (South Africa bounding box) -----------------------------
# Define a coarse bounding box for South Africa to speed up raster processing.
# All CHIRPS rasters are cropped to this extent before administrative aggregation.
ZAF_ext <- ext(15, 33, -35, -22)

#### 2) Load CHIRPS GeoTIFFs (1990–2017) --------------------------------------
# Read all monthly CHIRPS GeoTIFF files and keep only the 1990–2017 window.
# Filenames are expected to include "YYYY.MM", which is used to (i) filter years
# and (ii) name raster layers for later reshaping to a panel dataset.
files_all <- list.files("wdata", pattern = "tif$", full.names = TRUE)

# Extract year from filenames (expects "...YYYY.MM..." somewhere in name)
file_years <- as.integer(gsub(".*(\\d{4})\\.\\d{2}.*", "\\1", files_all))

# Keep only 1990–2017 window
files <- files_all[file_years >= 1990 & file_years <= 2017]

if (length(files) == 0) {
  stop("No CHIRPS TIF files found for 1990–2017 in 'wdata' directory.")
}

# Read + crop each raster, then stack
rain_list <- lapply(files, function(f) {
  r <- rast(f)
  crop(r, ZAF_ext)
})

rain_stack <- rast(rain_list)

# Name layers as "YYYY.MM"
names(rain_stack) <- gsub(".*(\\d{4}\\.\\d{2}).*", "\\1", files)

#### 3) Clean invalid rainfall ------------------------------------------------
# CHIRPS rainfall is non-negative by construction; negative values are treated
# as invalid and set to missing.
rain_stack[rain_stack < 0] <- NA

#### 4) District boundaries and extraction -----------------------------------
# Aggregate gridded rainfall to District Municipalities (2011 boundary).
# The output is a wide table where each column corresponds to a raster month,
# later reshaped into a district × month panel aligned with NIDS district codes.
dc_shp <- vect("data/shp/MDBDistrictMunicipalBoundary2011.gdb")
dc_shp <- project(dc_shp, crs(rain_stack))

# Stable ID for joining extracted values back to attributes
dc_shp$ID <- seq_len(nrow(dc_shp))

# District-average rainfall for each raster layer (month)
# (mean across grid cells intersecting each district polygon)
rain_dc_wide <- terra::extract(
  rain_stack,
  dc_shp,
  fun   = mean,
  na.rm = TRUE,
  ID    = TRUE
) |>
  as.data.frame()

# Bring district code from shapefile attributes
# (used for merges with NIDS district identifiers downstream)
dc_attrs <- as.data.frame(dc_shp)

rain_dc_wide <- rain_dc_wide |>
  left_join(dc_attrs, by = "ID") |>
  rename(w5_dc2011 = DistrictMunicipalityCode)

#### 5) Reshape to district × month panel -------------------------------------
# Reshape wide district-month rainfall into a long panel:
# one row per (district, year, month), with rainfall in millimeters.
# The "YYYY.MM" layer names are converted into a proper date index.
rain_dc_monthly <- rain_dc_wide |>
  select(w5_dc2011, matches("^\\d{4}\\.\\d{2}$")) |>
  pivot_longer(
    cols      = matches("^\\d{4}\\.\\d{2}$"),
    names_to  = "ym",
    values_to = "rain_mm"
  ) |>
  mutate(
    date  = as.Date(paste0(ym, ".01"), format = "%Y.%m.%d"),
    year  = year(date),
    month = month(date)
  ) |>
  arrange(w5_dc2011, year, month)

dir.create("output", showWarnings = FALSE)
saveRDS(rain_dc_monthly, "output/rain_district_monthly_mm.rds")

################################################################################
#### 6) Ward boundaries (MDB 2011) --------------------------------------------
# Construct the Ward × month rainfall panel for Figure 7 only.
# This panel is not used in the main estimation; it is used to illustrate that
# finer administrative units exhibit (potentially) greater dispersion in shocks
# relative to district-level aggregation.
#
# Because ward attribute schemas can vary, the script searches for a plausible
# ward-code field and falls back to the first field containing "ward" or "code"
# if common candidates are not found.

# Read Ward boundaries (File Geodatabase)
ward_shp <- terra::vect("data/shp/MDBWard2011.gdb")
ward_shp <- terra::project(ward_shp, terra::crs(rain_stack))

# Stable internal ID for extraction join
ward_shp$ID <- seq_len(nrow(ward_shp))

# Identify a Ward code field robustly
ward_attrs <- as.data.frame(ward_shp)

cand_ward_fields <- c(
  "WardCode",
  "WARD_CODE",
  "WARDNO",
  "WARD",
  "CODE",
  "CAT_B"
)

ward_code_field <- cand_ward_fields[cand_ward_fields %in% names(ward_attrs)][1]

if (is.na(ward_code_field)) {
  hit <- grep("ward|code", names(ward_attrs), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) {
    stop("Could not find a Ward code field. Inspect names(ward_attrs).")
  } else {
    ward_code_field <- hit[1]
    message("Ward code field not found in common candidates. Using fallback: ", ward_code_field)
  }
} else {
  message("Using Ward code field: ", ward_code_field)
}

# Create clean Ward code
ward_shp$ward_code <- as.character(ward_shp[[ward_code_field]])


#### 8) Ward–average rainfall -------------------------------------------------
# Aggregate gridded rainfall to wards (2011 boundary), analogous to the district
# extraction step. The result is reshaped to a ward × month panel.
rain_ward_wide <- terra::extract(
  rain_stack,
  ward_shp,
  fun   = mean,
  na.rm = TRUE,
  ID    = TRUE
) |>
  as.data.frame() |>
  left_join(
    ward_attrs |>
      transmute(
        ID = ID,
        ward_code = as.character(.data[[ward_code_field]])
      ),
    by = "ID"
  )

#### 9) Reshape to Ward × month panel -----------------------------------------
# Convert the wide ward extraction output into a long panel indexed by
# (ward_code, year, month).
rain_ward_monthly <- rain_ward_wide |>
  select(ward_code, matches("^\\d{4}\\.\\d{2}$")) |>
  pivot_longer(
    cols      = matches("^\\d{4}\\.\\d{2}$"),
    names_to  = "ym",
    values_to = "rain_mm"
  ) |>
  mutate(
    date  = as.Date(paste0(ym, ".01"), format = "%Y.%m.%d"),
    year  = lubridate::year(date),
    month = lubridate::month(date)
  ) |>
  arrange(ward_code, year, month)

# Sanity check
# Warn if any ward codes are missing after the join (should be rare and indicates
# an attribute/join issue).
if (any(is.na(rain_ward_monthly$ward_code))) {
  warning("Some Ward codes are NA after join. Check ward_code_field.")
}


#### 10) Distribution comparison: District vs Ward ----------------------------
# For Figure 7: compare the dispersion of standardized rainfall shocks at the
# district vs ward level.
#
# Steps:
# (1) Compute a national monthly benchmark (mean and sd) using district rainfall.
# (2) Standardize rainfall for each unit-month: (rain_mm - nat_mean) / nat_sd.
# (3) Pool district and ward shocks and plot their distributions.
#
# Note: The benchmark is defined at the national-by-month level to keep the
# comparison on a common scale across spatial units.

# National monthly mean & sd (benchmark)
national_month_mean <- rain_dc_monthly |>
  group_by(year, month) |>
  summarise(
    nat_mean = mean(rain_mm, na.rm = TRUE),
    nat_sd   = sd(rain_mm, na.rm = TRUE),
    .groups  = "drop"
  )

# District-level standardized shock
dc_shock <- rain_dc_monthly |>
  left_join(national_month_mean, by = c("year", "month")) |>
  mutate(
    shock = (rain_mm - nat_mean) / nat_sd
  ) |>
  transmute(shock, level = "District")

# Ward-level standardized shock
ward_shock <- rain_ward_monthly |>
  left_join(national_month_mean, by = c("year", "month")) |>
  mutate(
    shock = (rain_mm - nat_mean) / nat_sd
  ) |>
  transmute(shock, level = "Ward")

# Figure7
plot_shock <- bind_rows(dc_shock, ward_shock) |>
  filter(is.finite(shock))

ggplot(plot_shock, aes(x = shock, color = level)) +
  geom_density()


