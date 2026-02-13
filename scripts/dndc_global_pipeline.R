# DNDC Global Prep Pipeline (0.25°) — CHIRPS + SPAM + Stubs for SoilGrids & N-dep
# Author: KGML/DNDC pipeline
# Target directory: G:/My Drive/Research/simulation/main_dndc/simulate_dndc (from all_file_paths_main_dndc.txt)
# Notes:
#  - Avoids APIs; uses direct HTTP downloads where possible.
#  - CHIRPS v2.0 daily precipitation (0.25°) is downloaded per year (1981–present).
#  - SPAM 2020 rice masks (irrigated & rainfed) are expected to exist locally.
#  - SoilGrids & N-deposition downloads are provided as HTTP stubs; adjust as needed.
#  - Produces global 0.25° grid with GridID, CHIRPS stacking sanity files,
#    and rice area aggregated to grid cells.
#  - Extend with temperature/wind/radiation/humidity downloads (non-API) if you have sources.
#  - After climate+soil are ready, map to DNDC GIS & climate library format per the manual.
#
# HOW TO RUN
#   Sys.setenv(DNDC_SIM_DIR = "G:/My Drive/Research/simulation/main_dndc/simulate_dndc")  # optional override
#   source("dndc_global_pipeline.R")
#   # Check console messages for paths of outputs.
#
suppressPackageStartupMessages({
  pkgs <- c("terra","sf","data.table","dplyr","arrow","glue","readr","tibble","tidyr")
  miss <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
  library(terra); library(sf); library(data.table); library(dplyr)
  library(arrow); library(glue); library(readr); library(tibble); library(tidyr)
})

options(timeout = max(3600, getOption("timeout", 60)))
terraOptions(progress = 1, memfrac = 0.7)

# ------------------------------------------------------------------
# 0) Directories & basic settings (paths from all_file_paths_main_dndc.txt)
# ------------------------------------------------------------------
proj_dir <- Sys.getenv("DNDC_SIM_DIR", unset = "G:/My Drive/Research/simulation/main_dndc/simulate_dndc")
raw_dir  <- file.path(proj_dir, "data", "raw")
proc_dir <- file.path(proj_dir, "data", "processed")
grid_dir <- file.path(proj_dir, "grids")
look_dir <- file.path(proj_dir, "lookups")
clim_dir <- file.path(proj_dir, "Lib_clim")  # DNDC-style climate library root
gis_dir  <- file.path(proj_dir, "GIS")       # DNDC GIS text files
for (d in c(proj_dir, raw_dir, proc_dir, grid_dir, look_dir, clim_dir, gis_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# SPAM rice masks (under simulate_dndc/rice mask/rice_mask_global)
spam_rice_i <- file.path(proj_dir, "rice mask", "rice_mask_global", "spam2020_V2r0_global_H_RICE_I.tif")
spam_rice_a <- file.path(proj_dir, "rice mask", "rice_mask_global", "spam2020_V2r0_global_H_RICE_A.tif")

# Global 0.25° grid (WGS84)
res_deg    <- 0.25
crs_epsg   <- "EPSG:4326"
world_bbox <- c(-180, -90, 180, 90)  # xmin, ymin, xmax, ymax

# CHIRPS v2.0 (0.25°) daily precipitation
chirps_res   <- "p25"
chirps_years <- 1981:2024
chirps_base  <- glue("https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/{chirps_res}")
chirps_dir   <- file.path(raw_dir, glue("chirps_{chirps_res}")); dir.create(chirps_dir, recursive = TRUE, showWarnings = FALSE)

# SoilGrids stubs — 0–5 cm mean layers (COG VRTs). Adjust/extend as needed.
# (Large downloads; consider regional clipping or using tiles.)
soil_dir <- file.path(raw_dir, "SoilGrids"); dir.create(soil_dir, recursive = TRUE, showWarnings = FALSE)
soil_urls <- c(
  soc  = "https://files.isric.org/soilgrids/latest/data/soc/soc_0-5cm_mean.vrt",
  clay = "https://files.isric.org/soilgrids/latest/data/clay/clay_0-5cm_mean.vrt",
  ph   = "https://files.isric.org/soilgrids/latest/data/phh2o/phh2o_0-5cm_mean.vrt",
  bd   = "https://files.isric.org/soilgrids/latest/data/bdod/bdod_0-5cm_mean.vrt"
)
# N-deposition stub (example; replace with your preferred product/URL):
n_dep_dir <- file.path(raw_dir, "N_dep"); dir.create(n_dep_dir, recursive = TRUE, showWarnings = FALSE)
n_dep_url <- "https://example.org/path/to/n_dep_kgN_ha_yr_global.nc"   # TODO: replace
n_dep_nc  <- file.path(n_dep_dir, "n_dep_kgN_ha_yr_global.nc")

# ------------------------------------------------------------------
# 1) Build global 0.25° grid and save
# ------------------------------------------------------------------
message("[1/8] Building global 0.25° grid…")
extent_obj   <- ext(world_bbox[1], world_bbox[3], world_bbox[2], world_bbox[4])
grid_template<- rast(ext = extent_obj, resolution = res_deg, crs = crs_epsg)

grid_sv <- as.polygons(grid_template, dissolve = FALSE)
grid_sf <- st_as_sf(grid_sv); rm(grid_sv)
grid_sf$grid_id <- seq_len(nrow(grid_sf))

centers <- st_centroid(grid_sf)
coords  <- st_coordinates(centers)
grid_sf$lon <- coords[,1]; grid_sf$lat <- coords[,2]
rm(centers, coords)

gpkg_path   <- file.path(grid_dir, sprintf("global_grid_%.2f.gpkg", res_deg))
if (file.exists(gpkg_path)) unlink(gpkg_path)
st_write(grid_sf, gpkg_path, layer = "grid", quiet = TRUE)
arrow::write_parquet(st_drop_geometry(grid_sf), file.path(grid_dir, sprintf("global_grid_%.2f.parquet", res_deg)))

message("  → Grid saved: ", gpkg_path)

# ------------------------------------------------------------------
# 2) Download CHIRPS daily NetCDFs (per-year) & stack
# ------------------------------------------------------------------
message("[2/8] Downloading CHIRPS v2.0 daily (", chirps_res, ") per-year…")
chirps_url_y <- function(y) glue("{chirps_base}/chirps-v2.0.{y}.days_{chirps_res}.nc")
chirps_file_y<- function(y) file.path(chirps_dir, glue("chirps-v2.0.{y}.days_{chirps_res}.nc"))

for (y in chirps_years) {
  url <- chirps_url_y(y); dst <- chirps_file_y(y)
  if (file.exists(dst) && file.size(dst) > 1e6) { message(glue("  ✓ {basename(dst)}")); next }
  message(glue("  ⇣ {basename(dst)}"))
  try(utils::download.file(url, dst, mode = "wb", quiet = FALSE), silent = TRUE)
}

message("[3/8] Loading CHIRPS stack…")
nc_files <- list.files(chirps_dir, pattern = "\\.nc$", full.names = TRUE)
if (!length(nc_files)) stop("No CHIRPS NetCDF files found in ", chirps_dir)
chirps_r <- try(rast(nc_files), silent = TRUE)
if (inherits(chirps_r, "try-error")) stop("terra::rast failed to open CHIRPS NetCDF files.")

chirps_r <- crop(chirps_r, extent_obj)
quick_tif <- file.path(proc_dir, "chirps_quicklook_layer1.tif")
writeRaster(chirps_r[[1]], quick_tif, overwrite = TRUE)

layer_time <- try(time(chirps_r), silent = TRUE)
idx_df <- tibble(layer = seq_len(nlyr(chirps_r)),
                 date  = if (inherits(layer_time,"Date")||inherits(layer_time,"POSIXct")) layer_time else NA)
readr::write_csv(idx_df, file.path(proc_dir, "chirps_layer_index.csv"))

message("  → CHIRPS sanity: ", quick_tif)

# ------------------------------------------------------------------
# 3) SPAM rice masks → rice area per grid cell
# ------------------------------------------------------------------
message("[4/8] Aggregating SPAM rice masks to 0.25° grid…")
if (!file.exists(spam_rice_i) || !file.exists(spam_rice_a)) {
  warning("SPAM rice masks not found. Please place the TIFFs at:",
          "\n  ", spam_rice_i, "\n  ", spam_rice_a)
} else {
  r_i <- rast(spam_rice_i)
  r_a <- rast(spam_rice_a)
  # Project to grid CRS if needed
  if (!identical(crs(r_i), crs(grid_template))) r_i <- project(r_i, crs(grid_template))
  if (!identical(crs(r_a), crs(grid_template))) r_a <- project(r_a, crs(grid_template))
  # Aggregate to 0.25° by sum of hectares
  fact_i <- round(res(grid_template)[1] / res(r_i)[1])
  fact_a <- round(res(grid_template)[1] / res(r_a)[1])
  r_i_agg <- if (fact_i > 1) aggregate(r_i, fact = fact_i, fun = sum, na.rm = TRUE) else r_i
  r_a_agg <- if (fact_a > 1) aggregate(r_a, fact = fact_a, fun = sum, na.rm = TRUE) else r_a
  # Extract to grid polygons (mean is fine if units are ha/pixel and pixels align; else exactextractr)
  vals_i <- terra::extract(r_i_agg, vect(grid_sf), fun = sum, na.rm = TRUE)
  vals_a <- terra::extract(r_a_agg, vect(grid_sf), fun = sum, na.rm = TRUE)
  rice_dt <- data.table(grid_id = grid_sf$grid_id,
                        rice_irr_ha = vals_i[,2],
                        rice_rf_ha  = vals_a[,2])
  rice_dt[, rice_total_ha := fifelse(is.na(rice_irr_ha),0,rice_irr_ha) + fifelse(is.na(rice_rf_ha),0,rice_rf_ha)]
  arrow::write_parquet(rice_dt, file.path(look_dir, "rice_area_grid_025.parquet"))
  message("  → Rice lookup: ", file.path(look_dir, "rice_area_grid_025.parquet"))
}

# ------------------------------------------------------------------
# 4) (Stub) SoilGrids download & extraction (0–5 cm)
# ------------------------------------------------------------------
message("[5/8] SoilGrids stubs — download VRTs (optional large)…")
for (nm in names(soil_urls)) {
  url <- soil_urls[[nm]]
  dst <- file.path(soil_dir, basename(url))
  if (file.exists(dst) && file.size(dst) > 1000) { message("  ✓ ", basename(dst)); next }
  message("  ⇣ ", basename(dst))
  try(utils::download.file(url, dst, mode = "wb", quiet = FALSE), silent = TRUE)
}

soil_lookup_path <- file.path(look_dir, "soil_top_0_5cm_lookup.parquet")
if (all(file.exists(file.path(soil_dir, basename(soil_urls))))) {
  message("  → Extracting 0–5 cm mean values to grid…")
  soil_dt <- data.table(grid_id = grid_sf$grid_id, lon = grid_sf$lon, lat = grid_sf$lat)
  # Note: VRTs reference COG mosaics; terra can read directly if GDAL supports it.
  for (nm in names(soil_urls)) {
    vpath <- file.path(soil_dir, basename(soil_urls[[nm]]))
    r <- try(rast(vpath), silent = TRUE)
    if (inherits(r,"try-error")) { warning("Failed to open ", vpath); next }
    r <- crop(r, extent_obj)
    vals <- terra::extract(r, vect(grid_sf), fun = mean, na.rm = TRUE)
    soil_dt[[nm]] <- vals[,2]
  }
  # Unit adjustments (check actual units on your side):
  # soc (g/kg) → % ; clay (g/kg) → % ; ph (x10) → pH ; bdod (cg/cm3) → g/cm3
  if ("soc" %in% names(soil_dt))  soil_dt[, soc_pct  := soc/100]
  if ("clay"%in% names(soil_dt))  soil_dt[, clay_pct := clay/10]
  if ("ph"  %in% names(soil_dt))  soil_dt[, ph_val   := ph/10]
  if ("bd"  %in% names(soil_dt))  soil_dt[, bd_gcm3  := bd/100]
  arrow::write_parquet(soil_dt, soil_lookup_path)
  message("  → Soil lookup: ", soil_lookup_path)
} else {
  warning("SoilGrids VRTs not fully available. Skipping extraction.")
}

# ------------------------------------------------------------------
# 5) (Stub) N-deposition download & conversion to mg N/L
# ------------------------------------------------------------------
message("[6/8] N-deposition stub…")
if (!file.exists(n_dep_nc)) {
  message("  ⇣ Attempting to download N-dep NetCDF… (please replace URL if needed)")
  try(utils::download.file(n_dep_url, n_dep_nc, mode = "wb", quiet = FALSE), silent = TRUE)
}
n_dep_lookup_path <- file.path(look_dir, "n_dep_lookup.parquet")
if (file.exists(n_dep_nc)) {
  # Read N-dep (kg N/ha/yr) and compute mg/L using annual precip
  ndep_r <- try(rast(n_dep_nc), silent = TRUE)
  if (!inherits(ndep_r,"try-error")) {
    ndep_r <- crop(ndep_r, extent_obj)
    ndep_vals <- terra::extract(ndep_r, vect(grid_sf), fun = mean, na.rm = TRUE)
    n_dep_dt  <- data.table(grid_id = grid_sf$grid_id, n_kgN_ha_yr = ndep_vals[,2])
    # Compute precip_mm_year from CHIRPS dates (sum daily per year then average or choose a baseline year)
    # For now, compute a climatological mean precip using a sample year (e.g., 2000) if present:
    idx_csv <- file.path(proc_dir, "chirps_layer_index.csv")
    if (file.exists(idx_csv)) {
      idx_df <- readr::read_csv(idx_csv, show_col_types = FALSE)
      if (any(!is.na(idx_df$date))) {
        # pick a target year with full coverage (e.g., 2000) if present
        target_year <- 2000
        layers <- idx_df$layer[lubridate::year(idx_df$date) == target_year]
        if (length(layers)) {
          message("  → Summing CHIRPS for year ", target_year, " to compute precip_mm_year…")
          pr2000 <- sum(chirps_r[[layers]], na.rm = TRUE)
          pvals  <- terra::extract(pr2000, vect(grid_sf), fun = mean, na.rm = TRUE)
          precip_y <- data.table(grid_id = grid_sf$grid_id, precip_mm_year = pvals[,2])
          n_dep_dt <- merge(n_dep_dt, precip_y, by = "grid_id", all.x = TRUE)
          n_dep_dt[, N_dep_mgL := 100 * n_kgN_ha_yr / precip_mm_year]
        }
      }
    }
    arrow::write_parquet(n_dep_dt, n_dep_lookup_path)
    message("  → N-dep lookup: ", n_dep_lookup_path)
  } else {
    warning("Could not read N-dep NetCDF. Skipping conversion.")
  }
} else {
  warning("N-dep NetCDF not found; please provide a valid URL or local file.")
}

# ------------------------------------------------------------------
# 6) Climate library scaffolding (DNDC Station files; precip only demo)
# ------------------------------------------------------------------
message("[7/8] Climate library scaffolding (precip-only demo)…")
# Assign climate_file_id = grid_id for simplicity (or a separate mapping)
clim_index <- data.table(grid_id = grid_sf$grid_id,
                         climate_file_id = grid_sf$grid_id,
                         lon = grid_sf$lon, lat = grid_sf$lat)
arrow::write_parquet(clim_index, file.path(look_dir, "climate_index.parquet"))

# Build a small demo for a single year to show file format (expand later)
demo_year <- 2000
yr_idx <- readr::read_csv(file.path(proc_dir, "chirps_layer_index.csv"), show_col_types = FALSE)
yr_layers <- yr_idx$layer[lubridate::year(yr_idx$date) == demo_year]
if (length(yr_layers)) {
  dir.create(file.path(clim_dir, as.character(demo_year)), showWarnings = FALSE, recursive = TRUE)
  # Example: write for first 3 stations only (expand to all later)
  demo_ids <- head(clim_index$climate_file_id, 3)
  for (cid in demo_ids) {
    # subset grid cell polygon
    g <- grid_sf[grid_sf$grid_id == cid, ]
    # extract daily precip (mm); convert to cm if needed by DNDC
    v <- terra::extract(chirps_r[[yr_layers]], vect(g), fun = mean, na.rm = TRUE)
    # v is 1 x nDays; pivot to table
    p <- as.numeric(v[1, -1])
    dd <- data.table(
      doy = seq_along(p),
      tmin = NA_real_,
      tmax = NA_real_,
      precip = p / 10,  # mm → cm/day
      wind = NA_real_,
      rh   = NA_real_,
      rad  = NA_real_
    )
    # Write Station file
    f <- file.path(clim_dir, as.character(demo_year), sprintf("Station_%05d.txt", cid))
    con <- file(f, "w"); on.exit(close(con), add = TRUE)
    writeLines(sprintf("Station_%05d  %.4f  %.4f", cid, g$lat[1], g$lon[1]), con)
    write.table(dd, con, sep = " ", row.names = FALSE, col.names = FALSE)
    close(con)
  }
  message("  → Demo station files written in ", file.path(clim_dir, as.character(demo_year)))
} else {
  message("  • Skipping climate library demo (no dated layers found in CHIRPS index).")
}

# ------------------------------------------------------------------
# 7) DNDC ClimateSoil skeleton (requires soil & N-dep to be complete)
# ------------------------------------------------------------------
message("[8/8] ClimateSoil skeleton (fill once soil & N-dep complete)…")
clim_soil_path <- file.path(gis_dir, "Global_ClimateSoil.txt")
# Attempt to load soil & n-dep lookups if present
soil_path <- file.path(look_dir, "soil_top_0_5cm_lookup.parquet")
ndep_path <- file.path(look_dir, "n_dep_lookup.parquet")
if (file.exists(soil_path) && file.exists(ndep_path)) {
  soil_dt <- arrow::read_parquet(soil_path)
  ndep_dt <- arrow::read_parquet(ndep_path)
  cs <- merge(clim_index, soil_dt, by = c("grid_id","lon","lat"), all.x = TRUE)
  cs <- merge(cs, ndep_dt[, .(grid_id, N_dep_mgL)], by = "grid_id", all.x = TRUE)
  # Create ±10% ranges for required fields (SOC, clay, pH, bulk density)
  cs <- cs %>%
    mutate(
      SOCmax  = soc_pct  * 1.10, SOCmin  = soc_pct  * 0.90,
      CLAYmax = clay_pct * 1.10, CLAYmin = clay_pct * 0.90,
      PHmax   = ph_val   * 1.10, PHmin   = ph_val   * 0.90,
      DENSmax = bd_gcm3  * 1.10, DENSmin = bd_gcm3  * 0.90,
      slope   = 0, salinity = 0, column = 0, row = 0, region = 0
    ) %>%
    select(grid_id, lon, lat, climate_file_id, N_dep_mgL,
           SOCmax, SOCmin, CLAYmax, CLAYmin, PHmax, PHmin, DENSmax, DENSmin,
           slope, salinity, column, row, region)
  fwrite(cs, clim_soil_path, sep = " ", col.names = FALSE)
  message("  → ClimateSoil (skeleton) written: ", clim_soil_path)
} else {
  message("  • Soil or N-dep lookup missing; ClimateSoil not written. Re-run after completing those steps.")
}

message("\nAll done.\nKey outputs:\n  • Grid: ", gpkg_path,
        "\n  • CHIRPS downloads: ", chirps_dir,
        "\n  • CHIRPS quicklook: ", quick_tif,
        "\n  • Rice area lookup (if SPAM present): ", file.path(look_dir, "rice_area_grid_025.parquet"),
        "\n  • Soil/N-dep lookups (if downloaded): ", soil_lookup_path, " / ", n_dep_lookup_path,
        "\n  • Climate index parquet: ", file.path(look_dir, "climate_index.parquet"),
        "\n  • Demo climate library year: ", file.path(clim_dir, "2000"))
