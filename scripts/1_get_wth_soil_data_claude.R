##### 1_get_wth_soil_data.R - DNDC VERSION #####
# NASA POWER-based implementation following ORYZA pattern
# Downloads climate + soil data for Myanmar DNDC simulations

suppressPackageStartupMessages({
  library(nasapower)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(geodata)
})

# --- PATHS ---
path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"

dir.create(path, FALSE, TRUE)
setwd(path)
raw_dir <- file.path("data", "raw")
power_dir <- file.path(raw_dir, "weather", "power")
dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)

# --- CONFIG ---
# DNDC requires: T2M_MAX, T2M_MIN, T2M (avg), PRECTOTCORR, WS2M, ALLSKY_SFC_SW_DWN, RH2M
vars <- c("T2M_MAX", "T2M_MIN", "T2M", "PRECTOTCORR", "WS2M", "ALLSKY_SFC_SW_DWN", "RH2M")
years <- 2005:2024
ext <- terra::ext(91.5, 101.5, 8, 29)  # xmin, xmax, ymin, ymax (Myanmar)
bbox <- c(xmin(ext), ymin(ext), xmax(ext), ymax(ext))  # lon_min, lat_min, lon_max, lat_max for nasapower

# Build the target grid (POWER native grid is 0.5-degree)
lon_seq <- seq(xmin(ext), xmax(ext), by = 0.5)
lat_seq <- seq(ymin(ext), ymax(ext), by = 0.5)
lat_seq_desc <- sort(lat_seq, decreasing = TRUE)

# Helper: turn a POWER dataframe for a single day into a raster layer
.make_layer <- function(df, value_col){
  r <- terra::rast(ncols = length(lon_seq), nrows = length(lat_seq_desc),
                   xmin = xmin(ext), xmax = xmax(ext), ymin = ymin(ext), ymax = ymax(ext),
                   crs = "EPSG:4326")
  df_ord <- df |>
    dplyr::arrange(dplyr::desc(LAT), LON)
  terra::values(r) <- df_ord[[value_col]]
  r
}

# Helper: fetch a single year's worth for one parameter using tiles
# NASA POWER limits: 
#   - Max 100 points per request (10° x 10° at 0.5° resolution)
#   - Minimum 2° x 2° per request
.fetch_year <- function(parm, yr){
  d1 <- paste0(yr, "-01-01"); d2 <- paste0(yr, "-12-31")
  
  # Use 4° x 4° tiles (8x8 points = 64 points, safe margin)
  tile_size <- 4.0
  min_size <- 2.0   # NASA POWER minimum
  
  # Generate tile grid
  lon_min <- xmin(ext)
  lon_max <- xmax(ext)
  lat_min <- ymin(ext)
  lat_max <- ymax(ext)
  
  # Create tile boundaries ensuring minimum 2° coverage
  lon_starts <- seq(lon_min, lon_max, by = tile_size)
  lat_starts <- seq(lat_min, lat_max, by = tile_size)
  
  # Adjust last tiles if they would be too small
  # If remaining space < 2°, extend the second-to-last tile
  if (length(lon_starts) > 1) {
    last_lon_extent <- lon_max - lon_starts[length(lon_starts)]
    if (last_lon_extent < min_size && last_lon_extent > 0) {
      # Remove last start, extend previous tile
      lon_starts <- lon_starts[-length(lon_starts)]
    }
  }
  
  if (length(lat_starts) > 1) {
    last_lat_extent <- lat_max - lat_starts[length(lat_starts)]
    if (last_lat_extent < min_size && last_lat_extent > 0) {
      # Remove last start, extend previous tile
      lat_starts <- lat_starts[-length(lat_starts)]
    }
  }
  
  all_tiles_data <- list()
  tile_count <- 0
  total_tiles <- length(lon_starts) * length(lat_starts)
  
  for (lon_start in lon_starts) {
    for (lat_start in lat_starts) {
      tile_count <- tile_count + 1
      
      # Calculate tile end
      lon_end <- min(lon_start + tile_size, lon_max)
      lat_end <- min(lat_start + tile_size, lat_max)
      
      # Verify minimum size
      tile_lon_size <- lon_end - lon_start
      tile_lat_size <- lat_end - lat_start
      
      message(sprintf("    Tile %d/%d [%.1f-%.1f, %.1f-%.1f] (%.1f° x %.1f°)", 
                      tile_count, total_tiles,
                      lon_start, lon_end, lat_start, lat_end,
                      tile_lon_size, tile_lat_size))
      
      # Double check minimum size
      if (tile_lon_size < min_size || tile_lat_size < min_size) {
        message(sprintf("      WARNING: Tile smaller than 2° minimum - skipping"))
        next
      }
      
      tile_bbox <- c(lon_start, lat_start, lon_end, lat_end)
      
      # Download with retry logic
      max_retries <- 3
      success <- FALSE
      
      for (attempt in 1:max_retries) {
        tryCatch({
          dat_tile <- nasapower::get_power(
            community = "ag",
            lonlat = tile_bbox,
            pars = parm,
            dates = c(d1, d2),
            temporal_api = "daily"
          )
          all_tiles_data[[length(all_tiles_data) + 1]] <- dat_tile
          success <- TRUE
          
          # Add delay to avoid rate limiting
          Sys.sleep(2)
          break
          
        }, error = function(e) {
          err_msg <- as.character(e)
          
          if (grepl("429", err_msg) || grepl("rate limit", err_msg, ignore.case = TRUE)) {
            wait_time <- 60 * attempt
            message(sprintf("      Rate limited. Waiting %d seconds... (attempt %d/%d)", 
                            wait_time, attempt, max_retries))
            Sys.sleep(wait_time)
            if (attempt == max_retries) {
              stop("Max retries exceeded due to rate limiting. Wait 1 hour and try again.")
            }
          } else {
            message(sprintf("      ERROR: %s", err_msg))
            if (attempt == max_retries) {
              stop(paste("Failed to download tile after", max_retries, "attempts:", err_msg))
            }
            Sys.sleep(10)
          }
        })
      }
      
      if (!success) {
        warning(sprintf("Tile %d failed to download", tile_count))
      }
    }
  }
  
  # Combine all tiles
  if (length(all_tiles_data) == 0) {
    stop("No tiles were successfully downloaded")
  }
  
  message(sprintf("    Successfully downloaded %d/%d tiles", 
                  length(all_tiles_data), total_tiles))
  
  dat <- dplyr::bind_rows(all_tiles_data)
  dat <- dplyr::rename(dat, DATE = YYYYMMDD)
  dat$DATE <- as.Date(dat$DATE)
  
  # Remove duplicates (from overlapping tiles)
  dat <- dat[!duplicated(dat[, c("LAT", "LON", "DATE")]), ]
  
  # Create raster layers
  days <- sort(unique(dat$DATE))
  layers <- vector("list", length(days))
  for (i in seq_along(days)){
    dd <- dat[dat$DATE == days[i], c("LAT", "LON", parm)]
    layers[[i]] <- .make_layer(dd, parm)
  }
  
  r <- terra::rast(layers)
  terra::time(r) <- days
  names(r) <- paste0(parm, "_", format(days, "%Y%j"))
  r
}

# Helper: load if exists; otherwise fetch all years and write a single NetCDF per variable
.get_or_load_var <- function(parm){
  fn_raw <- file.path(power_dir, sprintf("%s-2005_2024-91.5x101.5x8x29.nc", parm))
  if (file.exists(fn_raw)){
    message("Found existing ", basename(fn_raw), "; loading…")
    return(terra::rast(fn_raw))
  }
  
  message("Downloading ", parm, " (this will take 10-20 minutes per year with rate limiting)")
  message("Progress: 0/", length(years), " years")
  
  rs <- list()
  for (i in seq_along(years)) {
    yy <- years[i]
    message(sprintf("  Year %d (%d/%d)", yy, i, length(years)))
    rs[[i]] <- .fetch_year(parm, yy)
  }
  
  message("  Combining years and writing NetCDF...")
  R <- terra::rast(rs)
  terra::writeCDF(R, fn_raw, varname = parm, longname = parm, overwrite = TRUE)
  message("  ✓ Complete: ", basename(fn_raw))
  R
}

# -----------------------------------------------------------------------------
# 1) Download (or load existing) raw POWER variables
# -----------------------------------------------------------------------------
message("=== Downloading DNDC Climate Variables ===")
message("NOTE: NASA POWER API has rate limits and 100-point request limits")
message("      Myanmar grid will be downloaded in tiles with 2-second delays")
message("      Expected time: 3-6 hours for all variables (2005-2024)")
message("      Downloads are resume-safe - rerun if interrupted")
message("")
R_TMAX     <- .get_or_load_var("T2M_MAX")            # °C
R_TMIN     <- .get_or_load_var("T2M_MIN")            # °C
R_TAVG     <- .get_or_load_var("T2M")                # °C (daily average)
R_PREC     <- .get_or_load_var("PRECTOTCORR")        # mm day-1
R_WIND     <- .get_or_load_var("WS2M")               # m s-1
R_ALLSKY   <- .get_or_load_var("ALLSKY_SFC_SW_DWN")  # kWh m-2 day-1
R_RH       <- .get_or_load_var("RH2M")               # % relative humidity

# -----------------------------------------------------------------------------
# 2) Create DNDC-unit derived products
# -----------------------------------------------------------------------------
message("=== Converting to DNDC Units ===")

# Tmax: already in °C, just copy
fn_tmax <- file.path(power_dir, "tmax-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_tmax)){
  terra::writeCDF(R_TMAX, fn_tmax, varname = "tmax", longname = "Daily Max Temperature", unit = "C", overwrite = TRUE)
} else {
  R_TMAX <- terra::rast(fn_tmax)
}

# Tmin: already in °C, just copy
fn_tmin <- file.path(power_dir, "tmin-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_tmin)){
  terra::writeCDF(R_TMIN, fn_tmin, varname = "tmin", longname = "Daily Min Temperature", unit = "C", overwrite = TRUE)
} else {
  R_TMIN <- terra::rast(fn_tmin)
}

# Tavg: already in °C, just copy
fn_tavg <- file.path(power_dir, "tavg-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_tavg)){
  terra::writeCDF(R_TAVG, fn_tavg, varname = "tavg", longname = "Daily Avg Temperature", unit = "C", overwrite = TRUE)
} else {
  R_TAVG <- terra::rast(fn_tavg)
}

# Precipitation: mm to cm (divide by 10)
fn_prec <- file.path(power_dir, "prec-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_prec)){
  prec <- R_PREC / 10  # mm -> cm
  terra::writeCDF(prec, fn_prec, varname = "prec", longname = "Daily Precipitation", unit = "cm", overwrite = TRUE)
} else {
  prec <- terra::rast(fn_prec)
}

# Wind: already in m/s, just copy
fn_wind <- file.path(power_dir, "wind-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_wind)){
  terra::writeCDF(R_WIND, fn_wind, varname = "wind", longname = "Wind Speed", unit = "m/s", overwrite = TRUE)
} else {
  R_WIND <- terra::rast(fn_wind)
}

# Solar radiation: kWh m-2 d-1 to MJ m-2 d-1 (multiply by 3.6)
fn_srad <- file.path(power_dir, "srad-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_srad)){
  srad <- R_ALLSKY * 3.6  # kWh m-2 d-1 -> MJ m-2 d-1
  terra::writeCDF(srad, fn_srad, varname = "srad", longname = "Solar Radiation", unit = "MJ/m2/day", overwrite = TRUE)
} else {
  srad <- terra::rast(fn_srad)
}

# Relative humidity: already in %, just copy
fn_rhum <- file.path(power_dir, "rhum-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_rhum)){
  terra::writeCDF(R_RH, fn_rhum, varname = "rhum", longname = "Relative Humidity", unit = "%", overwrite = TRUE)
} else {
  R_RH <- terra::rast(fn_rhum)
}

message("✓ Climate data ready in DNDC units")

# -----------------------------------------------------------------------------
# 3) Elevation (for reference, using same approach as ORYZA)
# -----------------------------------------------------------------------------
message("=== Processing Elevation ===")
elv_tif <- file.path(raw_dir, "elevation.tif")
if (!file.exists(elv_tif)){
  elv <- geodata::elevation_30s("Myanmar", path = raw_dir)
  elv <- terra::resample(elv, R_TMAX[[1]], "average", filename = elv_tif)
} else {
  elv <- terra::rast(elv_tif)
}
message("✓ Elevation ready")

# -----------------------------------------------------------------------------
# 4) Soil (using geodata like ORYZA, but DNDC variables)
# -----------------------------------------------------------------------------
message("=== Processing Soil Data ===")
soil_tif <- file.path(raw_dir, "soil.tif")
soil_agg_tif <- file.path(raw_dir, "soil_agg.tif")

if (!file.exists(soil_agg_tif)){
  if (!file.exists(soil_tif)){
    # DNDC needs: bdod, clay, sand, silt, soc, phh2o
    vars_soil <- c("bdod", "clay", "sand", "silt", "soc", "phh2o")
    depths <- c(5, 15, 30)  # 0-5cm, 5-15cm, 15-30cm
    soil <- geodata::soil_world(vars_soil, depths, stat = "mean", vsi = TRUE)
    soil <- terra::crop(soil, ext, filename = soil_tif)
  } else {
    soil <- terra::rast(soil_tif)
  }
  
  # Convert SoilGrids units to DNDC units
  # Bulk density: cg/cm³ to g/cm³ (divide by 100)
  soil[["bdod_0-5cm"]] <- soil[["bdod_0-5cm"]] / 100
  soil[["bdod_5-15cm"]] <- soil[["bdod_5-15cm"]] / 100
  soil[["bdod_15-30cm"]] <- soil[["bdod_15-30cm"]] / 100
  
  # Clay, sand, silt: g/kg to fraction (divide by 1000)
  soil[["clay_0-5cm"]] <- soil[["clay_0-5cm"]] / 1000
  soil[["sand_0-5cm"]] <- soil[["sand_0-5cm"]] / 1000
  soil[["silt_0-5cm"]] <- soil[["silt_0-5cm"]] / 1000
  
  # SOC: dg/kg to fraction (divide by 1000)
  soil[["soc_0-5cm"]] <- soil[["soc_0-5cm"]] / 1000
  
  # pH: pH*10 to pH (divide by 10)
  soil[["phh2o_0-5cm"]] <- soil[["phh2o_0-5cm"]] / 10
  
  # Calculate porosity from bulk density (using particle density = 2.65 g/cm³)
  soil[["poros_0-5cm"]] <- 1 - (soil[["bdod_0-5cm"]] / 2.65)
  
  # Calculate texture class from sand/silt/clay
  clay_pct <- soil[["clay_0-5cm"]] * 100
  sand_pct <- soil[["sand_0-5cm"]] * 100
  silt_pct <- soil[["silt_0-5cm"]] * 100
  
  # USDA texture classification (simplified for DNDC)
  texture <- clay_pct
  texture[] <- 4  # Default to loam
  
  # Apply texture rules
  texture[clay_pct > 40 & sand_pct > 45] <- 10  # Sandy clay
  texture[clay_pct > 40 & silt_pct > 40] <- 11  # Silty clay
  texture[clay_pct > 40 & sand_pct <= 45 & silt_pct <= 40] <- 12  # Clay
  texture[clay_pct >= 27 & clay_pct <= 40 & sand_pct > 20 & sand_pct < 45] <- 7  # Sandy clay loam
  texture[clay_pct >= 27 & clay_pct <= 40 & silt_pct > 40] <- 9  # Silty clay loam
  texture[clay_pct >= 27 & clay_pct <= 40 & sand_pct >= 20 & sand_pct <= 45 & silt_pct <= 40] <- 8  # Clay loam
  texture[clay_pct < 27 & clay_pct >= 7 & sand_pct < 50 & silt_pct >= 28] <- 5  # Silt loam
  texture[clay_pct < 27 & clay_pct >= 7 & silt_pct >= 80] <- 6  # Silt
  texture[clay_pct < 27 & clay_pct >= 7 & sand_pct >= 23 & sand_pct < 52 & silt_pct < 50] <- 4  # Loam
  texture[clay_pct < 27 & clay_pct >= 7 & sand_pct >= 50] <- 3  # Sandy loam
  texture[clay_pct < 7 & sand_pct >= 85] <- 1  # Sand
  texture[clay_pct < 7 & sand_pct >= 70 & sand_pct < 85] <- 2  # Loamy sand
  
  names(texture) <- "texture_0-5cm"
  soil <- c(soil, texture)
  
  # Resample to climate grid
  soil2 <- terra::resample(soil, elv, "average", filename = soil_agg_tif, overwrite = TRUE)
} else {
  soil2 <- terra::rast(soil_agg_tif)
}
message("✓ Soil data ready in DNDC units")

# -----------------------------------------------------------------------------
# 5) Rice mask and cells table
# -----------------------------------------------------------------------------
message("=== Creating Rice Cells Table ===")
cells_rds <- file.path("data", "cells.rds")

if (!file.exists(cells_rds)){
  # Load Myanmar admin boundary
  aoi <- geodata::gadm("Myanmar", level = 1, path = raw_dir)
  
  # Create base grid from elevation
  r <- rast(elv)
  r <- mask(init(r, "cell"), aoi, touches = TRUE)
  
  # Load SPAM rice area
  rice_rast <- geodata::crop_spam(crop = "rice", var = "area", raw_dir) |> 
    terra::crop(aoi, mask = TRUE)
  rice_rast <- resample(rice_rast[[1]], r, "sum") > 0
  r <- mask(r, rice_rast, maskvalue = FALSE)
  cells <- data.frame(r)[,1]
  
  # Get xy coordinates and elevation
  xy <- data.frame(xyFromCell(r, cells))
  relv <- rast(elv_tif)
  xy$elevation <- round(relv[cells])
  xy$cell <- cells
  
  # Extract soil properties for each cell
  soil_data <- terra::extract(soil2, cells)
  xy <- cbind(xy, soil_data[, -1])  # Remove ID column
  
  saveRDS(xy, cells_rds)
  message(sprintf("✓ Created cells table with %d rice cells", nrow(xy)))
} else {
  xy <- readRDS(cells_rds)
  message(sprintf("✓ Loaded existing cells table with %d rice cells", nrow(xy)))
}

# Print summary
message("\n=== SUMMARY ===")
message(sprintf("Climate data: 2005-2024 (%d years)", length(years)))
message(sprintf("Variables: %s", paste(c("tmax", "tmin", "tavg", "prec", "wind", "srad", "rhum"), collapse = ", ")))
message(sprintf("Grid resolution: 0.5° x 0.5°"))
message(sprintf("Rice cells: %d", nrow(xy)))
message(sprintf("Cells table saved: %s", cells_rds))
message("\n✓ Step 1 complete. Ready for Step 2 (make_climate_files.R)")

##### END OF: 1_get_wth_soil_data.R #####