##### 1_get_wth_soil_data.R - DNDC VERSION #####
# NASA POWER-based implementation following ORYZA pattern
# Downloads climate + soil data for Myanmar DNDC simulations

suppressPackageStartupMessages({
  library(nasapower)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(geodata)
  library(httr)  # For timeout configuration
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
  
  # Configure longer timeout for slow server responses (5 minutes)
  httr::set_config(httr::timeout(300))
  on.exit(httr::reset_config())  # Reset after function completes
  
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
      
      # Download with retry logic and increased timeout
      max_retries <- 5  # Increased from 3
      success <- FALSE
      
      for (attempt in 1:max_retries) {
        tryCatch({
          # Set timeout to 5 minutes (configured via httr above)
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
          Sys.sleep(3)  # Increased to 3 seconds
          break
          
        }, error = function(e) {
          err_msg <- as.character(e)
          
          if (grepl("Timeout", err_msg) || grepl("timed out", err_msg, ignore.case = TRUE)) {
            wait_time <- 30 * attempt  # 30s, 60s, 90s, 120s, 150s
            message(sprintf("      Timeout. Waiting %d seconds... (attempt %d/%d)", 
                            wait_time, attempt, max_retries))
            Sys.sleep(wait_time)
            if (attempt == max_retries) {
              message(sprintf("      Max retries exceeded. Skipping this tile."))
              message(sprintf("      This year may be incomplete. Rerun script to retry."))
              return(NULL)  # Skip this tile but continue
            }
          } else if (grepl("429", err_msg) || grepl("rate limit", err_msg, ignore.case = TRUE)) {
            wait_time <- 120 * attempt  # 2 min, 4 min, 6 min
            message(sprintf("      Rate limited. Waiting %d seconds... (attempt %d/%d)", 
                            wait_time, attempt, max_retries))
            Sys.sleep(wait_time)
            if (attempt == max_retries) {
              stop("Max retries exceeded due to rate limiting. Wait 1 hour and try again.")
            }
          } else {
            message(sprintf("      ERROR: %s", err_msg))
            if (attempt == max_retries) {
              message(sprintf("      Max retries exceeded. Skipping this tile."))
              return(NULL)  # Skip this tile but continue
            }
            Sys.sleep(20)
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

# Helper: load if exists AND is valid; otherwise fetch all years and write a single NetCDF per variable
.get_or_load_var <- function(parm){
  fn_raw <- file.path(power_dir, sprintf("%s-2005_2024-91.5x101.5x8x29.nc", parm))
  
  # Check if file exists and is valid
  if (file.exists(fn_raw)){
    message("Found existing ", basename(fn_raw), "; validating...")
    
    # Try to load - if it fails, file is corrupted
    valid <- tryCatch({
      test_load <- terra::rast(fn_raw)
      TRUE
    }, error = function(e) {
      message("  File is corrupted or unreadable. Will re-download.")
      FALSE
    })
    
    if (valid) {
      message("  File is valid; loading...")
      return(terra::rast(fn_raw))
    } else {
      # Delete corrupted file
      file.remove(fn_raw)
      message("  Deleted corrupted file")
    }
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
  message("✓ Tmax converted")
} else {
  message("✓ Tmax already exists")
}

# Tmin: already in °C, just copy
fn_tmin <- file.path(power_dir, "tmin-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_tmin)){
  terra::writeCDF(R_TMIN, fn_tmin, varname = "tmin", longname = "Daily Min Temperature", unit = "C", overwrite = TRUE)
  message("✓ Tmin converted")
} else {
  message("✓ Tmin already exists")
}

# Tavg: already in °C, just copy
fn_tavg <- file.path(power_dir, "tavg-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_tavg)){
  terra::writeCDF(R_TAVG, fn_tavg, varname = "tavg", longname = "Daily Avg Temperature", unit = "C", overwrite = TRUE)
  message("✓ Tavg converted")
} else {
  message("✓ Tavg already exists")
}

# Precipitation: mm to cm (divide by 10)
fn_prec <- file.path(power_dir, "prec-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_prec)){
  prec_file <- file.path(power_dir, "PRECTOTCORR-2005_2024-91.5x101.5x8x29.nc")
  
  if (!file.exists(prec_file)) {
    message("✗ PRECTOTCORR file missing")
    stop("Missing PRECTOTCORR file - rerun script")
  }
  
  R_PREC <- tryCatch({
    terra::rast(prec_file)
  }, error = function(e) {
    message("✗ PRECTOTCORR file is corrupted")
    file.remove(prec_file)
    stop("Corrupted PRECTOTCORR file - deleted. Please rerun script.")
  })
  
  prec <- R_PREC / 10  # mm -> cm
  terra::writeCDF(prec, fn_prec, varname = "prec", longname = "Daily Precipitation", unit = "cm", overwrite = TRUE)
  rm(R_PREC, prec)
  message("✓ Precipitation converted (mm → cm)")
} else {
  message("✓ Precipitation already exists")
}

# Wind: already in m/s, just copy
fn_wind <- file.path(power_dir, "wind-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_wind)){
  ws2m_file <- file.path(power_dir, "WS2M-2005_2024-91.5x101.5x8x29.nc")
  
  if (!file.exists(ws2m_file)) {
    message("✗ WS2M file missing - need to re-download")
    message("  Rerun script to download WS2M")
    stop("Missing WS2M file")
  }
  
  # Validate file before loading
  R_WIND <- tryCatch({
    terra::rast(ws2m_file)
  }, error = function(e) {
    message("✗ WS2M file is corrupted")
    message("  Deleting corrupted file...")
    file.remove(ws2m_file)
    message("  Rerun script to re-download WS2M")
    stop("Corrupted WS2M file - deleted. Please rerun script.")
  })
  
  terra::writeCDF(R_WIND, fn_wind, varname = "wind", longname = "Wind Speed", unit = "m/s", overwrite = TRUE)
  rm(R_WIND)
  message("✓ Wind converted")
} else {
  message("✓ Wind already exists")
}

# Solar radiation: kWh m-2 d-1 to MJ m-2 d-1 (multiply by 3.6)
fn_srad <- file.path(power_dir, "srad-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_srad)){
  allsky_file <- file.path(power_dir, "ALLSKY_SFC_SW_DWN-2005_2024-91.5x101.5x8x29.nc")
  
  if (!file.exists(allsky_file)) {
    message("✗ ALLSKY_SFC_SW_DWN file missing")
    stop("Missing ALLSKY file - rerun script")
  }
  
  R_ALLSKY <- tryCatch({
    terra::rast(allsky_file)
  }, error = function(e) {
    message("✗ ALLSKY_SFC_SW_DWN file is corrupted")
    file.remove(allsky_file)
    stop("Corrupted ALLSKY file - deleted. Please rerun script.")
  })
  
  srad <- R_ALLSKY * 3.6  # kWh m-2 d-1 -> MJ m-2 d-1
  terra::writeCDF(srad, fn_srad, varname = "srad", longname = "Solar Radiation", unit = "MJ/m2/day", overwrite = TRUE)
  rm(R_ALLSKY, srad)
  message("✓ Solar radiation converted (kWh → MJ)")
} else {
  message("✓ Solar radiation already exists")
}

# Relative humidity: already in %, just copy
fn_rhum <- file.path(power_dir, "rhum-2005_2024-91.5x101.5x8x29.nc")
if (!file.exists(fn_rhum)){
  rh2m_file <- file.path(power_dir, "RH2M-2005_2024-91.5x101.5x8x29.nc")
  
  if (!file.exists(rh2m_file)) {
    message("✗ RH2M file missing")
    stop("Missing RH2M file - rerun script")
  }
  
  R_RH <- tryCatch({
    terra::rast(rh2m_file)
  }, error = function(e) {
    message("✗ RH2M file is corrupted")
    file.remove(rh2m_file)
    stop("Corrupted RH2M file - deleted. Please rerun script.")
  })
  
  terra::writeCDF(R_RH, fn_rhum, varname = "rhum", longname = "Relative Humidity", unit = "%", overwrite = TRUE)
  rm(R_RH)
  message("✓ Relative humidity converted")
} else {
  message("✓ Relative humidity already exists")
}

message("✓ Climate data ready in DNDC units")

# -----------------------------------------------------------------------------
# 3) Elevation (for reference, using same approach as ORYZA)
# -----------------------------------------------------------------------------
message("=== Processing Elevation ===")
elv_tif <- file.path(raw_dir, "elevation.tif")
if (!file.exists(elv_tif)){
  # Load reference grid from saved tmax file
  ref_tmax <- terra::rast(file.path(power_dir, "T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))
  
  elv <- geodata::elevation_30s("Myanmar", path = raw_dir)
  elv <- terra::resample(elv, ref_tmax[[1]], "average", filename = elv_tif)
  rm(ref_tmax, elv)
  message("✓ Elevation processed")
} else {
  message("✓ Elevation already exists")
}
message("✓ Elevation ready")

# -----------------------------------------------------------------------------
# 4) Soil (using geodata like ORYZA, but DNDC variables)
# -----------------------------------------------------------------------------
message("=== Processing Soil Data ===")
soil_tif <- file.path(raw_dir, "soil.tif")
soil_agg_tif <- file.path(raw_dir, "soil_agg.tif")

if (!file.exists(soil_agg_tif)){
  # Load reference elevation for resampling
  elv_ref <- terra::rast(elv_tif)
  
  if (!file.exists(soil_tif)){
    # DNDC needs: bdod, clay, sand, silt, soc, phh2o
    vars_soil <- c("bdod", "clay", "sand", "silt", "soc", "phh2o")
    depths <- c(5, 15, 30)  # 0-5cm, 5-15cm, 15-30cm
    message("  Downloading soil data from SoilGrids...")
    soil <- geodata::soil_world(vars_soil, depths, stat = "mean", vsi = TRUE)
    
    # Debug: show what layer names we actually got
    message("  Available soil layers: ", paste(names(soil), collapse = ", "))
    
    soil <- terra::crop(soil, ext, filename = soil_tif)
  } else {
    message("  Loading existing soil data...")
    soil <- terra::rast(soil_tif)
    message("  Available soil layers: ", paste(names(soil), collapse = ", "))
    
    # Check if this is a partially processed file (has mixed layer names)
    layer_names <- names(soil)
    has_mean_suffix <- any(grepl("_mean$", layer_names))
    has_converted <- any(!grepl("_mean$", layer_names) & grepl("bdod|clay|sand|silt|soc|phh2o", layer_names))
    
    if (has_mean_suffix && has_converted) {
      message("  ⚠️  Detected partially processed soil file (mixed layer names)")
      message("  Deleting and re-downloading for clean data...")
      file.remove(soil_tif)
      
      # Re-download
      vars_soil <- c("bdod", "clay", "sand", "silt", "soc", "phh2o")
      depths <- c(5, 15, 30)
      soil <- geodata::soil_world(vars_soil, depths, stat = "mean", vsi = TRUE)
      soil <- terra::crop(soil, ext, filename = soil_tif)
      message("  ✓ Fresh soil data downloaded")
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Dynamic layer lookup - robust to _mean suffix or not
  # ---------------------------------------------------------------------------
  .get_soil_layer <- function(r, prefix) {
    n <- names(r)
    if (prefix %in% n)                          return(r[[prefix]])
    if (paste0(prefix, "_mean") %in% n)         return(r[[paste0(prefix, "_mean")]])
    matches <- grep(paste0("^", prefix), n, value = TRUE)
    if (length(matches) > 0)                    return(r[[matches[1]]])
    stop(sprintf("Cannot find layer '%s'. Available: %s", prefix, paste(n, collapse=", ")))
  }
  
  message("  Converting units...")
  
  # Bulk density: cg/cm3 to g/cm3 (divide by 100)
  soil[["bdod_0-5cm"]]   <- .get_soil_layer(soil, "bdod_0-5cm")   / 100
  soil[["bdod_5-15cm"]]  <- .get_soil_layer(soil, "bdod_5-15cm")  / 100
  soil[["bdod_15-30cm"]] <- .get_soil_layer(soil, "bdod_15-30cm") / 100
  
  # Clay, sand, silt: g/kg to fraction (divide by 1000)
  soil[["clay_0-5cm"]]   <- .get_soil_layer(soil, "clay_0-5cm")   / 1000
  soil[["sand_0-5cm"]]   <- .get_soil_layer(soil, "sand_0-5cm")   / 1000
  soil[["silt_0-5cm"]]   <- .get_soil_layer(soil, "silt_0-5cm")   / 1000
  
  # SOC: dg/kg to fraction (divide by 1000)
  soil[["soc_0-5cm"]]    <- .get_soil_layer(soil, "soc_0-5cm")    / 1000
  
  # pH: pH*10 to pH (divide by 10)
  soil[["pH_0-5cm"]]     <- .get_soil_layer(soil, "phh2o_0-5cm")  / 10
  
  message("  Units converted")
  
  # Porosity from bulk density (particle density = 2.65 g/cm3)
  soil[["poros_0-5cm"]] <- 1 - (soil[["bdod_0-5cm"]] / 2.65)
  
  # Texture classification (USDA) using values to avoid extent issues
  clay_vals <- terra::values(soil[["clay_0-5cm"]]) * 100
  sand_vals <- terra::values(soil[["sand_0-5cm"]]) * 100
  silt_vals <- terra::values(soil[["silt_0-5cm"]]) * 100
  
  texture_vals <- rep(4L, length(clay_vals))  # Default: loam
  texture_vals[clay_vals > 40 & sand_vals > 45]                                             <- 10L  # Sandy clay
  texture_vals[clay_vals > 40 & silt_vals > 40]                                             <- 11L  # Silty clay
  texture_vals[clay_vals > 40 & sand_vals <= 45 & silt_vals <= 40]                          <- 12L  # Clay
  texture_vals[clay_vals >= 27 & clay_vals <= 40 & sand_vals > 20 & sand_vals < 45]         <-  7L  # Sandy clay loam
  texture_vals[clay_vals >= 27 & clay_vals <= 40 & silt_vals > 40]                          <-  9L  # Silty clay loam
  texture_vals[clay_vals >= 27 & clay_vals <= 40 & sand_vals >= 20 & sand_vals <= 45 & silt_vals <= 40] <- 8L  # Clay loam
  texture_vals[clay_vals < 27 & clay_vals >= 7 & sand_vals < 50 & silt_vals >= 28]          <-  5L  # Silt loam
  texture_vals[clay_vals < 27 & clay_vals >= 7 & silt_vals >= 80]                           <-  6L  # Silt
  texture_vals[clay_vals < 27 & clay_vals >= 7 & sand_vals >= 23 & sand_vals < 52 & silt_vals < 50] <- 4L  # Loam
  texture_vals[clay_vals < 27 & clay_vals >= 7 & sand_vals >= 50]                           <-  3L  # Sandy loam
  texture_vals[clay_vals < 7  & sand_vals >= 85]                                            <-  1L  # Sand
  texture_vals[clay_vals < 7  & sand_vals >= 70 & sand_vals < 85]                           <-  2L  # Loamy sand
  
  texture <- soil[["clay_0-5cm"]]
  terra::values(texture) <- texture_vals
  names(texture) <- "texture_0-5cm"
  soil <- c(soil, texture)
  rm(clay_vals, sand_vals, silt_vals, texture_vals, texture)
  
  
  # Resample to climate grid (using elevation as reference)
  message("  Resampling soil to climate grid...")
  soil2 <- terra::resample(soil, elv_ref, "average", filename = soil_agg_tif, overwrite = TRUE)
  rm(elv_ref)
  message("✓ Soil data processed")
} else {
  message("✓ Soil data already exists")
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
  elv_for_cells <- terra::rast(elv_tif)
  r <- terra::rast(elv_for_cells)
  r <- mask(init(r, "cell"), aoi, touches = TRUE)
  
  # Load SPAM rice area
  rice_rast <- geodata::crop_spam(crop = "rice", var = "area", raw_dir) |> 
    terra::crop(aoi, mask = TRUE)
  rice_rast <- resample(rice_rast[[1]], r, "sum") > 0
  r <- mask(r, rice_rast, maskvalue = FALSE)
  cells <- data.frame(r)[,1]
  
  # Get xy coordinates and elevation
  xy <- data.frame(xyFromCell(r, cells))
  xy$elevation <- round(elv_for_cells[cells])
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