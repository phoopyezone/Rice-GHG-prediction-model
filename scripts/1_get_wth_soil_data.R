##### 1_get_wth_soil_data.R #####
#
# Download-only DNDC climate and soil inputs for Myanmar.
# Follows the Oryza downloader style (NASA POWER + geodata SoilGrids),
# but writes outputs using DNDC variable names and units.
#
suppressPackageStartupMessages({
  library(nasapower)
  library(terra)
  library(dplyr)
  library(sf)
  library(geodata)
})

# -----------------------------------------------------------------------------
# 0) Paths and setup
# -----------------------------------------------------------------------------
sim_root <- Sys.getenv(
  "DNDC_SIM_DIR",
  unset = "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
)

raw_weather_dir <- file.path(sim_root, "data", "raw", "weather", "power")
raw_soil_dir <- file.path(sim_root, "data", "raw", "soil")
dir.create(raw_weather_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_soil_dir, recursive = TRUE, showWarnings = FALSE)

# Myanmar boundary to derive download extent
mmr_shp <- file.path(
  sim_root,
  "boundaries_shapefiles",
  "mmr_polbnda2_adm1_250k_mimu_1.shp"
)
if (!file.exists(mmr_shp)) {
  stop("Myanmar boundary shapefile not found: ", mmr_shp)
}

mmr <- sf::st_read(mmr_shp, quiet = TRUE)
bb <- sf::st_bbox(mmr)

# Snap to POWER 0.5-degree grid (same pattern as Oryza workflow)
xmin_dl <- floor(as.numeric(bb["xmin"]) * 2) / 2
xmax_dl <- ceiling(as.numeric(bb["xmax"]) * 2) / 2
ymin_dl <- floor(as.numeric(bb["ymin"]) * 2) / 2
ymax_dl <- ceiling(as.numeric(bb["ymax"]) * 2) / 2

ext_dl <- terra::ext(xmin_dl, xmax_dl, ymin_dl, ymax_dl)
bbox_power <- c(ymin_dl, xmin_dl, ymax_dl, xmax_dl) # lat_min, lon_min, lat_max, lon_max

years <- 1995:2024
lon_seq <- seq(xmin_dl, xmax_dl, by = 0.5)
lat_seq <- seq(ymin_dl, ymax_dl, by = 0.5)
lat_seq_desc <- sort(lat_seq, decreasing = TRUE)

message("Download extent: ", paste(round(c(xmin_dl, xmax_dl, ymin_dl, ymax_dl), 2), collapse = ", "))
message("Simulation root: ", sim_root)

# -----------------------------------------------------------------------------
# 1) Climate download (NASA POWER) mapped to DNDC variables
# -----------------------------------------------------------------------------
# Mapping:
# - POWER variable      -> DNDC variable name, unit conversion to DNDC unit
climate_map <- data.frame(
  power_var = c("T2M_MAX", "T2M_MIN", "T2M", "PRECTOTCORR", "WS2M", "ALLSKY_SFC_SW_DWN", "RH2M"),
  dndc_var = c(
    "Daily_max_air_temperature",
    "Daily_min_air_temperature",
    "Daily_avg_air_temperature",
    "Daily_precipitation",
    "Daily_avg_wind_speed",
    "Daily_solar_radiation",
    "Daily_relative_humidity"
  ),
  conv_to_dndc = c(1, 1, 1, 0.1, 1, 3.6, 1),
  stringsAsFactors = FALSE
)

.make_layer <- function(df, value_col, ext_obj, lon_vals, lat_vals_desc) {
  r <- terra::rast(
    ncols = length(lon_vals), nrows = length(lat_vals_desc),
    xmin = terra::xmin(ext_obj), xmax = terra::xmax(ext_obj),
    ymin = terra::ymin(ext_obj), ymax = terra::ymax(ext_obj),
    crs = "EPSG:4326"
  )
  df_ord <- df |>
    dplyr::arrange(dplyr::desc(LAT), LON)
  terra::values(r) <- df_ord[[value_col]]
  r
}

.fetch_year <- function(parm, yr, bbox_vals, ext_obj, lon_vals, lat_vals_desc) {
  d1 <- paste0(yr, "-01-01")
  d2 <- paste0(yr, "-12-31")
  dat <- nasapower::get_power(
    community = "AG",
    temporal_average = "DAILY",
    pars = parm,
    bbox = bbox_vals,
    dates = c(d1, d2)
  )
  dat <- dplyr::rename(dat, DATE = YYYYMMDD)
  dat$DATE <- as.Date(dat$DATE)

  days <- sort(unique(dat$DATE))
  layers <- vector("list", length(days))
  for (i in seq_along(days)) {
    dd <- dat[dat$DATE == days[i], c("LAT", "LON", parm)]
    layers[[i]] <- .make_layer(dd, parm, ext_obj, lon_vals, lat_vals_desc)
  }

  r <- terra::rast(layers)
  terra::time(r) <- days
  names(r) <- paste0(parm, "_", format(days, "%Y%j"))
  r
}

.get_or_download_power <- function(parm) {
  fn_raw <- file.path(raw_weather_dir, sprintf("%s-1995_2024-myanmar.nc", parm))
  if (file.exists(fn_raw)) {
    message("Found existing POWER file: ", basename(fn_raw))
    return(terra::rast(fn_raw))
  }

  message("Downloading POWER variable: ", parm)
  rs <- lapply(years, function(yy) {
    message("  - Year ", yy)
    .fetch_year(parm, yy, bbox_power, ext_dl, lon_seq, lat_seq_desc)
  })
  R <- terra::rast(rs)
  terra::writeCDF(R, fn_raw, varname = parm, longname = parm, overwrite = TRUE)
  R
}

for (i in seq_len(nrow(climate_map))) {
  power_var <- climate_map$power_var[i]
  dndc_var <- climate_map$dndc_var[i]
  conv <- climate_map$conv_to_dndc[i]

  raw_r <- .get_or_download_power(power_var)

  fn_dndc <- file.path(raw_weather_dir, sprintf("%s-1995_2024-myanmar.nc", dndc_var))
  if (!file.exists(fn_dndc)) {
    message("Writing DNDC climate variable: ", dndc_var)
    dndc_r <- raw_r * conv
    terra::writeCDF(dndc_r, fn_dndc, varname = dndc_var, longname = dndc_var, overwrite = TRUE)
  } else {
    message("Found existing DNDC climate file: ", basename(fn_dndc))
  }
}

# -----------------------------------------------------------------------------
# 2) Soil download (SoilGrids via geodata) mapped to DNDC variables
# -----------------------------------------------------------------------------
# Download top-layer source variables needed to build selected DNDC inputs.
soil_src_vars <- c("bdod", "clay", "soc", "phh2o", "sand", "silt")
soil_depths <- c(5) # top layer only

soil_src_tif <- file.path(raw_soil_dir, "soilgrids_top_0_5cm_source.tif")
if (!file.exists(soil_src_tif)) {
  message("Downloading SoilGrids source layers (top 0-5 cm)...")
  soil_src <- geodata::soil_world(soil_src_vars, soil_depths, stat = "mean", path = raw_soil_dir)
  soil_src <- terra::crop(soil_src, ext_dl)
  terra::writeRaster(soil_src, soil_src_tif, overwrite = TRUE)
} else {
  message("Found existing SoilGrids source file: ", basename(soil_src_tif))
  soil_src <- terra::rast(soil_src_tif)
}

# Convert/write DNDC soil variables from downloaded source layers.
soil_outputs <- list(
  Bulk_density = list(src = "bdod_0-5cm", conv = 1 / 100),   # cg/cm3 -> g/cm3
  Clay_fraction = list(src = "clay_0-5cm", conv = 1 / 1000), # g/kg -> fraction
  Top_layer_SOC = list(src = "soc_0-5cm", conv = 1 / 1000),  # g/kg (scaled) -> fraction
  pH = list(src = "phh2o_0-5cm", conv = 1 / 10),             # pH*10 -> pH
  sand_pct = list(src = "sand_0-5cm", conv = 1 / 10),        # g/kg -> %
  silt_pct = list(src = "silt_0-5cm", conv = 1 / 10)         # g/kg -> %
)

for (nm in names(soil_outputs)) {
  src_name <- soil_outputs[[nm]]$src
  conv <- soil_outputs[[nm]]$conv

  if (!(src_name %in% names(soil_src))) {
    warning("Soil layer missing in downloaded stack: ", src_name)
    next
  }

  out_tif <- file.path(raw_soil_dir, paste0(nm, "_myanmar.tif"))
  if (!file.exists(out_tif)) {
    message("Writing DNDC soil variable: ", nm)
    r <- soil_src[[src_name]] * conv
    terra::writeRaster(r, out_tif, overwrite = TRUE)
  } else {
    message("Found existing DNDC soil file: ", basename(out_tif))
  }
}

message("Download step complete.")
message("Climate outputs: ", raw_weather_dir)
message("Soil outputs: ", raw_soil_dir)

##### END OF FILE #####
# 1_get_wth_soil_data.R - COMPLETE VERSION WITH DOWNLOADS
# Myanmar-only DNDC data preparation: climate + soil acquisition and harmonization
# Uses NASA POWER for climate (no account needed), SoilGrids for soil

library(terra)
library(sf)
library(tidyverse)
library(nasapower)  # For climate data
library(httr)       # For SoilGrids downloads
library(lubridate)

# ============================================================================
# 0. CONFIGURATION & SETUP
# ============================================================================
base_dir <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"

# Create directory structure
dirs <- c(
  file.path(base_dir, "data/raw/climate"),
  file.path(base_dir, "data/raw/soil"),
  file.path(base_dir, "data/raw/mask"),
  file.path(base_dir, "data/processed/climate"),
  file.path(base_dir, "data/processed/soil"),
  file.path(base_dir, "data/processed/tables"),
  file.path(base_dir, "dndc_inputs/climate_files")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Simulation period
sim_years <- 2005:2024

# Myanmar boundary
mmr_boundary <- st_read(file.path(base_dir, "boundaries_shapefiles/mmr_polbnda2_adm1_250k_mimu_1.shp"))
mmr_bbox <- st_bbox(mmr_boundary)

cat("Myanmar Bounding Box:\n")
cat(sprintf("  Lat: %.2f to %.2f\n", mmr_bbox["ymin"], mmr_bbox["ymax"]))
cat(sprintf("  Lon: %.2f to %.2f\n\n", mmr_bbox["xmin"], mmr_bbox["xmax"]))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Function to download NASA POWER data for Myanmar
download_nasa_power <- function(parameter, year, bbox, output_file) {
  if (file.exists(output_file)) {
    cat(" [EXISTS]")
    return(TRUE)
  }
  
  tryCatch({
    # NASA POWER API call
    dates <- c(paste0(year, "-01-01"), paste0(year, "-12-31"))
    
    # Download point data for Myanmar grid
    lon_seq <- seq(bbox["xmin"], bbox["xmax"], by = 0.5)
    lat_seq <- seq(bbox["ymin"], bbox["ymax"], by = 0.5)
    
    # Create grid of points
    grid_points <- expand.grid(lon = lon_seq, lat = lat_seq)
    
    cat(sprintf(" [Downloading %d points]", nrow(grid_points)))
    
    # Download for each point (NASA POWER limitation)
    all_data <- list()
    for (i in 1:nrow(grid_points)) {
      pt_data <- get_power(
        community = "ag",
        lonlat = c(grid_points$lon[i], grid_points$lat[i]),
        pars = parameter,
        dates = dates,
        temporal_api = "daily"
      )
      pt_data$lon <- grid_points$lon[i]
      pt_data$lat <- grid_points$lat[i]
      all_data[[i]] <- pt_data
      
      if (i %% 10 == 0) cat(sprintf("\r    Progress: %d/%d", i, nrow(grid_points)))
    }
    
    # Combine all data
    combined_data <- bind_rows(all_data)
    
    # Convert to raster
    dates_seq <- seq(as.Date(dates[1]), as.Date(dates[2]), by = "day")
    n_days <- length(dates_seq)
    
    # Create raster stack
    r_list <- list()
    for (d in 1:n_days) {
      day_data <- combined_data %>%
        filter(YYYYMMDD == format(dates_seq[d], "%Y%m%d")) %>%
        select(lon, lat, value = all_of(parameter))
      
      # Convert to raster
      r <- rast(day_data, type = "xyz", crs = "EPSG:4326")
      r_list[[d]] <- r
    }
    
    # Stack all days
    r_stack <- rast(r_list)
    names(r_stack) <- paste0("day_", 1:n_days)
    
    # Save as NetCDF
    writeCDF(r_stack, output_file, overwrite = TRUE)
    cat(" ???\n")
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf(" [ERROR: %s]\n", e$message))
    return(FALSE)
  })
}

# Function to download SoilGrids data
download_soilgrids <- function(variable, depth = "0-5cm", bbox, output_file) {
  if (file.exists(output_file)) {
    cat(" [EXISTS]")
    return(TRUE)
  }
  
  tryCatch({
    # SoilGrids variable mapping
    sg_vars <- list(
      bd = "bdod",         # Bulk density
      clay = "clay",       # Clay content
      soc = "soc",         # Soil organic carbon
      ph = "phh2o",        # pH in H2O
      sand = "sand",       # Sand content
      silt = "silt"        # Silt content
    )
    
    sg_var <- sg_vars[[variable]]
    if (is.null(sg_var)) {
      cat(" [UNKNOWN VARIABLE]")
      return(FALSE)
    }
    
    # SoilGrids WCS endpoint
    base_url <- "https://maps.isric.org/mapserv"
    
    # Build WCS request
    wcs_request <- list(
      service = "WCS",
      version = "2.0.1",
      request = "GetCoverage",
      coverageid = paste0(sg_var, "_", depth, "_mean"),
      subset = paste0("Lat(", bbox["ymin"], ",", bbox["ymax"], ")"),
      subset = paste0("Long(", bbox["xmin"], ",", bbox["xmax"], ")"),
      format = "image/tiff"
    )
    
    cat(" [Downloading from SoilGrids]")
    
    # Download using httr
    response <- GET(base_url, query = wcs_request, timeout(300))
    
    if (status_code(response) == 200) {
      # Save to file
      writeBin(content(response, "raw"), output_file)
      cat(" ???\n")
      return(TRUE)
    } else {
      cat(sprintf(" [HTTP ERROR: %d]\n", status_code(response)))
      return(FALSE)
    }
    
  }, error = function(e) {
    cat(sprintf(" [ERROR: %s]\n", e$message))
    return(FALSE)
  })
}

# ============================================================================
# 1. DOWNLOAD AND PROCESS CLIMATE DATA (ESTABLISHES GRID)
# ============================================================================

cat("\n=== DOWNLOADING CLIMATE DATA ===\n")
cat("Using NASA POWER API (0.5?? resolution)\n\n")

# NASA POWER parameter mapping
nasa_params <- list(
  tmax = "T2M_MAX",
  tmin = "T2M_MIN",
  tavg = "T2M",
  prec = "PRECTOTCORR",
  wind = "WS2M",
  srad = "ALLSKY_SFC_SW_DWN",
  rhum = "RH2M"
)

# Conversion factors for DNDC units
conv_factors <- list(
  tmax = 1,      # ??C to ??C
  tmin = 1,      # ??C to ??C
  tavg = 1,      # ??C to ??C
  prec = 0.1,    # mm to cm
  wind = 1,      # m/s to m/s
  srad = 3.6,    # kWh/m??/day to MJ/m??/day
  rhum = 1       # % to %
)

for (year in sim_years) {
  cat(sprintf("\nYear %d:\n", year))
  
  for (vname in names(nasa_params)) {
    param <- nasa_params[[vname]]
    conv <- conv_factors[[vname]]
    
    cat(sprintf("  %s (%s)...", vname, param))
    
    # File paths
    raw_file <- file.path(base_dir, "data/raw/climate", 
                          sprintf("%s_%d_native.nc", vname, year))
    proc_file <- file.path(base_dir, "data/processed/climate",
                           sprintf("%s_%d_dndc.nc", vname, year))
    
    # Download if needed
    success <- download_nasa_power(param, year, mmr_bbox, raw_file)
    
    # Convert to DNDC units
    if (success && file.exists(raw_file)) {
      if (!file.exists(proc_file)) {
        r_native <- rast(raw_file)
        r_dndc <- r_native * conv
        writeCDF(r_dndc, proc_file, overwrite = TRUE)
        cat("  ??? Converted to DNDC units ???")
      } else {
        cat("  ??? Already processed ???")
      }
    }
    cat("\n")
  }
}

# ============================================================================
# 2. DOWNLOAD AND PROCESS SOIL DATA (HARMONIZE TO CLIMATE GRID)
# ============================================================================

cat("\n=== DOWNLOADING SOIL DATA ===\n")
cat("Using SoilGrids 250m\n\n")

# Get reference climate grid
ref_climate_file <- file.path(base_dir, "data/processed/climate", 
                              sprintf("tmax_%d_dndc.nc", sim_years[1]))

if (!file.exists(ref_climate_file)) {
  stop("ERROR: Reference climate file not found. Climate data must be downloaded first.")
}

ref_grid <- rast(ref_climate_file)[[1]]
cat(sprintf("Reference grid: %d ?? %d cells, %.4f?? resolution\n\n", 
            ncol(ref_grid), nrow(ref_grid), res(ref_grid)[1]))

# Soil variables
soil_downloads <- list(
  bd = list(sg_var = "bd", depth = "0-5cm", conv = 100),      # cg/cm?? to g/cm??
  clay = list(sg_var = "clay", depth = "0-5cm", conv = 0.001), # g/kg to fraction  
  soc = list(sg_var = "soc", depth = "0-5cm", conv = 0.001),   # dg/kg to fraction
  ph = list(sg_var = "ph", depth = "0-5cm", conv = 0.1),       # pH*10 to pH
  sand = list(sg_var = "sand", depth = "0-5cm", conv = 0.001),
  silt = list(sg_var = "silt", depth = "0-5cm", conv = 0.001)
)

for (sname in names(soil_downloads)) {
  sinfo <- soil_downloads[[sname]]
  cat(sprintf("%s...", sname))
  
  # Download raw
  raw_file <- file.path(base_dir, "data/raw/soil",
                        sprintf("%s_mmr_native.tif", sname))
  proc_file <- file.path(base_dir, "data/processed/soil",
                         sprintf("%s_dndc_on_climategrid.tif", sname))
  
  # Download
  download_soilgrids(sinfo$sg_var, sinfo$depth, mmr_bbox, raw_file)
  
  # Process if downloaded
  if (file.exists(raw_file) && !file.exists(proc_file)) {
    cat("  ??? Processing...")
    
    s_raw <- rast(raw_file)
    s_crop <- crop(s_raw, mmr_boundary)
    s_mask <- mask(s_crop, vect(mmr_boundary))
    
    # Apply conversion
    s_converted <- s_mask * sinfo$conv
    
    # Resample to climate grid
    s_final <- resample(s_converted, ref_grid, method = "bilinear")
    
    writeRaster(s_final, proc_file, overwrite = TRUE)
    cat(" ???\n")
  } else if (file.exists(proc_file)) {
    cat(" ??? Already processed ???\n")
  } else {
    cat("\n")
  }
}

# Calculate porosity from bulk density if not available
poros_file <- file.path(base_dir, "data/processed/soil", "poros_dndc_on_climategrid.tif")
bd_file <- file.path(base_dir, "data/processed/soil", "bd_dndc_on_climategrid.tif")

if (!file.exists(poros_file) && file.exists(bd_file)) {
  cat("Calculating porosity from bulk density...\n")
  bd <- rast(bd_file)
  # Porosity = 1 - (BD / 2.65), where 2.65 is particle density
  poros <- 1 - (bd / 2.65)
  writeRaster(poros, poros_file, overwrite = TRUE)
  cat("  Porosity calculated ???\n")
}

# Calculate soil texture ID from sand/silt/clay
texture_file <- file.path(base_dir, "data/processed/soil", "texture_dndc_on_climategrid.tif")
clay_file <- file.path(base_dir, "data/processed/soil", "clay_dndc_on_climategrid.tif")
sand_file <- file.path(base_dir, "data/processed/soil", "sand_dndc_on_climategrid.tif")
silt_file <- file.path(base_dir, "data/processed/soil", "silt_dndc_on_climategrid.tif")

if (!file.exists(texture_file) && all(file.exists(c(clay_file, sand_file, silt_file)))) {
  cat("Calculating soil texture class...\n")
  
  clay <- rast(clay_file) * 100  # Convert back to %
  sand <- rast(sand_file) * 100
  silt <- rast(silt_file) * 100
  
  # USDA texture classification (simplified)
  # 1=sand, 2=loamy sand, 3=sandy loam, 4=loam, 5=silt loam, 
  # 6=silt, 7=sandy clay loam, 8=clay loam, 9=silty clay loam, 10=sandy clay, 11=silty clay, 12=clay
  
  texture <- clay  # Start with clay raster structure
  texture[] <- 4   # Default to loam
  
  # Clay (> 40% clay)
  texture[clay > 40 & sand > 45] <- 10  # Sandy clay
  texture[clay > 40 & silt > 40] <- 11  # Silty clay
  texture[clay > 40 & sand <= 45 & silt <= 40] <- 12  # Clay
  
  # Clay loam (27-40% clay)
  texture[clay >= 27 & clay <= 40 & sand > 20 & sand < 45] <- 7  # Sandy clay loam
  texture[clay >= 27 & clay <= 40 & silt > 40] <- 9  # Silty clay loam
  texture[clay >= 27 & clay <= 40 & sand >= 20 & sand <= 45 & silt <= 40] <- 8  # Clay loam
  
  # Loam (7-27% clay)
  texture[clay < 27 & clay >= 7 & sand < 50 & silt >= 28] <- 5  # Silt loam
  texture[clay < 27 & clay >= 7 & silt >= 80] <- 6  # Silt
  texture[clay < 27 & clay >= 7 & sand >= 23 & sand < 52 & silt < 50] <- 4  # Loam
  texture[clay < 27 & clay >= 7 & sand >= 50] <- 3  # Sandy loam
  
  # Sand (< 7% clay)
  texture[clay < 7 & sand >= 85] <- 1  # Sand
  texture[clay < 7 & sand >= 70 & sand < 85] <- 2  # Loamy sand
  
  writeRaster(texture, texture_file, overwrite = TRUE, datatype = "INT1U")
  cat("  Soil texture class calculated ???\n")
}

# ============================================================================
# 3. PREPARE RICE MASK ON CLIMATE GRID
# ============================================================================

cat("\n=== PREPARING RICE MASK ===\n")

rice_raw_file <- file.path(base_dir, "rice mask/rice_mask_mm/spam2020_ha_rice_mmr.tif")
rice_mask_file <- file.path(base_dir, "data/raw/mask/spam2020_rice_mmr_on_climategrid.tif")

if (file.exists(rice_raw_file)) {
  cat("Processing SPAM 2020 rice mask...\n")
  
  rice <- rast(rice_raw_file)
  rice_crop <- crop(rice, mmr_boundary)
  rice_resampled <- resample(rice_crop, ref_grid, method = "bilinear")
  rice_mask <- rice_resampled > 0
  
  writeRaster(rice_mask, rice_mask_file, overwrite = TRUE)
  cat(sprintf("  Rice mask created: %d rice cells ???\n", global(rice_mask, "sum", na.rm = TRUE)$sum))
} else {
  cat("ERROR: Rice mask file not found!\n")
  cat(sprintf("Expected: %s\n", rice_raw_file))
}

# ============================================================================
# 4. CREATE CELLS TABLE
# ============================================================================

cat("\n=== CREATING CELLS TABLE ===\n")

if (file.exists(rice_mask_file)) {
  rice_mask <- rast(rice_mask_file)
  rice_area_file <- file.path(base_dir, "rice mask/rice_mask_mm/spam2020_ha_rice_mmr.tif")
  rice_area <- rast(rice_area_file) %>% 
    crop(mmr_boundary) %>% 
    resample(ref_grid, method = "bilinear")
  
  # Extract rice cells
  cells_df <- as.data.frame(rice_mask, xy = TRUE, cells = TRUE) %>%
    filter(!is.na(lyr.1) & lyr.1 > 0) %>%
    rename(cell_id = cell, lon = x, lat = y) %>%
    select(cell_id, lon, lat)
  
  # Add rice area
  rice_values <- terra::extract(rice_area, cells_df[, c("lon", "lat")])
  cells_df$rice_area <- rice_values[[2]]
  
  # Add ADM1 information
  cells_sf <- st_as_sf(cells_df, coords = c("lon", "lat"), crs = st_crs(mmr_boundary))
  cells_with_adm <- st_join(cells_sf, mmr_boundary["ADM1_EN"])
  cells_df$adm1 <- cells_with_adm$ADM1_EN
  
  # Save cells table
  cells_file <- file.path(base_dir, "data/processed/tables/cells_rice_0p5deg_mmr.rds")
  saveRDS(cells_df, cells_file)
  
  cat(sprintf("  Cells table created: %d rice cells\n", nrow(cells_df)))
  cat(sprintf("  Saved to: %s\n", cells_file))
  
  # Print summary by region
  cat("\n  Rice cells by region:\n")
  print(table(cells_df$adm1))
  
} else {
  cat("ERROR: Cannot create cells table - rice mask not found\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n=== DATA PREPARATION COMPLETE ===\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

cat("Downloaded and processed:\n")
cat(sprintf("  ??? Climate data: %d years ?? 7 variables\n", length(sim_years)))
cat("  ??? Soil data: 6+ variables\n")
cat("  ??? Rice mask: Myanmar rice cultivation areas\n")
cat(sprintf("  ??? Cells table: %d grid cells\n", nrow(cells_df)))
cat("\n")
cat("Next step: Run 2_make_climate_file.R to generate DNDC climate files\n")
cat("\n")