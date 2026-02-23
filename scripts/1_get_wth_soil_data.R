##### 1_get_wth_soil_data.R - DNDC VERSION #####
# Succinct, reproducible workflow (same steps, outputs, and file names).

suppressPackageStartupMessages({
  library(nasapower); library(terra); library(dplyr); library(geodata); library(httr)
})

# ---- Paths / config ----
path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
dir.create(path, FALSE, TRUE); setwd(path)
raw_dir <- file.path("data", "raw"); power_dir <- file.path(raw_dir, "weather", "power")
dir.create(power_dir, recursive = TRUE, showWarnings = FALSE); dir.create("data", showWarnings = FALSE)
years <- 2005:2024; ext0 <- ext(91.5, 101.5, 8, 29)
lon_seq <- seq(xmin(ext0), xmax(ext0), by = 0.5); lat_seq <- sort(seq(ymin(ext0), ymax(ext0), by = 0.5), TRUE)

# ---- Helpers ----
.make_layer <- function(df, val){
  r <- rast(ncols = length(lon_seq), nrows = length(lat_seq), xmin = xmin(ext0), xmax = xmax(ext0),
            ymin = ymin(ext0), ymax = ymax(ext0), crs = "EPSG:4326")
  values(r) <- (df |> arrange(desc(LAT), LON))[[val]]; r
}

.fetch_year <- function(parm, yr){
  httr::set_config(httr::timeout(300)); on.exit(httr::reset_config(), add = TRUE)
  d1 <- paste0(yr, "-01-01"); d2 <- paste0(yr, "-12-31"); tile <- 4; min_tile <- 2
  lon0s <- seq(xmin(ext0), xmax(ext0), by = tile); lat0s <- seq(ymin(ext0), ymax(ext0), by = tile)
  if (length(lon0s) > 1 && xmax(ext0) - tail(lon0s, 1) < min_tile) lon0s <- head(lon0s, -1)
  if (length(lat0s) > 1 && ymax(ext0) - tail(lat0s, 1) < min_tile) lat0s <- head(lat0s, -1)
  out <- list()
  for (x0 in lon0s) for (y0 in lat0s){
    x1 <- min(x0 + tile, xmax(ext0)); y1 <- min(y0 + tile, ymax(ext0))
    if ((x1 - x0) < min_tile || (y1 - y0) < min_tile) next
    bb <- c(x0, y0, x1, y1)
    for (k in 1:5){
      ok <- tryCatch({
        out[[length(out) + 1]] <- get_power("ag", lonlat = bb, pars = parm, dates = c(d1, d2), temporal_api = "daily")
        Sys.sleep(3); TRUE
      }, error = function(e){
        m <- as.character(e)
        if (grepl("Timeout|timed out", m, TRUE)) Sys.sleep(30 * k) else if (grepl("429|rate limit", m, TRUE)) Sys.sleep(120 * k) else Sys.sleep(20)
        FALSE
      })
      if (ok) break
    }
  }
  if (!length(out)) stop("No tiles downloaded for ", parm, " ", yr)
  dat <- bind_rows(out) |> rename(DATE = YYYYMMDD); dat$DATE <- as.Date(dat$DATE)
  dat <- dat[!duplicated(dat[, c("LAT", "LON", "DATE")]), ]; days <- sort(unique(dat$DATE))
  rs <- lapply(days, \(d) .make_layer(dat[dat$DATE == d, c("LAT", "LON", parm)], parm))
  r <- rast(rs); time(r) <- days; names(r) <- paste0(parm, "_", format(days, "%Y%j")); r
}

.get_or_load_var <- function(parm){
  fn <- file.path(power_dir, sprintf("%s-2005_2024-91.5x101.5x8x29.nc", parm))
  if (file.exists(fn)) return(rast(fn))
  R <- rast(lapply(years, \(yy) .fetch_year(parm, yy))); writeCDF(R, fn, varname = parm, longname = parm, overwrite = TRUE); R
}

write_if_missing <- function(x, fn, vn, ln, un){ if (!file.exists(fn)) writeCDF(x, fn, varname = vn, longname = ln, unit = un, overwrite = TRUE) }
getL <- function(r, nm){ n <- names(r); if (nm %in% n) r[[nm]] else if (paste0(nm, "_mean") %in% n) r[[paste0(nm, "_mean")]] else {m <- grep(paste0("^", nm), n, TRUE); if (length(m)) r[[m[1]]] else stop("Missing soil layer: ", nm)} }

# 1) POWER climate (download/load)
message("=== POWER climate (load/download) ===")
R_TMAX <- .get_or_load_var("T2M_MAX"); R_TMIN <- .get_or_load_var("T2M_MIN"); R_TAVG <- .get_or_load_var("T2M")
R_PREC <- .get_or_load_var("PRECTOTCORR"); R_WIND <- .get_or_load_var("WS2M"); R_ALLSKY <- .get_or_load_var("ALLSKY_SFC_SW_DWN"); R_RH <- .get_or_load_var("RH2M")

# 2) DNDC units (write-once)
message("=== Converting to DNDC Units ===")
write_if_missing(R_TMAX, file.path(power_dir, "tmax-2005_2024-91.5x101.5x8x29.nc"), "tmax", "Daily Max Temperature", "C")
write_if_missing(R_TMIN, file.path(power_dir, "tmin-2005_2024-91.5x101.5x8x29.nc"), "tmin", "Daily Min Temperature", "C")
write_if_missing(R_TAVG, file.path(power_dir, "tavg-2005_2024-91.5x101.5x8x29.nc"), "tavg", "Daily Avg Temperature", "C")
write_if_missing(R_PREC/10, file.path(power_dir, "prec-2005_2024-91.5x101.5x8x29.nc"), "prec", "Daily Precipitation", "cm")
write_if_missing(R_WIND, file.path(power_dir, "wind-2005_2024-91.5x101.5x8x29.nc"), "wind", "Wind Speed", "m/s")
write_if_missing(R_ALLSKY*3.6, file.path(power_dir, "srad-2005_2024-91.5x101.5x8x29.nc"), "srad", "Solar Radiation", "MJ/m2/day")
write_if_missing(R_RH, file.path(power_dir, "rhum-2005_2024-91.5x101.5x8x29.nc"), "rhum", "Relative Humidity", "%")

# 3) Elevation
message("=== Processing Elevation ===")
elv_tif <- file.path(raw_dir, "elevation.tif")
if (!file.exists(elv_tif)){
  ref <- rast(file.path(power_dir, "T2M_MAX-2005_2024-91.5x101.5x8x29.nc"))[[1]]
  writeRaster(resample(geodata::elevation_30s("Myanmar", path = raw_dir), ref, "average"), elv_tif, overwrite = TRUE)
}

# 4) Soil
message("=== Processing Soil Data ===")
soil_tif <- file.path(raw_dir, "soil.tif"); soil_agg_tif <- file.path(raw_dir, "soil_agg.tif")
if (!file.exists(soil_agg_tif)){
  soil <- if (!file.exists(soil_tif)) crop(geodata::soil_world(c("bdod","clay","sand","silt","soc","phh2o"), c(5,15,30), stat = "mean", vsi = TRUE), ext0, filename = soil_tif) else rast(soil_tif)
  soil[["bdod_0-5cm"]] <- getL(soil, "bdod_0-5cm")/100; soil[["bdod_5-15cm"]] <- getL(soil, "bdod_5-15cm")/100; soil[["bdod_15-30cm"]] <- getL(soil, "bdod_15-30cm")/100
  soil[["clay_0-5cm"]] <- getL(soil, "clay_0-5cm")/1000; soil[["sand_0-5cm"]] <- getL(soil, "sand_0-5cm")/1000; soil[["silt_0-5cm"]] <- getL(soil, "silt_0-5cm")/1000
  soil[["soc_0-5cm"]] <- getL(soil, "soc_0-5cm")/1000; soil[["pH_0-5cm"]] <- getL(soil, "phh2o_0-5cm")/10; soil[["poros_0-5cm"]] <- 1 - (soil[["bdod_0-5cm"]]/2.65)
  clay <- values(soil[["clay_0-5cm"]]) * 100; sand <- values(soil[["sand_0-5cm"]]) * 100; silt <- values(soil[["silt_0-5cm"]]) * 100; tex <- rep(4L, length(clay))
  tex[clay > 40 & sand > 45] <- 10L; tex[clay > 40 & silt > 40] <- 11L; tex[clay > 40 & sand <= 45 & silt <= 40] <- 12L; tex[clay >= 27 & clay <= 40 & sand > 20 & sand < 45] <- 7L
  tex[clay >= 27 & clay <= 40 & silt > 40] <- 9L; tex[clay >= 27 & clay <= 40 & sand >= 20 & sand <= 45 & silt <= 40] <- 8L; tex[clay < 27 & clay >= 7 & sand < 50 & silt >= 28] <- 5L
  tex[clay < 27 & clay >= 7 & silt >= 80] <- 6L; tex[clay < 27 & clay >= 7 & sand >= 23 & sand < 52 & silt < 50] <- 4L; tex[clay < 27 & clay >= 7 & sand >= 50] <- 3L
  tex[clay < 7 & sand >= 85] <- 1L; tex[clay < 7 & sand >= 70 & sand < 85] <- 2L
  tx <- soil[["clay_0-5cm"]]; values(tx) <- tex; names(tx) <- "texture_0-5cm"; soil <- c(soil, tx)
  writeRaster(resample(soil, rast(elv_tif), "average"), soil_agg_tif, overwrite = TRUE)
}
soil2 <- rast(soil_agg_tif)

# 5) Rice mask + cells table (ORYZA-consistent logic)
message("=== Creating Rice Cells Table ===")
cells_rds <- file.path("data", "cells.rds")
if (!file.exists(cells_rds)){
  aoi <- geodata::gadm("Myanmar", level = 1, path = raw_dir); elv <- rast(elv_tif)
  r <- mask(init(elv, "cell"), aoi, touches = TRUE)
  rice <- geodata::crop_spam(crop = "rice", var = "area", path = raw_dir) |> crop(aoi, mask = TRUE)
  r <- mask(r, resample(rice[[1]], r, "sum") > 0, maskvalue = FALSE)
  cells <- data.frame(r, na.rm = TRUE)[,1]; xy <- as.data.frame(xyFromCell(r, cells)); xy$elevation <- round(elv[cells]); xy$cell <- cells
  xy <- cbind(xy, extract(soil2, cells)[, -1, drop = FALSE]); saveRDS(xy, cells_rds)
} else xy <- readRDS(cells_rds)

message("\n=== SUMMARY ===")
message(sprintf("Climate data: 2005-2024 (%d years)", length(years)))
message("Variables: tmax, tmin, tavg, prec, wind, srad, rhum")
message("Grid resolution: 0.5° x 0.5°")
message(sprintf("Rice cells: %d", nrow(xy)))
message(sprintf("Cells table saved: %s", cells_rds))
message("\n✓ Step 1 complete. Ready for Step 2 (make_climate_files.R)")
##### END OF: 1_get_wth_soil_data.R #####
