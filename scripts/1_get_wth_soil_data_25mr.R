##### 1_get_wth_soil_data.R - DNDC VERSION #####

suppressPackageStartupMessages({
  library(nasapower); library(terra); library(dplyr); library(geodata); library(httr)
})

# ---- Paths / config ----
path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
dir.create(path, FALSE, TRUE); setwd(path)
raw_dir <- file.path("data", "raw"); power_dir <- file.path(raw_dir, "weather", "power")
soil_dir <- file.path(raw_dir, "soil")  # NEW: Dedicated soil folder
dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(soil_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# FIXED: Use proper 0.5° grid boundaries aligned to standard grid
years <- 2005:2024
ext0 <- ext(92, 101, 8.5, 29)  # Aligned to 0.5° intervals (not 91.5, 101.5, 8, 29)

# ---- Helpers ----
.make_layer <- function(df, val){
  r <- rast(xmin = xmin(ext0), xmax = xmax(ext0),
            ymin = ymin(ext0), ymax = ymax(ext0), 
            resolution = 0.5, crs = "EPSG:4326")
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
  fn <- file.path(power_dir, sprintf("%s-2005_2024-92x101x8.5x29.nc", parm))
  if (file.exists(fn)) {
    message(sprintf("Loading cached %s", basename(fn)))
    return(rast(fn))
  }
  R <- rast(lapply(years, \(yy) .fetch_year(parm, yy))); writeCDF(R, fn, varname = parm, longname = parm, overwrite = TRUE); R
}

write_if_missing <- function(x, fn, vn, ln, un){ if (!file.exists(fn)) writeCDF(x, fn, varname = vn, longname = ln, unit = un, overwrite = TRUE) }
getL <- function(r, nm){ n <- names(r); if (nm %in% n) r[[nm]] else if (paste0(nm, "_mean") %in% n) r[[paste0(nm, "_mean")]] else {m <- grep(paste0("^", nm), n, TRUE); if (length(m)) r[[m[1]]] else stop("Missing soil layer: ", nm)} }

.check_xy <- function(xy){
  if (min(xy[["pH_0-5cm"]], na.rm = TRUE) < 2 || max(xy[["pH_0-5cm"]], na.rm = TRUE) > 10) stop("Bad pH values in cells.rds")
  if (min(xy[["clay_0-5cm"]], na.rm = TRUE) <= 0 || max(xy[["clay_0-5cm"]], na.rm = TRUE) >= 1) stop("Bad clay fraction values in cells.rds")
  if (min(xy[["bdod_5-15cm"]], na.rm = TRUE) < 0.5 || max(xy[["bdod_5-15cm"]], na.rm = TRUE) > 2.0) stop("Bad bulk density values in cells.rds")
  if (min(xy[["poros_0-5cm"]], na.rm = TRUE) < 0.2 || max(xy[["poros_0-5cm"]], na.rm = TRUE) > 0.8) stop("Bad porosity values in cells.rds")
  invisible(xy)
}

# 1) POWER climate (download/load)
message("=== POWER climate (load/download) ===")
R_TMAX <- .get_or_load_var("T2M_MAX"); R_TMIN <- .get_or_load_var("T2M_MIN"); R_TAVG <- .get_or_load_var("T2M")
R_PREC <- .get_or_load_var("PRECTOTCORR"); R_WIND <- .get_or_load_var("WS2M"); R_ALLSKY <- .get_or_load_var("ALLSKY_SFC_SW_DWN"); R_RH <- .get_or_load_var("RH2M")

# 2) DNDC units (write-once)
message("=== Converting to DNDC Units ===")
write_if_missing(R_TMAX, file.path(power_dir, "tmax-2005_2024-92x101x8.5x29.nc"), "tmax", "Daily Max Temperature", "C")
write_if_missing(R_TMIN, file.path(power_dir, "tmin-2005_2024-92x101x8.5x29.nc"), "tmin", "Daily Min Temperature", "C")
write_if_missing(R_TAVG, file.path(power_dir, "tavg-2005_2024-92x101x8.5x29.nc"), "tavg", "Daily Avg Temperature", "C")
write_if_missing(R_PREC/10, file.path(power_dir, "prec-2005_2024-92x101x8.5x29.nc"), "prec", "Daily Precipitation", "cm")
write_if_missing(R_WIND, file.path(power_dir, "wind-2005_2024-92x101x8.5x29.nc"), "wind", "Wind Speed", "m/s")
write_if_missing(R_ALLSKY*3.6, file.path(power_dir, "srad-2005_2024-92x101x8.5x29.nc"), "srad", "Solar Radiation", "MJ/m2/day")
write_if_missing(R_RH, file.path(power_dir, "rhum-2005_2024-92x101x8.5x29.nc"), "rhum", "Relative Humidity", "%")

# 3) Elevation
message("=== Processing Elevation ===")
elv_tif <- file.path(raw_dir, "elevation.tif")
if (file.exists(elv_tif)) {
  message("Removing old elevation.tif to regenerate with aligned grid")
  file.remove(elv_tif)
}
ref <- rast(file.path(power_dir, "T2M_MAX-2005_2024-92x101x8.5x29.nc"))[[1]]
writeRaster(resample(geodata::elevation_30s("Myanmar", path = raw_dir), ref, "average"), elv_tif, overwrite = TRUE)

# 4) Soil
message("=== Processing Soil Data ===")
soil_tif <- file.path(soil_dir, "soil.tif")
soil_agg_tif <- file.path(soil_dir, "soil_agg.tif")

if (file.exists(soil_agg_tif)) {
  message("Removing old soil_agg.tif to regenerate with aligned grid")
  file.remove(soil_agg_tif)
}

soil <- crop(geodata::soil_world(c("bdod","clay","sand","silt","soc","phh2o"), c(5,15,30), stat = "mean", vsi = TRUE), ext0)
writeRaster(soil, soil_tif, overwrite = TRUE)
# Convert to DNDC-ready units for script 3
soil[["bdod_0-5cm"]]  <- getL(soil, "bdod_0-5cm")          # already g/cm3
soil[["bdod_5-15cm"]] <- getL(soil, "bdod_5-15cm")
soil[["bdod_15-30cm"]] <- getL(soil, "bdod_15-30cm")

# Keep raw texture components in % for classification
clay_pct_r <- getL(soil, "clay_0-5cm")                    # already %
sand_pct_r <- getL(soil, "sand_0-5cm")                    # already %
silt_pct_r <- getL(soil, "silt_0-5cm")                    # already %

# DNDC-ready continuous variables
soil[["clay_0-5cm"]] <- clay_pct_r / 100                  # % -> fraction
soil[["sand_0-5cm"]] <- sand_pct_r / 100                  # % -> fraction
soil[["silt_0-5cm"]] <- silt_pct_r / 100                  # % -> fraction
soil[["soc_0-5cm"]]  <- getL(soil, "soc_0-5cm") / 1000    # g/kg -> kg/kg
soil[["pH_0-5cm"]]   <- getL(soil, "phh2o_0-5cm")         # already pH
soil[["poros_0-5cm"]] <- 1 - (soil[["bdod_0-5cm"]] / 2.65)

# Texture classification (uses % values, not fractions)
clay_pct <- values(clay_pct_r)
sand_pct <- values(sand_pct_r)
silt_pct <- values(silt_pct_r)

tex <- rep(4L, length(clay_pct))
tex[clay_pct > 40 & sand_pct > 45] <- 10L
tex[clay_pct > 40 & silt_pct > 40] <- 11L
tex[clay_pct > 40 & sand_pct <= 45 & silt_pct <= 40] <- 12L
tex[clay_pct >= 27 & clay_pct <= 40 & sand_pct > 20 & sand_pct < 45] <- 7L
tex[clay_pct >= 27 & clay_pct <= 40 & silt_pct > 40] <- 9L
tex[clay_pct >= 27 & clay_pct <= 40 & sand_pct >= 20 & sand_pct <= 45 & silt_pct <= 40] <- 8L
tex[clay_pct < 27 & clay_pct >= 7 & sand_pct < 50 & silt_pct >= 28] <- 5L
tex[clay_pct < 27 & clay_pct >= 7 & silt_pct >= 80] <- 6L
tex[clay_pct < 27 & clay_pct >= 7 & sand_pct >= 23 & sand_pct < 52 & silt_pct < 50] <- 4L
tex[clay_pct < 27 & clay_pct >= 7 & sand_pct >= 50] <- 3L
tex[clay_pct < 7 & sand_pct >= 85] <- 1L
tex[clay_pct < 7 & sand_pct >= 70 & sand_pct < 85] <- 2L

rtx <- soil[["bdod_0-5cm"]]
values(rtx) <- tex
names(rtx) <- "texture_0-5cm"
soil <- c(soil, rtx)

cont_names <- c("bdod_0-5cm", "bdod_5-15cm", "bdod_15-30cm",
                "clay_0-5cm", "sand_0-5cm", "silt_0-5cm",
                "soc_0-5cm", "pH_0-5cm", "poros_0-5cm")
soil_cont <- soil[[cont_names]]
soil_tex  <- soil[["texture_0-5cm"]]

soil_resampled <- c(
  resample(soil_cont, rast(elv_tif), "average"),
  resample(soil_tex,  rast(elv_tif), "near")
)

writeRaster(soil_resampled, soil_agg_tif, overwrite = TRUE)
soil2 <- rast(soil_agg_tif)
# 5) Rice mask + cells table
message("=== Creating Rice Cells Table ===")
cells_rds <- file.path("data", "cells.rds")
aoi <- geodata::gadm("Myanmar", level = 1, path = raw_dir)
elv <- rast(elv_tif)

r <- mask(init(elv, "cell"), aoi, touches = TRUE)
rice <- geodata::crop_spam(crop = "rice", var = "area", path = raw_dir) |> crop(aoi, mask = TRUE)
r <- mask(r, resample(rice[[1]], r, "sum") > 500, maskvalue = FALSE)

cells <- data.frame(r, na.rm = TRUE)[, 1]
xy <- as.data.frame(xyFromCell(r, cells))
xy$elevation <- round(elv[cells])
xy$cell <- cells

# extract soil values; do NOT drop the first column because extract() is already
# returning only raster layers for numeric cell indices here
xy <- cbind(xy, extract(soil2, cells))

xy <- xy[complete.cases(xy[c("pH_0-5cm", "clay_0-5cm", "bdod_0-5cm")]), ]

# sanity check that script 3 will find the expected DNDC-ready columns
stopifnot(all(c("bdod_0-5cm", "clay_0-5cm", "sand_0-5cm", "silt_0-5cm",
                "soc_0-5cm", "pH_0-5cm", "poros_0-5cm", "texture_0-5cm") %in% names(xy)))

.check_xy(xy)
saveRDS(xy, cells_rds)

message("\n=== SUMMARY ===")
message(sprintf("Climate data: 2005-2024 (%d years)", length(years)))
message("Variables: tmax, tmin, tavg, prec, wind, srad, rhum")
message("Grid resolution: 0.5° x 0.5°")
message(sprintf("Extent: X=[%.1f, %.1f], Y=[%.1f, %.1f]", xmin(ext0), xmax(ext0), ymin(ext0), ymax(ext0)))
message(sprintf("Rice cells: %d", nrow(xy)))
message(sprintf("Cells table saved: %s", cells_rds))
message(sprintf("pH range: %.2f - %.2f", min(xy[["pH_0-5cm"]], na.rm = TRUE), max(xy[["pH_0-5cm"]], na.rm = TRUE)))
message(sprintf("Clay range: %.3f - %.3f", min(xy[["clay_0-5cm"]], na.rm = TRUE), max(xy[["clay_0-5cm"]], na.rm = TRUE)))
message(sprintf("Bulk density range: %.3f - %.3f", min(xy[["bdod_0-5cm"]], na.rm = TRUE), max(xy[["bdod_0-5cm"]], na.rm = TRUE)))
message(sprintf("Porosity range: %.3f - %.3f", min(xy[["poros_0-5cm"]], na.rm = TRUE), max(xy[["poros_0-5cm"]], na.rm = TRUE)))

# Verify grid alignment
x_unique <- sort(unique(xy$x))
y_unique <- sort(unique(xy$y))
x_spacing <- unique(round(diff(x_unique), 6))
y_spacing <- unique(round(diff(y_unique), 6))
message(sprintf("\nGrid verification:"))
message(sprintf("  X spacing: %.6f° (expected: 0.5)", x_spacing[1]))
message(sprintf("  Y spacing: %.6f° (expected: 0.5)", y_spacing[1]))
if (all(abs(x_spacing - 0.5) < 0.001) && all(abs(y_spacing - 0.5) < 0.001)) {
  message("  Status: ✅ Grid aligned to 0.5°")
} else {
  message("  Status: ⚠️ Mixed spacing detected")
}

message("\n✓ Step 1 complete. Ready for Step 2 (make_climate_files.R)")