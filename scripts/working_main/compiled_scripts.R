# ---- 1_get_wth_soil_data_25mr.R ----
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

# ---- 10_figures_maps_monsoon.R ----
###############################################################################
# Script 9+10 — Myanmar CH4 Figures & Spatial Maps | MONSOON ONLY
# DOY 166–309 (Jun 15 – Nov 5)  |  13 scenarios  |  263 grid cells
#
# INPUTS : enriched_daily_panel.rds
# OUTPUTS: PNG figures -> FIGURE_DIR (daily/flux) and MAP_DIR (spatial)
#
# FLUX FIGURES (F1–F7)
#   F1  All 13 scenarios — monsoon daily CH4, single panel
#   F2  Faceted by scenario category vs baseline
#   F3  CH4 components (prod / oxid / net) — Baseline, AWD, Hi-Residue
#   F4  Interannual variability 2005-2024 (3 key scenarios)
#   F5  Heatmap: scenario × DOY (7-day rolling mean)
#   F6  Flooding window overlay (mechanism view)
#   F7  Monsoon cumulative CH4 — scenarios ranked by season-end total
#
# SPATIAL MAPS (M1–M6)
#   M1  Baseline mean monsoon CH4
#   M2  All scenarios — small-multiple tile maps
#   M3  Delta vs baseline (diverging)
#   M4  AWD mitigation — absolute + % reduction
#   M5  Coefficient of variation (interannual variability)
#   M6  Optimal (lowest-emitting) scenario per cell
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
###############################################################################
# Script 9+10 — Myanmar CH4 Figures & Spatial Maps | MONSOON ONLY
# DOY 166–309 (Jun 15 – Nov 5)  |  13 scenarios  |  263 grid cells
#
# INPUTS : enriched_daily_panel.rds
# OUTPUTS: PNG figures -> FIGURE_DIR (daily/flux) and MAP_DIR (spatial)
#
# FLUX FIGURES (F1–F7)
#   F1  All 13 scenarios — monsoon daily CH4, single panel
#   F2  Faceted by scenario category vs baseline
#   F3  CH4 components (prod / oxid / net) — Baseline, AWD, Hi-Residue
#   F4  Interannual variability 2005-2024 (3 key scenarios)
#   F5  Heatmap: scenario × DOY (7-day rolling mean)
#   F6  Flooding window overlay (mechanism view)
#   F7  Monsoon cumulative CH4 — scenarios ranked by season-end total
#
# SPATIAL MAPS (M1–M6)
#   M1  Baseline mean monsoon CH4
#   M2  All scenarios — small-multiple tile maps
#   M3  Delta vs baseline (diverging)
#   M4  AWD mitigation — absolute + % reduction
#   M5  Coefficient of variation (interannual variability)
#   M6  Optimal (lowest-emitting) scenario per cell
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =============================================================================
# 0. PATHS & CONSTANTS
# =============================================================================

PANEL_RDS  <- "C:/DNDC/simulate_dndc/dndc_outputs/enriched_daily_panel.rds"
FIGURE_DIR <- "C:/DNDC/simulate_dndc/figures/monsoon"
MAP_DIR    <- "C:/DNDC/simulate_dndc/figures/monsoon/maps"
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAP_DIR,    showWarnings = FALSE, recursive = TRUE)

MONSOON    <- c(166L, 309L)   # DOY window
TILE_SIZE  <- 0.48            # slightly < 0.5° to leave gaps between tiles

stopifnot("enriched_daily_panel.rds not found" = file.exists(PANEL_RDS))

# =============================================================================
# 1. METADATA
# =============================================================================

SCENE_ORDER <- c("S000","S001","S002","S003","S004","S005","S006",
                 "S008","S009","S010","S011","S012","S013")

SCENE_SHORT <- c(
  S000="Baseline",   S001="Inj.Fert",   S002="Med-Hi Fert",
  S003="Ext.Fert",   S004="Split App",  S005="Short Flood",
  S006="AWD",        S008="No-Till",    S009="Deep Plow",
  S010="Cont.Irrig", S011="Rainfed",    S012="Low Residue",
  S013="Hi Residue"
)

SCENE_COLS <- c(
  S000="#2C3E50", S001="#1565C0", S002="#42A5F5", S003="#00BCD4",
  S004="#80DEEA", S005="#2E7D32", S006="#66BB6A", S008="#E65100",
  S009="#FF8F00", S010="#AB47BC", S011="#CE93D8", S012="#C62828",
  S013="#EF9A9A"
)

CAT_COLS <- c(
  baseline="#2C3E50", fertilizer="#1565C0", water="#2E7D32",
  tillage="#E65100",  irrigation="#AB47BC",  residue="#C62828"
)

# Month x-axis within monsoon window
MON_BREAKS <- c(166, 182, 213, 244, 274, 305, 309)
MON_LABELS <- c("Jun 15","Jul","Aug","Sep","Oct","Nov","Nov 5")

# =============================================================================
# 2. SHARED HELPERS
# =============================================================================

theme_ch4 <- function(b = 11) {
  theme_bw(base_size = b) +
    theme(plot.title = element_text(face = "bold", size = b + 2),
          plot.subtitle = element_text(size = b - 1, colour = "grey40"),
          plot.caption  = element_text(size = b - 2, colour = "grey50", hjust = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          strip.background = element_rect(fill = "grey93", colour = "grey70"),
          strip.text = element_text(face = "bold", size = b - 1))
}

theme_map <- function(base = 10) {
  b <- base
  theme_bw(base_size = b) +
    theme(plot.title = element_text(face = "bold", size = b + 2),
          plot.subtitle = element_text(size = b - 1, colour = "grey40"),
          plot.caption  = element_text(size = b - 2, colour = "grey50", hjust = 0),
          panel.grid.major = element_line(colour = "grey88", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          legend.key.width = unit(0.5, "cm"),
          axis.text  = element_text(size = b - 2),
          axis.title = element_blank(),
          strip.background = element_rect(fill = "grey90", colour = "grey70"),
          strip.text = element_text(face = "bold", size = b - 2))
}

save_png <- function(p, name, dir = FIGURE_DIR, w = 14, h = 8, dpi = 300) {
  path <- file.path(dir, name)
  ggsave(path, p, width = w, height = h, dpi = dpi, device = "png", bg = "white")
  message(sprintf("  [PNG] %s  (%d KB)", name,
                  round(file.info(path)$size / 1024)))
  invisible(path)
}

ch4_fill_scale <- function(name = "Monsoon CH4\n(kg C/ha)", ...) {
  scale_fill_gradientn(
    colours = c("#FFFDE7","#FFF176","#FFB300","#E65100","#B71C1C"),
    name = name, labels = label_number(accuracy = 1), ...)
}

mon_x <- function() {
  scale_x_continuous(breaks = MON_BREAKS, labels = MON_LABELS, expand = c(0.01, 0))
}

# =============================================================================
# 3. LOAD & AGGREGATE — MONSOON ONLY
# =============================================================================

message("Loading panel ...")
panel <- readRDS(PANEL_RDS)
setnames(panel,
         old = c("CH4-flux","CH4-prod.","CH4-oxid.","CH4-pool"),
         new = c("CH4_flux","CH4_prod","CH4_oxid","CH4_pool"),
         skip_absent = TRUE)

stopifnot("CH4_flux missing" = "CH4_flux" %in% names(panel),
          "'x' missing"      = "x"        %in% names(panel),
          "'y' missing"      = "y"        %in% names(panel))

panel <- panel[Day >= MONSOON[1] & Day <= MONSOON[2]]
message(sprintf("  Monsoon panel: %s rows", format(nrow(panel), big.mark = ",")))

# Optional CH4_prod / CH4_oxid guard
if (!"CH4_prod" %in% names(panel)) panel[, CH4_prod := NA_real_]
if (!"CH4_oxid" %in% names(panel)) panel[, CH4_oxid := NA_real_]

# ── Agg A: scene × DOY ────────────────────────────────────────────────────────
message("Aggregating: scene × Day ...")
agg <- panel[, .(
  ch4_med  = median(CH4_flux, na.rm = TRUE),
  ch4_mean = mean(  CH4_flux, na.rm = TRUE),
  ch4_q25  = quantile(CH4_flux, 0.25, na.rm = TRUE),
  ch4_q75  = quantile(CH4_flux, 0.75, na.rm = TRUE),
  ch4_q10  = quantile(CH4_flux, 0.10, na.rm = TRUE),
  ch4_q90  = quantile(CH4_flux, 0.90, na.rm = TRUE),
  prod_med  = median(CH4_prod, na.rm = TRUE),
  oxid_med  = median(CH4_oxid, na.rm = TRUE)
), by = .(scene_id, scene_category, Day)]

agg[, scene_id    := factor(scene_id, levels = SCENE_ORDER)]
agg[, scene_short := factor(SCENE_SHORT[as.character(scene_id)],
                            levels = SCENE_SHORT[SCENE_ORDER])]
agg[, scene_category := factor(scene_category,
                               levels = c("baseline","fertilizer","water",
                                          "tillage","irrigation","residue"))]

# ── Agg B: scene × year × DOY (3 focal scenarios only) ────────────────────────
message("Aggregating: scene × year × Day ...")
agg_yr <- panel[scene_id %in% c("S000","S006","S013"), .(
  ch4_med = median(CH4_flux, na.rm = TRUE),
  ch4_q25 = quantile(CH4_flux, 0.25, na.rm = TRUE),
  ch4_q75 = quantile(CH4_flux, 0.75, na.rm = TRUE)
), by = .(scene_id, cal_year, Day)]

# ── Agg C: cell-level seasonal totals (for maps) ──────────────────────────────
message("Aggregating: cell × scene × year totals ...")
cell_ann <- panel[, .(
  ann_CH4 = sum(CH4_flux, na.rm = TRUE)
), by = .(cell_id, x, y, scene_id, cal_year)]

map_data <- cell_ann[, .(
  mean_CH4 = mean(ann_CH4, na.rm = TRUE),
  sd_CH4   = sd(ann_CH4,   na.rm = TRUE),
  cv_CH4   = sd(ann_CH4, na.rm = TRUE) / mean(ann_CH4, na.rm = TRUE) * 100,
  n_years  = .N
), by = .(cell_id, x, y, scene_id)]

map_data[, scene_id := as.character(scene_id)]
map_data[, scene_short := factor(SCENE_SHORT[scene_id], levels = SCENE_SHORT[SCENE_ORDER])]

# Delta vs S000
base_map <- map_data[scene_id == "S000", .(cell_id, x, y, base_CH4 = mean_CH4)]
map_data  <- merge(map_data, base_map, by = c("cell_id","x","y"), all.x = TRUE)
map_data[, delta_CH4 := mean_CH4 - base_CH4]
map_data[, pct_change := delta_CH4 / base_CH4 * 100]

rm(panel, cell_ann); gc()
message("  Panel freed. Aggregation complete.\n")

# =============================================================================
# 4. MYANMAR BORDER
# =============================================================================

mmr_border <- tryCatch({
  if (!requireNamespace("maps", quietly = TRUE))
    install.packages("maps", repos = "https://cloud.r-project.org")
  library(maps)
  bd <- map_data("world", region = "Myanmar")
  if (!nrow(bd)) map_data("world", region = "Burma") else bd
}, error = function(e) {
  message("  WARN: Myanmar border unavailable — maps will lack outline")
  data.frame(long = numeric(0), lat = numeric(0), group = integer(0))
})

x_range <- range(map_data$x, na.rm = TRUE)
y_range <- range(map_data$y, na.rm = TRUE)
xlim    <- c(max(92.0, x_range[1] - 1.5), min(102.0, x_range[2] + 1.5))
ylim    <- c(max( 9.5, y_range[1] - 1.0), min( 29.0, y_range[2] + 1.0))

base_map_layers <- function() list(
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = "grey95", colour = "grey50", linewidth = 0.35),
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE),
  scale_x_continuous(breaks = seq(92, 102, 2),
                     labels = paste0(seq(92, 102, 2), "°E")),
  scale_y_continuous(breaks = seq(10,  30,  4),
                     labels = paste0(seq(10,  30, 4), "°N"))
)

# =============================================================================
# F1 — All scenarios, monsoon daily CH4, single panel
# =============================================================================
message("F1: All-scenario monsoon flux ...")

p1 <- ggplot(agg, aes(x = Day, colour = scene_id, fill = scene_id)) +
  geom_ribbon(aes(ymin = ch4_q25, ymax = ch4_q75), alpha = 0.10, colour = NA) +
  geom_line(aes(y = ch4_med), linewidth = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_fill_manual(  values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  mon_x() +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(title    = "Monsoon Daily CH4 Flux — All 13 Scenarios",
       subtitle = "Median across 263 cells × 20 years. Ribbons = IQR (25th–75th pct). DOY 166–309.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Source: DNDC 9.5 batch simulation | Myanmar irrigated rice") +
  theme_ch4() + theme(legend.position = "right")

save_png(p1, "F1_monsoon_CH4_all_scenarios.png")

# =============================================================================
# F2 — Faceted by category, baseline reference in each facet
# =============================================================================
message("F2: Faceted by category ...")

baseline_ref <- agg[scene_id == "S000", .(Day, ch4_baseline = ch4_med)]

p2 <- ggplot(agg, aes(x = Day, colour = scene_id, fill = scene_id)) +
  geom_line(data = baseline_ref, aes(x = Day, y = ch4_baseline),
            inherit.aes = FALSE, colour = "grey50", linewidth = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = ch4_q25, ymax = ch4_q75), alpha = 0.18, colour = NA) +
  geom_line(aes(y = ch4_med), linewidth = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_fill_manual(  values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  mon_x() +
  facet_wrap(~ scene_category, ncol = 3, scales = "free_y") +
  labs(title    = "Monsoon CH4 Flux by Management Category",
       subtitle = "Dashed grey = Baseline (S000). Ribbons = IQR. Scales free within category.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Source: DNDC 9.5 batch simulation | Myanmar irrigated rice") +
  theme_ch4() + theme(legend.position = "bottom")

save_png(p2, "F2_monsoon_CH4_by_category.png", w = 15, h = 10)

# =============================================================================
# F3 — CH4 components: production / oxidation / net (Baseline, AWD, Hi-Residue)
# =============================================================================
message("F3: CH4 components ...")

focal <- c(S000 = "Baseline (S000)", S006 = "AWD (S006)", S013 = "High Residue (S013)")
comp_long <- melt(
  agg[scene_id %in% names(focal),
      .(scene_id, Day,
        "Net flux"   = ch4_med,
        "Production" = prod_med,
        "Oxidation"  = oxid_med)],
  id.vars = c("scene_id","Day"), variable.name = "component", value.name = "ch4"
)
comp_long[, scene_label := focal[as.character(scene_id)]]
comp_long[, component   := factor(component, levels = c("Production","Oxidation","Net flux"))]

p3 <- ggplot(comp_long, aes(x = Day, y = ch4, colour = component,
                            linetype = component, linewidth = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line() +
  scale_colour_manual(  values = c(Production="#C62828", Oxidation="#1565C0", "Net flux"="#2C3E50"), name = NULL) +
  scale_linetype_manual(values = c(Production="solid",   Oxidation="dashed",  "Net flux"="solid"),  name = NULL) +
  scale_linewidth_manual(values = c(Production=0.9,      Oxidation=0.7,       "Net flux"=1.2),      name = NULL) +
  mon_x() +
  facet_wrap(~ scene_label, ncol = 3) +
  labs(title    = "CH4 Production, Oxidation, and Net Flux — Monsoon",
       subtitle = "Red solid = gross production | Blue dashed = oxidation | Dark = net flux",
       x = "Day of year", y = "CH4 (kg C/ha/day)",
       caption = "Source: DNDC 9.5 | 263 cells × 20 years median") +
  theme_ch4() + theme(legend.position = "bottom")

save_png(p3, "F3_monsoon_CH4_components.png", w = 15, h = 6)

# =============================================================================
# F4 — Interannual variability 2005-2024 (3 scenarios)
# =============================================================================
message("F4: Interannual variability ...")

YR_PAL <- colorRampPalette(c("#1565C0","#66BB6A","#C62828"))(20)
names(YR_PAL) <- as.character(2005:2024)
agg_yr[, scene_label := focal[as.character(scene_id)]]

p4 <- ggplot(agg_yr, aes(x = Day, colour = factor(cal_year))) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(aes(y = ch4_med), linewidth = 0.5, alpha = 0.85) +
  scale_colour_manual(values = YR_PAL, name = "Year") +
  mon_x() +
  facet_wrap(~ scene_label, ncol = 3, scales = "free_y") +
  labs(title    = "Interannual Variability of Monsoon CH4 (2005-2024)",
       subtitle = "Each line = one year's median across 263 cells. Blue→green→red: early→mid→late years.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Source: DNDC 9.5 | Baseline, AWD, High Residue") +
  theme_ch4() + theme(legend.key.height = unit(0.35, "cm"),
                      legend.text = element_text(size = 8))

save_png(p4, "F4_monsoon_CH4_interannual.png", w = 16, h = 6)

# =============================================================================
# F5 — Heatmap: scenario × DOY (7-day rolling mean)
# =============================================================================
message("F5: Heatmap scene × DOY ...")

agg_sm <- copy(agg)
agg_sm[, ch4_7d := frollmean(ch4_mean, 7, align = "center", na.rm = TRUE), by = scene_id]

# Order scenarios by peak monsoon emission
peak_ord <- agg[, .(peak = max(ch4_med, na.rm = TRUE)), by = scene_id][order(-peak)]$scene_id
agg_sm[, scene_ord := factor(
  scene_id,
  levels = rev(as.character(peak_ord)),
  labels = rev(paste0(as.character(peak_ord), "  ", SCENE_SHORT[as.character(peak_ord)]))
)]

p5 <- ggplot(agg_sm, aes(x = Day, y = scene_ord, fill = ch4_7d)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2980B9", mid = "#F9E79F", high = "#C0392B",
                       midpoint = 0.05, na.value = "grey85",
                       name = "Mean CH4\n(kg C/ha/d)",
                       labels = label_number(accuracy = 0.01)) +
  mon_x() +
  labs(title    = "CH4 Flux Heatmap — Scenario × Monsoon Day",
       subtitle = "7-day rolling mean | Scenarios ranked by peak emission (high→low top→bottom)",
       x = "Day of year", y = NULL,
       caption = "Source: DNDC 9.5 | 263 cells × 20 years") +
  theme_ch4() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey60", fill = NA),
        axis.text.y = element_text(size = 9))

save_png(p5, "F5_monsoon_CH4_heatmap.png", w = 13, h = 7)

# =============================================================================
# F6 — All scenarios with flooding window overlay
# =============================================================================
message("F6: CH4 flux with flooding windows ...")

flood_ann <- data.frame(
  xmin = c(185, 247), xmax = c(226, 290),
  label = c("Flood period 1", "Flood period 2")
)

p6 <- ggplot(agg, aes(x = Day, colour = scene_id, fill = scene_id)) +
  annotate("rect", xmin = flood_ann$xmin, xmax = flood_ann$xmax,
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#1565C0") +
  annotate("text", x = (flood_ann$xmin + flood_ann$xmax) / 2, y = Inf,
           label = flood_ann$label, vjust = 1.4, size = 2.8,
           colour = "#1565C0", fontface = "italic") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_ribbon(aes(ymin = ch4_q10, ymax = ch4_q90), alpha = 0.05, colour = NA) +
  geom_ribbon(aes(ymin = ch4_q25, ymax = ch4_q75), alpha = 0.12, colour = NA) +
  geom_line(aes(y = ch4_med), linewidth = 0.75) +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_fill_manual(  values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  mon_x() +
  labs(title    = "Monsoon CH4 Flux with Flooding Windows",
       subtitle = "Blue shading = monsoon flooding periods (S000 baseline). IQR + 10-90th pct ribbons.",
       x = "Day of year", y = "CH4 flux (kg C/ha/day)",
       caption = "Flooding activates methanogenesis; AWD reduces emission within same window | DNDC 9.5") +
  theme_ch4() + theme(legend.position = "right")

save_png(p6, "F6_monsoon_CH4_flooding_overlay.png", w = 16, h = 8)

# =============================================================================
# F7 — Cumulative monsoon CH4 by scenario (ranked)
# =============================================================================
message("F7: Cumulative monsoon CH4 ...")

agg_cum <- copy(agg)[order(scene_id, Day)]
agg_cum[, cum_ch4 := cumsum(ch4_med), by = scene_id]

# Final totals for ranking
finals <- agg_cum[, .(total = max(cum_ch4, na.rm = TRUE)), by = .(scene_id, scene_short)]
setorder(finals, -total)
agg_cum[, scene_short := factor(scene_short, levels = finals$scene_short)]

p7 <- ggplot(agg_cum, aes(x = Day, y = cum_ch4, colour = scene_id)) +
  geom_line(linewidth = 0.9) +
  geom_text(data = finals,
            aes(x = 310, y = total, label = sprintf("%s  %.0f", scene_short, total),
                colour = scene_id),
            hjust = 0, size = 2.8, inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_manual(values = SCENE_COLS, labels = SCENE_SHORT, name = "Scenario") +
  scale_x_continuous(breaks = MON_BREAKS, labels = MON_LABELS,
                     expand = c(0.01, 0), limits = c(166, 360)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  labs(title    = "Cumulative Monsoon CH4 Emission by Scenario",
       subtitle = "Season-end totals (kg C/ha) labelled at right. Scenarios diverge with management intensity.",
       x = "Day of year", y = "Cumulative CH4 (kg C/ha)",
       caption = "Median across 263 cells × 20 years | Source: DNDC 9.5") +
  theme_ch4() + theme(legend.position = "none")

save_png(p7, "F7_monsoon_CH4_cumulative_ranked.png")

# =============================================================================
# M1 — Baseline monsoon CH4 tile map
# =============================================================================
message("\nM1: Baseline monsoon CH4 map ...")

m1 <- ggplot() +
  base_map_layers() +
  geom_tile(data = map_data[scene_id == "S000"],
            aes(x = x, y = y, fill = mean_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  ch4_fill_scale(name = "Mean monsoon\nCH4 (kg C/ha)") +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title    = "Mean Monsoon CH4 — Myanmar Rice Paddies",
       subtitle = "Baseline (S000) | 20-year mean (2005-2024) | DOY 166–309 | 0.5° resolution",
       caption  = "Source: DNDC 9.5 | SoilGrids | NASA POWER") +
  theme_map()

save_png(m1, "M1_monsoon_CH4_baseline.png", dir = MAP_DIR, w = 9, h = 12)

# =============================================================================
# M2 — All 13 scenarios small-multiple tile maps
# =============================================================================
message("M2: All-scenario monsoon maps ...")

ch4_lim <- range(map_data$mean_CH4, na.rm = TRUE)
m2_data <- map_data[scene_id %in% SCENE_ORDER]
m2_data[, scene_short := factor(SCENE_SHORT[scene_id], levels = SCENE_SHORT[SCENE_ORDER])]

m2 <- ggplot() +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = "grey93", colour = "grey55", linewidth = 0.25) +
  geom_tile(data = m2_data, aes(x = x, y = y, fill = mean_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  ch4_fill_scale(limits = ch4_lim) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.3) +
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_wrap(~ scene_short, ncol = 4) +
  scale_x_continuous(breaks = seq(94, 100, 3), labels = paste0(seq(94, 100, 3), "°E")) +
  scale_y_continuous(breaks = seq(12, 28,  4), labels = paste0(seq(12, 28,  4), "°N")) +
  labs(title    = "Mean Monsoon CH4 by Scenario — Myanmar Rice Paddies",
       subtitle = "20-year mean | DOY 166–309 | Same colour scale across all panels",
       caption  = "Source: DNDC 9.5 batch simulation") +
  theme_map(base = 8) +
  theme(legend.position = "right", strip.text = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 6))

save_png(m2, "M2_monsoon_CH4_all_scenarios.png", dir = MAP_DIR, w = 18, h = 20)

# =============================================================================
# M3 — Delta vs baseline (diverging)
# =============================================================================
message("M3: Monsoon delta vs baseline ...")

m3_data <- map_data[scene_id != "S000" & scene_id %in% SCENE_ORDER]
m3_data[, scene_short := factor(SCENE_SHORT[scene_id],
                                levels = SCENE_SHORT[SCENE_ORDER[SCENE_ORDER != "S000"]])]
dlim <- max(abs(m3_data$delta_CH4), na.rm = TRUE)

m3 <- ggplot() +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = "grey93", colour = "grey55", linewidth = 0.25) +
  geom_tile(data = m3_data, aes(x = x, y = y, fill = delta_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828",
                       midpoint = 0, limits = c(-dlim, dlim),
                       name = "ΔMonsoon CH4\n(kg C/ha)",
                       labels = label_number(accuracy = 1)) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.3) +
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_wrap(~ scene_short, ncol = 4) +
  scale_x_continuous(breaks = seq(94, 100, 3), labels = paste0(seq(94, 100, 3), "°E")) +
  scale_y_continuous(breaks = seq(12, 28,  4), labels = paste0(seq(12, 28,  4), "°N")) +
  labs(title    = "Monsoon CH4 Change vs Baseline by Scenario",
       subtitle = "Blue = reduction | Red = increase | 20-year mean | DOY 166–309",
       caption  = "Source: DNDC 9.5 | Myanmar irrigated rice") +
  theme_map(base = 8) +
  theme(legend.position = "right", strip.text = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 6))

save_png(m3, "M3_monsoon_CH4_delta_baseline.png", dir = MAP_DIR, w = 18, h = 20)

# =============================================================================
# M4 — AWD mitigation potential (S000 - S006)
# =============================================================================
message("M4: AWD mitigation ...")

awd <- merge(
  map_data[scene_id == "S000", .(cell_id, x, y, ch4_base = mean_CH4)],
  map_data[scene_id == "S006", .(cell_id, ch4_awd = mean_CH4)],
  by = "cell_id"
)
awd[, abs_red := ch4_base - ch4_awd]
awd[, pct_red := abs_red / ch4_base * 100]

m4_abs <- ggplot() + base_map_layers() +
  geom_tile(data = awd, aes(x = x, y = y, fill = abs_red),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(colours = c("#E3F2FD","#64B5F6","#1565C0","#0D47A1"),
                       name = "CH4 reduction\n(kg C/ha)",
                       labels = label_number(accuracy = 1)) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "AWD — Absolute CH4 Reduction", subtitle = "S000 − S006 | Monsoon") +
  theme_map()

m4_pct <- ggplot() + base_map_layers() +
  geom_tile(data = awd, aes(x = x, y = y, fill = pct_red),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(colours = c("#E8F5E9","#81C784","#2E7D32","#1B5E20"),
                       name = "CH4 reduction\n(%)",
                       labels = label_percent(scale = 1, accuracy = 1)) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "AWD — % CH4 Reduction", subtitle = "Darker = greater reduction") +
  theme_map()

m4 <- m4_abs + m4_pct +
  plot_annotation(
    title   = "AWD Mitigation Potential — Monsoon Season",
    caption = "Left: absolute (kg C/ha) | Right: % vs baseline | Source: DNDC 9.5",
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )

save_png(m4, "M4_AWD_mitigation_monsoon.png", dir = MAP_DIR, w = 18, h = 12)

# =============================================================================
# M5 — Coefficient of variation (interannual variability) — baseline
# =============================================================================
message("M5: Interannual CV map ...")

m5 <- ggplot() + base_map_layers() +
  geom_tile(data = map_data[scene_id == "S000"],
            aes(x = x, y = y, fill = cv_CH4),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(
    colours = c("#F3E5F5","#CE93D8","#7B1FA2","#4A148C"),
    name = "CV (%)\ninterannual",
    labels = label_number(accuracy = 1)
  ) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title    = "Interannual Variability of Monsoon CH4 — Baseline",
       subtitle = "Coefficient of variation (SD/mean × 100%) | 2005-2024 | Higher = more variable",
       caption  = "Source: DNDC 9.5 | Variability driven by climate variability") +
  theme_map()

save_png(m5, "M5_monsoon_CH4_CV_baseline.png", dir = MAP_DIR, w = 9, h = 12)

# =============================================================================
# M6 — Optimal (lowest-emitting) scenario per cell
# =============================================================================
message("M6: Optimal scenario map ...")

opt <- map_data[, .(
  best_scene = scene_id[which.min(mean_CH4)],
  best_CH4   = min(mean_CH4, na.rm = TRUE),
  base_CH4   = mean_CH4[scene_id == "S000"],
  max_red_pct= (mean_CH4[scene_id == "S000"] - min(mean_CH4, na.rm = TRUE)) /
    mean_CH4[scene_id == "S000"] * 100
), by = .(cell_id, x, y)]

opt[, scene_cat := fcase(
  best_scene == "S000", "baseline",
  best_scene %in% c("S001","S002","S003","S004"), "fertilizer",
  best_scene %in% c("S005","S006"), "water",
  best_scene %in% c("S008","S009"), "tillage",
  best_scene %in% c("S010","S011"), "irrigation",
  best_scene %in% c("S012","S013"), "residue",
  default = "baseline"
)]
opt[, scene_cat := factor(scene_cat, levels = names(CAT_COLS))]

m6_cat <- ggplot() + base_map_layers() +
  geom_tile(data = opt, aes(x = x, y = y, fill = scene_cat),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_manual(values = CAT_COLS, name = "Best category", drop = FALSE) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "Best Management Category", subtitle = "Category of lowest-emitting scenario per cell") +
  theme_map()

m6_pct <- ggplot() + base_map_layers() +
  geom_tile(data = opt, aes(x = x, y = y, fill = max_red_pct),
            width = TILE_SIZE, height = TILE_SIZE) +
  scale_fill_gradientn(
    colours = c("#FFF8E1","#FFCC02","#E65100","#B71C1C"),
    name = "Max reduction\n(%)", labels = label_number(accuracy = 1, suffix = "%")
  ) +
  geom_polygon(data = mmr_border, aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey30", linewidth = 0.5) +
  labs(title = "Maximum Achievable CH4 Reduction", subtitle = "% vs baseline | Best scenario per cell") +
  theme_map()

m6 <- m6_cat + m6_pct +
  plot_annotation(
    title   = "Optimal Management per Rice Cell — Monsoon Season",
    caption = "Left: category | Right: % reduction vs S000 | Source: DNDC 9.5",
    theme   = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.caption = element_text(size = 9, colour = "grey50"))
  )

save_png(m6, "M6_monsoon_optimal_scenario.png", dir = MAP_DIR, w = 18, h = 12)

# =============================================================================
# SUMMARY
# =============================================================================

all_png <- c(
  list.files(FIGURE_DIR, pattern = "\\.png$", full.names = TRUE),
  list.files(MAP_DIR,    pattern = "\\.png$", full.names = TRUE)
)

message("\n", strrep("=", 65))
message("COMPLETE — Monsoon figures saved")
message(strrep("=", 65))
message(sprintf("Run command    : Rscript 9_10_monsoon_figures.R"))
message(sprintf("Figures total  : %d PNG files\n", length(all_png)))
message("Flux figures -> ", FIGURE_DIR)
for (f in list.files(FIGURE_DIR, pattern = "\\.png$")) message(sprintf("  %s", f))
message("Spatial maps -> ", MAP_DIR)
for (f in list.files(MAP_DIR,    pattern = "\\.png$")) message(sprintf("  %s", f))
message(strrep("=", 65))

# ---- 13_monsoon_sequential_kgml.R ----
###############################################################################
# Script 13 — Monsoon-Only Sequential KGML Pipeline
#             LSTM + Transformer Encoder vs RF / XGBoost / KNN
#             Myanmar Rice Paddy CH4 Emission | Carbon Credit MRV
#
# Input   : C:/DNDC/simulate_dndc/dndc_outputs/enriched_daily_panel.rds
# Models  : LSTM, Transformer, Random Forest, XGBoost, KNN
# Season  : MONSOON ONLY (Day 166–309). Summer left entirely untouched.
#
# Changes in v1.1
# ─────────────────────────────────────────────────────────────────────────
#  - Added coro to package list (required by torch dataloader loops).
#  - Fixed flat feature matrix cbind: explicit as.matrix() on both sides;
#    replaced !is.finite() on data.frame with column-wise loop on matrix.
#  - Expanded STATIC_MGMT with full management factors + inline comments.
#
# Feature Policy
# ─────────────────────────────────────────────────────────────────────────
#  Sequence features (daily, observable):
#    Precipitation, Irrigation, Ponding
#
#  Static features (farm-observable or public data):
#    Soil   : clay, bulk_density, pH, SOC_static, porosity, sand, silt, elevation
#    Mgmt   : N fertiliser (total/urea/ammonium/n_apps/method),
#             flooding (total/FL1/FL2 days, AWD flag, irrigation code),
#             tillage (method code, n_passes),
#             residue (monsoon fraction, summer fraction)
#    Spatial: lat, lon
#
#  FORBIDDEN (runtime hard-stop — DNDC internal states/process variables):
#    CH4_flux/prod/oxid/pool, SOC/DOC/Microbe/Humads/Humus, NEE, NPP,
#    GrainC/LAI/TotalCropN/LeafC/StemC/RootC, Water_stress, N_stress,
#    Evaporation, Transpiration, IniSoilWater, EndSoilWater,
#    Leaching, Runoff, SoilHetResp, dSOC, LitterC, ManureC, DOC_leach
#
# Split Strategy
# ─────────────────────────────────────────────────────────────────────────
#  Spatial holdout : 20% cell_id → test (never seen during training)
#  Temporal split  : train 2005–2020 | val 2021–2024 (within train cells)
###############################################################################

cat("\n", strrep("=", 72), "\n")
cat("  Script 13 v1.1 — Monsoon Sequential KGML\n")
cat("  LSTM + Transformer vs RF / XGBoost / KNN\n")
cat(strrep("=", 72), "\n\n")

# =============================================================================
# 0. PACKAGES & CONFIG
# =============================================================================

pkgs_cran <- c("data.table","ggplot2","patchwork","scales",
               "ranger","xgboost","FNN","openxlsx","coro")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(sprintf("  Installing %s ...", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scales); library(ranger); library(xgboost)
  library(FNN); library(openxlsx); library(coro)
})

if (!requireNamespace("torch", quietly = TRUE))
  install.packages("torch", repos = "https://cloud.r-project.org")
if (!torch::torch_is_installed()) {
  message("  Installing torch backend (~500 MB, one-time only) ...")
  torch::install_torch()
}
library(torch)

set.seed(42)
torch::torch_manual_seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
PANEL_RDS <- "C:/DNDC/simulate_dndc/dndc_outputs/enriched_daily_panel.rds"
MODEL_DIR <- "C:/DNDC/simulate_dndc/models"
FIG_DIR   <- "C:/DNDC/simulate_dndc/figures/models"
OUT_DIR   <- "C:/DNDC/simulate_dndc/dndc_outputs"
for (d in c(MODEL_DIR, FIG_DIR, OUT_DIR)) dir.create(d, FALSE, TRUE)

# ── Season window ─────────────────────────────────────────────────────────────
MONSOON_DOY <- c(166L, 309L)
SEQ_LEN     <- MONSOON_DOY[2] - MONSOON_DOY[1] + 1L   # 144 days

# ── Split config ──────────────────────────────────────────────────────────────
TEST_CELL_FRAC <- 0.20
TRAIN_YEARS    <- 2005:2020
VAL_YEARS      <- 2021:2024

# ── LSTM hyperparameters ──────────────────────────────────────────────────────
LSTM_HIDDEN   <- 64L
LSTM_LAYERS   <- 2L
LSTM_DROPOUT  <- 0.25
LSTM_EPOCHS   <- 100L
LSTM_LR       <- 5e-4
LSTM_BATCH    <- 512L
LSTM_PATIENCE <- 15L
LSTM_CLIP     <- 1.0

# ── Transformer hyperparameters ───────────────────────────────────────────────
TF_D_MODEL  <- 64L      # must be divisible by TF_NHEAD
TF_NHEAD    <- 4L
TF_LAYERS   <- 2L
TF_FF_DIM   <- 256L
TF_DROPOUT  <- 0.20
TF_EPOCHS   <- 100L
TF_LR       <- 5e-4
TF_BATCH    <- 512L
TF_PATIENCE <- 15L
TF_CLIP     <- 1.0

device <- if (cuda_is_available()) "cuda" else "cpu"
message(sprintf("  Torch device: %s", device))


# =============================================================================
# 1. LOAD DAILY PANEL — MONSOON DAYS ONLY
# =============================================================================
cat("\n[1] Loading daily panel ...\n")

panel <- readRDS(PANEL_RDS)
setDT(panel)

.rename_safe <- function(dt, old, new) {
  present <- old %in% names(dt)
  if (any(present)) setnames(dt, old[present], new[present], skip_absent = TRUE)
}
.rename_safe(panel,
             old = c("CH4-flux","CH4-prod.","CH4-oxid.","CH4-pool",
                     "Litter-C","Manure-C","DOC-leach","Soil-heterotrophic-respiration"),
             new = c("CH4_flux","CH4_prod","CH4_oxid","CH4_pool",
                     "LitterC","ManureC","DOC_leach","SoilHetResp"))

message(sprintf("  Panel: %s rows × %d cols",
                format(nrow(panel), big.mark=","), ncol(panel)))
message("  Cols: ", paste(head(names(panel), 25), collapse=", "), " ...")

# Filter monsoon only — summer never enters any model object
monsoon_dt <- panel[Day >= MONSOON_DOY[1] & Day <= MONSOON_DOY[2]]
rm(panel); invisible(gc())

message(sprintf("  Monsoon subset: %s rows (Day %d–%d = %d days/year)",
                format(nrow(monsoon_dt), big.mark=","),
                MONSOON_DOY[1], MONSOON_DOY[2], SEQ_LEN))


# =============================================================================
# 2. FORBIDDEN VARIABLE GUARD
# =============================================================================
cat("\n[2] Enforcing feature policy ...\n")

FORBIDDEN_VARS <- c(
  # Target + CH4 process intermediates (direct leakage)
  "CH4_flux","CH4_prod","CH4_oxid","CH4_pool",
  # Carbon pool dynamics — DNDC internal states
  "SOC","dSOC","DOC","Microbe","Humads","Humus",
  "LitterC","ManureC","DOC_leach","SoilHetResp",
  # Ecosystem fluxes — DNDC computed
  "NEE","NPP",
  # Crop physiological state — simulated by DNDC crop growth module
  "GrainC","LAI","TotalCropN","LeafC","StemC","RootC",
  "Water_stress","N_stress",
  # Internal water balance — DNDC computed fluxes
  "Evaporation","Transpiration","IniSoilWater","EndSoilWater",
  "Leaching","Runoff"
)

.enforce_feature_policy <- function(feature_vec, context="") {
  hits <- intersect(feature_vec, FORBIDDEN_VARS)
  if (length(hits) > 0)
    stop(sprintf(
      "[FORBIDDEN FEATURE — %s] DNDC internal variables not permitted:\n  %s\n",
      context, paste(hits, collapse=", ")))
  invisible(TRUE)
}

message(sprintf("  %d forbidden variables registered.", length(FORBIDDEN_VARS)))


# =============================================================================
# 3. COLUMN RENAMING & FEATURE DECLARATION
# =============================================================================
cat("\n[3] Declaring features ...\n")

# Soil columns arrive hyphenated from cells.rds
SOIL_RENAME_MAP <- c(
  "clay_0-5cm"="clay_fraction", "bdod_0-5cm"="bulk_density",
  "pH_0-5cm"="soil_pH",         "soc_0-5cm"="SOC_static",
  "poros_0-5cm"="porosity",     "sand_0-5cm"="sand_fraction",
  "silt_0-5cm"="silt_fraction", "elevation"="elevation_m"
)
present_soil <- names(SOIL_RENAME_MAP)[names(SOIL_RENAME_MAP) %in% names(monsoon_dt)]
if (length(present_soil) > 0)
  setnames(monsoon_dt, present_soil, SOIL_RENAME_MAP[present_soil])

if (all(c("x","y") %in% names(monsoon_dt)) && !"lon" %in% names(monsoon_dt))
  setnames(monsoon_dt, c("x","y"), c("lon","lat"))

# ── Time-varying sequence features (daily, observable) ────────────────────────
SEQ_FEATURES <- c(
  "Precipitation",   # mm/day — NASA POWER rainfall (drives DNDC water balance)
  "Irrigation",      # mm/day — management schedule input from scene_registry
  "Ponding"          # mm     — water depth; set by flooding management + rainfall
)
.enforce_feature_policy(SEQ_FEATURES, "SEQ_FEATURES")

# ── Static soil features (SoilGrids 0–5 cm, public data) ─────────────────────
STATIC_SOIL <- c(
  "clay_fraction",   # g/g   — texture; controls CH4 diffusion and drainage
  "bulk_density",    # g/cm³ — affects porosity and aeration status
  "soil_pH",         # pH    — controls methanogen community activity
  "SOC_static",      # kg/kg — initial organic C substrate for methanogenesis
  "porosity",        # v/v   — saturated zone capacity; flooding persistence
  "sand_fraction",   # g/g   — drainage / water retention proxy
  "silt_fraction",   # g/g   — texture complement
  "elevation_m"      # m     — altitude proxy for temperature + drainage
)

# ── Static management features (from scene_registry; farm-observable) ─────────
STATIC_MGMT <- c(
  # ── Fertiliser ─────────────────────────────────────────────────────────────
  "fert_total_N_kgNha",    # kg N/ha  total N applied across all events
  "fert_urea_kgNha",       # kg N/ha  urea fraction of total N
  "fert_ammonium_kgNha",   # kg N/ha  ammonium sulphate fraction
  "fert_n_applications",   # integer  number of split application events
  "fert_method_code",      # code     0=broadcast, 1=incorporated, 2=injected
  # ── Water management ───────────────────────────────────────────────────────
  "total_flood_days",      # days     FL1 + FL2 combined flooding duration
  "FL1_flood_days",        # days     first flood event duration
  "FL2_flood_days",        # days     second flood event duration (if present)
  "awd_flag",              # binary   0=continuous flooding, 1=AWD
  "irrigation_control_code", # code   0=rainfed, 1=controlled, 2=continuous
  # ── Tillage ────────────────────────────────────────────────────────────────
  "till_method_code",      # code     0=no-till, 1=conventional, 2=deep plow
  "till_applications",     # integer  number of tillage passes per season
  # ── Residue management ─────────────────────────────────────────────────────
  "residue_monsoon_frac",  # fraction monsoon crop residue returned to soil
  "residue_summer_frac"    # fraction summer crop residue returned to soil
)

# ── Spatial context ───────────────────────────────────────────────────────────
STATIC_SPATIAL <- c(
  "lat",   # degrees N — temperature / solar radiation gradient
  "lon"    # degrees E — regional precipitation pattern proxy
)

STATIC_FEATURES <- c(STATIC_SOIL, STATIC_MGMT, STATIC_SPATIAL)
.enforce_feature_policy(STATIC_FEATURES, "STATIC_FEATURES")

TARGET <- "ann_CH4_kgCha"   # seasonal monsoon CH4 sum (kg C/ha)

# Resolve which declared columns are actually present
avail_seq    <- intersect(SEQ_FEATURES,    names(monsoon_dt))
avail_static <- intersect(STATIC_FEATURES, names(monsoon_dt))

if (length(avail_seq) == 0)
  stop("No sequence features found. Check Precipitation/Irrigation/Ponding in panel.")
if (!"CH4_flux" %in% names(monsoon_dt))
  stop("CH4_flux missing — run Script 8 to regenerate enriched_daily_panel.rds.")

message(sprintf("  Seq features   (%d/%d): %s",
                length(avail_seq), length(SEQ_FEATURES),
                paste(avail_seq, collapse=", ")))
message(sprintf("  Static features(%d/%d) — soil:%d  mgmt:%d  spatial:%d",
                length(avail_static), length(STATIC_FEATURES),
                sum(avail_static %in% STATIC_SOIL),
                sum(avail_static %in% STATIC_MGMT),
                sum(avail_static %in% STATIC_SPATIAL)))

missing_static <- setdiff(STATIC_FEATURES, names(monsoon_dt))
if (length(missing_static) > 0)
  message("  Not found (skipped): ", paste(missing_static, collapse=", "))


# =============================================================================
# 4. AGGREGATE TARGET + STATIC FEATURES
# =============================================================================
cat("\n[4] Aggregating target and static features ...\n")

monsoon_dt[, monsoon_day := Day - MONSOON_DOY[1] + 1L]

sample_agg <- monsoon_dt[, {
  ch4_sum     <- sum(CH4_flux, na.rm=TRUE)
  static_list <- lapply(avail_static, function(col) {
    v <- .SD[[col]]; if (length(v) > 0) v[[1L]] else NA_real_
  })
  names(static_list) <- avail_static
  c(list(ann_CH4_kgCha=ch4_sum), static_list)
}, by=.(cell_id, scene_id, cal_year),
.SDcols=c("CH4_flux", avail_static)]

setorder(sample_agg, cell_id, scene_id, cal_year)
sample_agg[, sample_idx := .I]
N_samples <- nrow(sample_agg)

message(sprintf("  Samples: %s (%d cells × %d scenes × %d years)",
                format(N_samples, big.mark=","),
                uniqueN(sample_agg$cell_id),
                uniqueN(sample_agg$scene_id),
                uniqueN(sample_agg$cal_year)))

monsoon_dt <- monsoon_dt[
  sample_agg[, .(cell_id, scene_id, cal_year, sample_idx)],
  on=.(cell_id, scene_id, cal_year)
]


# =============================================================================
# 5. BUILD 3D SEQUENCE ARRAY  [N, SEQ_LEN, n_seq_feat]
# =============================================================================
cat(sprintf("\n[5] Building sequence array [%d × %d × %d] ...\n",
            N_samples, SEQ_LEN, length(avail_seq)))

setorder(monsoon_dt, sample_idx, monsoon_day)

rows_check <- monsoon_dt[, .N, by=sample_idx]
n_short    <- sum(rows_check$N != SEQ_LEN)
if (n_short > 0)
  message(sprintf("  WARNING: %d samples have != %d days (padding with 0).",
                  n_short, SEQ_LEN))

n_seq_feat <- length(avail_seq)
X_seq_raw  <- array(0.0, dim=c(N_samples, SEQ_LEN, n_seq_feat),
                    dimnames=list(NULL, NULL, avail_seq))

n_total <- N_samples * SEQ_LEN
for (fi in seq_len(n_seq_feat)) {
  col      <- avail_seq[fi]
  raw_vals <- as.numeric(monsoon_dt[[col]])
  raw_vals[!is.finite(raw_vals)] <- 0.0
  n_have   <- length(raw_vals)
  if (n_have < n_total) raw_vals <- c(raw_vals, rep(0.0, n_total - n_have))
  X_seq_raw[,,fi] <- matrix(raw_vals[seq_len(n_total)],
                            nrow=N_samples, ncol=SEQ_LEN, byrow=TRUE)
}
message(sprintf("  Array: %.1f MB", object.size(X_seq_raw) / 1024^2))

# Static matrix [N, n_static] — explicitly numeric matrix
n_static_feat <- length(avail_static)
X_static_raw  <- matrix(0.0, nrow=N_samples, ncol=n_static_feat,
                        dimnames=list(NULL, avail_static))
for (col in avail_static) {
  vals <- as.numeric(sample_agg[[col]])
  vals[!is.finite(vals)] <- 0.0
  X_static_raw[, col] <- vals
}
message(sprintf("  Static matrix: %d × %d", N_samples, n_static_feat))

y_raw <- as.numeric(sample_agg[[TARGET]])
y_raw[!is.finite(y_raw)] <- 0.0


# =============================================================================
# 6. TRAIN / VAL / TEST SPLIT
# =============================================================================
cat("\n[6] Splitting data ...\n")

all_cells   <- sort(unique(sample_agg$cell_id))
n_test      <- max(1L, round(length(all_cells) * TEST_CELL_FRAC))
set.seed(42)
test_cells  <- sort(sample(all_cells, n_test))
train_cells <- setdiff(all_cells, test_cells)

idx_trainval <- which(sample_agg$cell_id %in% train_cells)
idx_test     <- which(sample_agg$cell_id %in% test_cells)
idx_train    <- idx_trainval[sample_agg$cal_year[idx_trainval] %in% TRAIN_YEARS]
idx_val      <- idx_trainval[sample_agg$cal_year[idx_trainval] %in% VAL_YEARS]

stopifnot(!any(sample_agg$cell_id[idx_train] %in% test_cells))
stopifnot(!any(sample_agg$cell_id[idx_val]   %in% test_cells))
stopifnot(length(intersect(idx_train, idx_val)) == 0L)

message(sprintf("  Cells  : total=%d | train=%d | test=%d (%.0f%% holdout)",
                length(all_cells), length(train_cells), length(test_cells),
                TEST_CELL_FRAC*100))
message(sprintf("  Samples: train=%d | val=%d | test=%d",
                length(idx_train), length(idx_val), length(idx_test)))


# =============================================================================
# 7. FEATURE SCALING  (fit ONLY on training data)
# =============================================================================
cat("\n[7] Fitting scalers on training set ...\n")

# Target: log1p + z-score
y_log      <- log1p(pmax(0.0, y_raw))
y_log_mean <- mean(y_log[idx_train])
y_log_sd   <- max(sd(y_log[idx_train]), 1e-8)
y_scaled   <- (y_log - y_log_mean) / y_log_sd
back_transform <- function(y_sc) expm1(y_sc * y_log_sd + y_log_mean)

# Sequence scaler (per-feature standardisation over training time steps)
seq_means <- vapply(seq_len(n_seq_feat),
                    function(fi) mean(X_seq_raw[idx_train,,fi]), numeric(1))
seq_sds   <- vapply(seq_len(n_seq_feat),
                    function(fi) max(sd(as.vector(X_seq_raw[idx_train,,fi])), 1e-8), numeric(1))
names(seq_means) <- names(seq_sds) <- avail_seq

.scale_seq <- function(arr, m, s) {
  for (fi in seq_along(m)) arr[,,fi] <- (arr[,,fi] - m[fi]) / s[fi]
  arr
}
X_seq_sc <- .scale_seq(X_seq_raw, seq_means, seq_sds)

# Static scaler
static_means <- colMeans(X_static_raw[idx_train,, drop=FALSE], na.rm=TRUE)
static_sds   <- pmax(apply(X_static_raw[idx_train,, drop=FALSE], 2, sd, na.rm=TRUE), 1e-8)
.scale_static <- function(mat, m, s) sweep(sweep(mat, 2, m, "-"), 2, s, "/")
X_static_sc <- .scale_static(X_static_raw, static_means, static_sds)
X_static_sc[!is.finite(X_static_sc)] <- 0.0

message("  Seq means: ", paste(names(seq_means),"=",round(seq_means,2), collapse=" | "))
message("  Seq sds  : ", paste(round(seq_sds,2), collapse=" | "))


# =============================================================================
# 8. FLAT FEATURE MATRIX  (for RF / XGBoost / KNN)
# =============================================================================
cat("\n[8] Building flat feature matrix for ML baselines ...\n")

# Summarise each daily sequence into 6 statistics
flat_blocks <- lapply(seq_len(n_seq_feat), function(fi) {
  mat  <- X_seq_raw[,,fi]          # [N, T] raw values
  feat <- avail_seq[fi]
  out  <- cbind(
    rowSums(mat,  na.rm=TRUE),
    rowMeans(mat, na.rm=TRUE),
    apply(mat, 1, max,      na.rm=TRUE),
    apply(mat, 1, sd,       na.rm=TRUE),
    apply(mat, 1, quantile, probs=0.75, na.rm=TRUE),
    mat[, SEQ_LEN]
  )
  colnames(out) <- paste0(feat, c("_sum","_mean","_max","_sd","_q75","_last"))
  out
})
# Both sides are explicitly matrices — no list coercion
flat_seq_mat <- do.call(cbind, flat_blocks)               # [N, 6 × n_seq_feat]
flat_all_raw <- cbind(flat_seq_mat, X_static_raw)         # [N, seq_stats + static]

# Replace non-finite values column by column (avoids type error on data.frames)
for (j in seq_len(ncol(flat_all_raw)))
  flat_all_raw[!is.finite(flat_all_raw[, j]), j] <- 0.0

.enforce_feature_policy(colnames(flat_all_raw), "ML flat feature matrix")

flat_means <- colMeans(flat_all_raw[idx_train,, drop=FALSE])
flat_sds   <- pmax(apply(flat_all_raw[idx_train,, drop=FALSE], 2, sd), 1e-8)
flat_sc    <- sweep(sweep(flat_all_raw, 2, flat_means, "-"), 2, flat_sds, "/")
for (j in seq_len(ncol(flat_sc)))
  flat_sc[!is.finite(flat_sc[, j]), j] <- 0.0

X_flat_train <- flat_sc[idx_train,, drop=FALSE]
X_flat_val   <- flat_sc[idx_val,,   drop=FALSE]
X_flat_test  <- flat_sc[idx_test,,  drop=FALSE]
y_train_raw  <- y_raw[idx_train]
y_val_raw    <- y_raw[idx_val]
y_test_raw   <- y_raw[idx_test]

n_feat_flat <- ncol(flat_sc)
message(sprintf("  Flat features: %d (%d seq-stats + %d static) | train:%d val:%d test:%d",
                n_feat_flat, ncol(flat_seq_mat), n_static_feat,
                length(idx_train), length(idx_val), length(idx_test)))


# =============================================================================
# 9. EVALUATION HELPER
# =============================================================================

evaluate_model <- function(y_actual, y_pred, model_name) {
  ya <- as.numeric(y_actual); yp <- as.numeric(y_pred)
  ok <- is.finite(ya) & is.finite(yp)
  ya <- ya[ok]; yp <- yp[ok]; n <- length(ya)
  ss_res <- sum((ya-yp)^2); ss_tot <- sum((ya-mean(ya))^2)
  r2    <- if (ss_tot < 1e-12) NA_real_ else 1 - ss_res/ss_tot
  rmse  <- sqrt(mean((ya-yp)^2))
  mae   <- mean(abs(ya-yp))
  mbe   <- mean(yp-ya)
  nrmse <- if (mean(ya) > 1e-6) rmse/mean(ya)*100 else NA_real_
  rpd   <- sd(ya) / max(rmse, 1e-12)
  data.table(Model=model_name,
             R2=round(r2,4), RMSE=round(rmse,3), MAE=round(mae,3),
             MBE=round(mbe,3), NRMSE_pct=round(nrmse,2),
             RPD=round(rpd,3), N=n)
}

results_list <- list()
pred_list    <- list()


# =============================================================================
# 10. TORCH DATASET & DATALOADERS
# =============================================================================
cat("\n[10] Building torch datasets ...\n")

.enforce_feature_policy(avail_seq,    "DL seq features")
.enforce_feature_policy(avail_static, "DL static features")

ch4_seq_dataset <- dataset(
  name = "CH4SeqDataset",
  initialize = function(seq_arr, static_mat, y_vec) {
    self$x_seq    <- torch_tensor(seq_arr,    dtype=torch_float())$to(device=device)
    self$x_static <- torch_tensor(static_mat, dtype=torch_float())$to(device=device)
    self$y        <- torch_tensor(y_vec,       dtype=torch_float())$to(device=device)
  },
  .getitem = function(i) list(
    x_seq    = self$x_seq[i,,],
    x_static = self$x_static[i,],
    y        = self$y[i]$unsqueeze(1L)
  ),
  .length = function() self$x_seq$size(1L)
)

.mk_ds <- function(idx)
  ch4_seq_dataset(
    seq_arr    = X_seq_sc[idx,,,  drop=FALSE],
    static_mat = X_static_sc[idx,, drop=FALSE],
    y_vec      = y_scaled[idx]
  )

train_ds <- .mk_ds(idx_train)
val_ds   <- .mk_ds(idx_val)
test_ds  <- .mk_ds(idx_test)

.mk_dl <- function(ds, batch, shuffle)
  dataloader(ds, batch_size=batch, shuffle=shuffle, drop_last=FALSE)


# =============================================================================
# 11. DL UTILITIES
# =============================================================================

# Sinusoidal PE → R matrix [T, d_model]
.make_sinusoidal_pe <- function(T, d_model) {
  pe <- matrix(0.0, nrow=T, ncol=d_model)
  for (pos in seq_len(T))
    for (i in seq(1, d_model-1, 2)) {
      freq       <- (pos-1) / (10000^((i-1L)/d_model))
      pe[pos, i] <- sin(freq)
      if (i+1L <= d_model) pe[pos, i+1L] <- cos(freq)
    }
  pe
}

.train_epoch <- function(model, dl, optimizer, clip) {
  model$train(); total <- 0.0; n <- 0L
  coro::loop(for (batch in dl) {
    optimizer$zero_grad()
    pred <- model(batch$x_seq, batch$x_static)
    loss <- nnf_mse_loss(pred, batch$y)
    loss$backward()
    nn_utils_clip_grad_norm_(model$parameters, max_norm=clip)
    optimizer$step()
    total <- total + loss$item(); n <- n + 1L
  })
  total / max(n, 1L)
}

.eval_epoch <- function(model, dl) {
  model$eval(); total <- 0.0; n <- 0L
  with_no_grad({
    coro::loop(for (batch in dl) {
      loss  <- nnf_mse_loss(model(batch$x_seq, batch$x_static), batch$y)
      total <- total + loss$item(); n <- n + 1L
    })
  })
  total / max(n, 1L)
}

.get_dl_preds <- function(model, dl) {
  model$eval(); preds <- numeric(0)
  with_no_grad({
    coro::loop(for (batch in dl) {
      out   <- model(batch$x_seq, batch$x_static)
      preds <- c(preds, as.numeric(out$squeeze(2L)$cpu()))
    })
  })
  preds
}

.fit_dl_model <- function(model, train_dl, val_dl, optimizer, scheduler,
                          n_epochs, patience, clip, model_name) {
  best_val <- Inf; pat_ctr <- 0L; best_state <- NULL
  t_hist   <- numeric(n_epochs); v_hist <- numeric(n_epochs); ep_run <- 0L
  
  for (ep in seq_len(n_epochs)) {
    t_loss <- .train_epoch(model, train_dl, optimizer, clip)
    v_loss <- .eval_epoch(model, val_dl)
    scheduler$step()
    t_hist[ep] <- t_loss; v_hist[ep] <- v_loss; ep_run <- ep
    
    if (ep %% 10L == 0L || ep == 1L)
      message(sprintf("    [%s] Ep%3d | train=%.5f  val=%.5f",
                      model_name, ep, t_loss, v_loss))
    
    if (v_loss < best_val - 1e-6) {
      best_val <- v_loss; best_state <- model$state_dict(); pat_ctr <- 0L
    } else {
      pat_ctr <- pat_ctr + 1L
      if (pat_ctr >= patience) {
        message(sprintf("    [%s] Early stop ep%d", model_name, ep)); break
      }
    }
  }
  if (!is.null(best_state)) model$load_state_dict(best_state)
  message(sprintf("    [%s] Best val=%.5f", model_name, best_val))
  list(model=model, train_hist=t_hist[seq_len(ep_run)],
       val_hist=v_hist[seq_len(ep_run)], n_epochs=ep_run)
}


# =============================================================================
# 12. LSTM MODEL
# =============================================================================
cat("\n[12] === LSTM ===\n")

lstm_net <- nn_module(
  "CH4_LSTM",
  initialize = function(n_seq, n_static, hidden, n_layers, dropout) {
    self$lstm <- nn_lstm(
      input_size=n_seq, hidden_size=hidden, num_layers=n_layers,
      batch_first=TRUE, dropout=if(n_layers>1L) dropout else 0.0
    )
    self$seq_norm  <- nn_layer_norm(hidden)
    self$seq_drop  <- nn_dropout(dropout)
    self$static_fc <- nn_sequential(nn_linear(n_static,32L), nn_relu(), nn_dropout(dropout))
    self$head      <- nn_sequential(nn_linear(hidden+32L,64L), nn_relu(),
                                    nn_dropout(dropout), nn_linear(64L,1L))
  },
  forward = function(x_seq, x_static) {
    out    <- self$lstm(x_seq)
    h_last <- out[[1]][,-1L,]         # last time step [B, hidden]
    h_last <- self$seq_drop(self$seq_norm(h_last))
    s      <- self$static_fc(x_static)
    self$head(torch_cat(list(h_last, s), dim=2L))
  }
)

lstm_model <- lstm_net(n_seq_feat, n_static_feat,
                       LSTM_HIDDEN, LSTM_LAYERS, LSTM_DROPOUT)$to(device=device)
lstm_opt   <- optim_adam(lstm_model$parameters, lr=LSTM_LR, weight_decay=1e-5)
lstm_sch   <- lr_step(lstm_opt, step_size=20L, gamma=0.5)

train_dl_lstm <- .mk_dl(train_ds, LSTM_BATCH, TRUE)
val_dl_lstm   <- .mk_dl(val_ds,   LSTM_BATCH, FALSE)
test_dl_lstm  <- .mk_dl(test_ds,  LSTM_BATCH, FALSE)

t0_lstm  <- proc.time()
lstm_fit <- .fit_dl_model(lstm_model, train_dl_lstm, val_dl_lstm,
                          lstm_opt, lstm_sch,
                          LSTM_EPOCHS, LSTM_PATIENCE, LSTM_CLIP, "LSTM")
elapsed_lstm <- (proc.time()-t0_lstm)[3]
message(sprintf("  LSTM: %.1f s | %d epochs", elapsed_lstm, lstm_fit$n_epochs))

lstm_pred_test  <- back_transform(.get_dl_preds(lstm_fit$model, test_dl_lstm))
lstm_pred_train <- back_transform(.get_dl_preds(lstm_fit$model,
                                                .mk_dl(train_ds, LSTM_BATCH, FALSE)))
results_list[["LSTM_train"]] <- evaluate_model(y_train_raw, lstm_pred_train, "LSTM (train)")
results_list[["LSTM_test"]]  <- evaluate_model(y_test_raw,  lstm_pred_test,  "LSTM (test)")
pred_list[["LSTM"]] <- data.table(actual=y_test_raw, predicted=lstm_pred_test,
                                  model="LSTM", residual=lstm_pred_test-y_test_raw)
torch_save(lstm_fit$model$state_dict(), file.path(MODEL_DIR,"lstm_monsoon_state.pt"))
message("  Saved: lstm_monsoon_state.pt")


# =============================================================================
# 13. TRANSFORMER ENCODER MODEL
# =============================================================================
cat("\n[13] === Transformer Encoder ===\n")

tf_has_batch_first <- tryCatch({
  tl <- nn_transformer_encoder_layer(d_model=8L, nhead=2L, batch_first=TRUE)
  rm(tl); TRUE
}, error=function(e) FALSE)
message(sprintf("  batch_first=%s (torch >= 0.10 required)", tf_has_batch_first))

pe_mat <- .make_sinusoidal_pe(SEQ_LEN, TF_D_MODEL)   # [T, d_model]

transformer_net <- nn_module(
  "CH4_Transformer",
  initialize = function(n_seq, n_static, d_model, nhead, n_layers,
                        ff_dim, dropout, pe_matrix) {
    self$input_proj <- nn_linear(n_seq, d_model)
    # Store PE as plain [T, d_model] tensor; unsqueeze to [1, T, d_model] in
    # forward() so it broadcasts correctly against [B, T, d_model].
    # Pre-permuting at init time produces wrong strides that cause a
    # "size of tensor a (B) must match b (T)" error at runtime.
    self$register_buffer(
      "pe",
      torch_tensor(pe_matrix, dtype=torch_float())   # [T, d_model]
    )
    enc_layer <- nn_transformer_encoder_layer(
      d_model=d_model, nhead=nhead, dim_feedforward=ff_dim,
      dropout=dropout, batch_first=TRUE
    )
    self$encoder   <- nn_transformer_encoder(enc_layer, num_layers=n_layers)
    self$static_fc <- nn_sequential(nn_linear(n_static,32L), nn_relu(), nn_dropout(dropout))
    self$pool_drop <- nn_dropout(dropout)
    self$head      <- nn_sequential(nn_linear(d_model+32L,64L), nn_relu(),
                                    nn_dropout(dropout), nn_linear(64L,1L))
  },
  forward = function(x_seq, x_static) {
    # x_seq   : [B, T, F_seq]
    # self$pe : [T, d_model] -> unsqueeze(1L) -> [1, T, d_model]
    #           broadcasts against [B, T, d_model] automatically
    h_proj <- self$input_proj(x_seq)                         # [B, T, d_model]
    h_enc  <- self$encoder(h_proj + self$pe$unsqueeze(1L))   # [B, T, d_model]
    h_pool <- self$pool_drop(h_enc$mean(dim=2L))             # [B, d_model]
    self$head(torch_cat(list(h_pool, self$static_fc(x_static)), dim=2L))
  }
)

tf_model <- transformer_net(n_seq_feat, n_static_feat,
                            TF_D_MODEL, TF_NHEAD, TF_LAYERS, TF_FF_DIM,
                            TF_DROPOUT, pe_matrix=pe_mat)$to(device=device)
tf_opt   <- optim_adam(tf_model$parameters, lr=TF_LR, weight_decay=1e-5)
tf_sch   <- lr_step(tf_opt, step_size=20L, gamma=0.5)

train_dl_tf <- .mk_dl(train_ds, TF_BATCH, TRUE)
val_dl_tf   <- .mk_dl(val_ds,   TF_BATCH, FALSE)
test_dl_tf  <- .mk_dl(test_ds,  TF_BATCH, FALSE)

t0_tf  <- proc.time()
tf_fit <- .fit_dl_model(tf_model, train_dl_tf, val_dl_tf,
                        tf_opt, tf_sch,
                        TF_EPOCHS, TF_PATIENCE, TF_CLIP, "Transformer")
elapsed_tf <- (proc.time()-t0_tf)[3]
message(sprintf("  Transformer: %.1f s | %d epochs", elapsed_tf, tf_fit$n_epochs))

tf_pred_test  <- back_transform(.get_dl_preds(tf_fit$model, test_dl_tf))
tf_pred_train <- back_transform(.get_dl_preds(tf_fit$model,
                                              .mk_dl(train_ds, TF_BATCH, FALSE)))
results_list[["TF_train"]] <- evaluate_model(y_train_raw, tf_pred_train, "Transformer (train)")
results_list[["TF_test"]]  <- evaluate_model(y_test_raw,  tf_pred_test,  "Transformer (test)")
pred_list[["Transformer"]] <- data.table(actual=y_test_raw, predicted=tf_pred_test,
                                         model="Transformer",
                                         residual=tf_pred_test-y_test_raw)
torch_save(tf_fit$model$state_dict(), file.path(MODEL_DIR,"transformer_monsoon_state.pt"))
message("  Saved: transformer_monsoon_state.pt")


# =============================================================================
# 14. RANDOM FOREST (ranger)
# =============================================================================
cat("\n[14] === Random Forest ===\n")
t0 <- proc.time()

rf_model <- ranger::ranger(
  x=X_flat_train, y=y_train_raw,
  num.trees=500L, mtry=max(1L, floor(sqrt(n_feat_flat))),
  min.node.size=5L, importance="impurity", seed=42L, num.threads=4L
)
rf_pred_train <- rf_model$predictions
rf_pred_test  <- predict(rf_model, data=X_flat_test)$predictions
message(sprintf("  %.1f s | OOB R2=%.4f", (proc.time()-t0)[3], rf_model$r.squared))

results_list[["RF_train"]] <- evaluate_model(y_train_raw, rf_pred_train, "RF (train OOB)")
results_list[["RF_test"]]  <- evaluate_model(y_test_raw,  rf_pred_test,  "RF (test)")
pred_list[["RF"]] <- data.table(actual=y_test_raw, predicted=rf_pred_test,
                                model="Random Forest", residual=rf_pred_test-y_test_raw)
rf_imp <- sort(rf_model$variable.importance, decreasing=TRUE)
saveRDS(rf_model, file.path(MODEL_DIR,"rf_monsoon.rds"))
message("  Saved: rf_monsoon.rds")


# =============================================================================
# 15. GRADIENT BOOSTING (XGBoost)
# =============================================================================
cat("\n[15] === Gradient Boosting (XGBoost) ===\n")
t0 <- proc.time()

xgb_dm_train <- xgb.DMatrix(data=X_flat_train, label=y_train_raw)
xgb_dm_val   <- xgb.DMatrix(data=X_flat_val,   label=y_val_raw)
xgb_dm_test  <- xgb.DMatrix(data=X_flat_test,  label=y_test_raw)

xgb_params <- list(
  booster="gbtree", objective="reg:squarederror",
  eta=0.05, max_depth=6L, subsample=0.80, colsample_bytree=0.80,
  min_child_weight=5L, lambda=1.0, alpha=0.1, eval_metric="rmse"
)

xgb_cv <- xgb.cv(params=xgb_params, data=xgb_dm_train,
                 nrounds=600L, nfold=5L, early_stopping_rounds=30L,
                 verbose=1L, print_every_n=100L, seed=42L)

best_rounds <- xgb_cv$best_iteration
if (is.null(best_rounds) || is.na(best_rounds) || best_rounds < 1L)
  best_rounds <- which.min(xgb_cv$evaluation_log$test_rmse_mean)
best_rounds <- as.integer(best_rounds)
message(sprintf("  Best rounds: %d", best_rounds))

xgb_model <- xgb.train(
  params=xgb_params, data=xgb_dm_train, nrounds=best_rounds,
  watchlist=list(train=xgb_dm_train, val=xgb_dm_val),
  print_every_n=100L, verbose=1L
)
xgb_pred_train <- predict(xgb_model, xgb_dm_train)
xgb_pred_test  <- predict(xgb_model, xgb_dm_test)
message(sprintf("  %.1f s", (proc.time()-t0)[3]))

results_list[["GBM_train"]] <- evaluate_model(y_train_raw, xgb_pred_train, "GBM (train CV)")
results_list[["GBM_test"]]  <- evaluate_model(y_test_raw,  xgb_pred_test,  "GBM (test)")
pred_list[["GBM"]] <- data.table(actual=y_test_raw, predicted=xgb_pred_test,
                                 model="Gradient Boosting",
                                 residual=xgb_pred_test-y_test_raw)
xgb_imp_raw <- xgb.importance(feature_names=colnames(X_flat_train), model=xgb_model)
xgb.save(xgb_model, file.path(MODEL_DIR,"xgb_monsoon.bin"))
message("  Saved: xgb_monsoon.bin")


# =============================================================================
# 16. K-NEAREST NEIGHBOURS (FNN)
# =============================================================================
cat("\n[16] === KNN ===\n")
t0 <- proc.time()

k_cands  <- c(3L,5L,7L,10L,15L,20L)
fold_ids <- sample(rep(1:5, length.out=nrow(X_flat_train)))
knn_cv_rmse <- vapply(k_cands, function(k) {
  mean(vapply(1:5, function(f) {
    Xtr <- X_flat_train[fold_ids!=f,, drop=FALSE]
    Xva <- X_flat_train[fold_ids==f,, drop=FALSE]
    ytr <- y_train_raw[fold_ids!=f]; yva <- y_train_raw[fold_ids==f]
    sqrt(mean((yva - FNN::knn.reg(Xtr,Xva,ytr,k)$pred)^2))
  }, numeric(1)))
}, numeric(1))

best_k <- k_cands[which.min(knn_cv_rmse)]
message(sprintf("  Best k=%d (CV RMSE=%.4f)", best_k, min(knn_cv_rmse)))

knn_pred_train <- FNN::knn.reg(X_flat_train, X_flat_train, y_train_raw, best_k)$pred
knn_pred_test  <- FNN::knn.reg(X_flat_train, X_flat_test,  y_train_raw, best_k)$pred
message(sprintf("  %.1f s", (proc.time()-t0)[3]))

knn_label <- sprintf("KNN k=%d", best_k)
results_list[["KNN_train"]] <- evaluate_model(y_train_raw, knn_pred_train,
                                              paste0(knn_label," (train)"))
results_list[["KNN_test"]]  <- evaluate_model(y_test_raw,  knn_pred_test,
                                              paste0(knn_label," (test)"))
pred_list[["KNN"]] <- data.table(actual=y_test_raw, predicted=knn_pred_test,
                                 model=knn_label, residual=knn_pred_test-y_test_raw)
saveRDS(list(X_train=X_flat_train, y_train=y_train_raw, k=best_k,
             flat_means=flat_means, flat_sds=flat_sds),
        file.path(MODEL_DIR,"knn_monsoon.rds"))
message("  Saved: knn_monsoon.rds")


# =============================================================================
# 17. COMPILE & RANK RESULTS
# =============================================================================
cat("\n[17] Compiling results ...\n")

results_all  <- rbindlist(results_list)
results_all[, Split     := fifelse(grepl("train|OOB|CV", Model, ignore.case=TRUE),
                                   "Train/OOB", "Test")]
results_all[, ModelName := sub(" \\(.*","", Model)]
test_results <- results_all[Split=="Test"][order(-R2)]

message("\nTEST SET PERFORMANCE:")
print(test_results[, .(Model, R2, RMSE, MAE, MBE, NRMSE_pct, RPD, N)])

fwrite(results_all,  file.path(OUT_DIR,"S13_metrics_all.csv"))
fwrite(test_results, file.path(OUT_DIR,"S13_metrics_test.csv"))

all_preds_dt <- rbindlist(pred_list, use.names=TRUE)
all_preds_dt[, scene_id := rep(sample_agg$scene_id[idx_test], length(pred_list))]
all_preds_dt[, cal_year := rep(sample_agg$cal_year[idx_test],  length(pred_list))]
all_preds_dt[, cell_id  := rep(sample_agg$cell_id[idx_test],   length(pred_list))]
fwrite(all_preds_dt, file.path(OUT_DIR,"S13_predictions_test.csv"))
message("  CSVs written.")


# =============================================================================
# 18. FIGURES (PNG ONLY)
# =============================================================================
cat("\n[18] Generating figures ...\n")

.theme_ml <- function(base=11)
  theme_bw(base_size=base) +
  theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"))

.save_png <- function(p, fname, w=12, h=8, dpi=300) {
  ggsave(file.path(FIG_DIR,fname), plot=p, width=w, height=h, dpi=dpi, device="png")
  message(sprintf("  Saved: %s", fname))
}

MODEL_COLS <- c("LSTM"="#C62828","Transformer"="#6A1B9A",
                "Random Forest"="#1565C0","Gradient Boosting"="#2E7D32")
MODEL_COLS[[knn_label]] <- "#E65100"

# ── F1: Performance bar chart ──────────────────────────────────────────────────
perf_long <- melt(test_results[, .(ModelName,R2,RMSE,MAE,NRMSE_pct)],
                  id.vars="ModelName", variable.name="Metric", value.name="Value")
perf_long[, Metric := factor(Metric, levels=c("R2","RMSE","MAE","NRMSE_pct"),
                             labels=c("R²","RMSE (kg C/ha)","MAE (kg C/ha)","NRMSE (%)"))]

p_perf <- ggplot(perf_long, aes(x=ModelName, y=Value, fill=ModelName)) +
  geom_col(alpha=0.85, width=0.65, colour="white", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.3f",Value)), vjust=-0.4, size=3.0, fontface="bold") +
  scale_fill_brewer(palette="Set2", guide="none") +
  facet_wrap(~Metric, scales="free_y", ncol=4) +
  labs(title=sprintf("Model Performance — Monsoon CH4 (Test Set, %d spatial holdout cells)",
                     length(test_cells)),
       x=NULL, y="Value",
       caption="R²: higher better | RMSE/MAE/NRMSE: lower better | Day 166–309") +
  .theme_ml() + theme(axis.text.x=element_text(angle=30,hjust=1))
.save_png(p_perf, "S13_F1_performance_comparison.png", w=15, h=6)

# ── F2: Predicted vs Actual ────────────────────────────────────────────────────
all_preds_plot <- copy(all_preds_dt)
all_preds_plot[, model_f := factor(model,
                                   levels=c("LSTM","Transformer","Random Forest","Gradient Boosting",knn_label))]
lim_max <- max(c(all_preds_plot$actual, all_preds_plot$predicted), na.rm=TRUE)*1.05

r2_ann <- data.table(
  model_f=factor(c("LSTM","Transformer","Random Forest","Gradient Boosting",knn_label),
                 levels=c("LSTM","Transformer","Random Forest","Gradient Boosting",knn_label)),
  R2_v   =c(results_list[["LSTM_test"]]$R2, results_list[["TF_test"]]$R2,
            results_list[["RF_test"]]$R2,   results_list[["GBM_test"]]$R2,
            results_list[["KNN_test"]]$R2),
  RMSE_v =c(results_list[["LSTM_test"]]$RMSE, results_list[["TF_test"]]$RMSE,
            results_list[["RF_test"]]$RMSE,   results_list[["GBM_test"]]$RMSE,
            results_list[["KNN_test"]]$RMSE)
)
r2_ann[, label := sprintf("R²=%.3f\nRMSE=%.1f", R2_v, RMSE_v)]

p_scatter <- ggplot(all_preds_plot, aes(x=actual, y=predicted, colour=model_f)) +
  geom_abline(slope=1, intercept=0, colour="grey40", linewidth=0.7, linetype="dashed") +
  geom_point(alpha=0.20, size=0.7) +
  geom_smooth(method="lm", se=FALSE, linewidth=1.0, formula=y~x) +
  geom_text(data=r2_ann, aes(x=lim_max*0.04,y=lim_max*0.96,label=label),
            inherit.aes=FALSE, hjust=0, vjust=1, size=3.0, fontface="bold", colour="black") +
  scale_colour_manual(values=MODEL_COLS, guide="none") +
  coord_fixed(xlim=c(0,lim_max), ylim=c(0,lim_max)) +
  facet_wrap(~model_f, ncol=3L) +
  labs(title="Predicted vs Actual — Monsoon CH4",
       subtitle="Spatial holdout | Dashed=1:1 | Solid=OLS fit",
       x="Actual CH4 (kg C/ha/season)", y="Predicted CH4 (kg C/ha/season)",
       caption="DNDC 9.5 | Myanmar 2005–2024") +
  .theme_ml()
.save_png(p_scatter, "S13_F2_predicted_vs_actual.png", w=14, h=9)

# ── F3: Residual distributions ────────────────────────────────────────────────
p_resid <- ggplot(all_preds_plot, aes(x=residual, fill=model_f, colour=model_f)) +
  geom_histogram(aes(y=after_stat(density)), bins=50, alpha=0.40, linewidth=0.15) +
  geom_density(alpha=0, linewidth=0.9) +
  geom_vline(xintercept=0, colour="grey25", linewidth=0.7, linetype="dashed") +
  scale_fill_manual(values=MODEL_COLS, guide="none") +
  scale_colour_manual(values=MODEL_COLS, guide="none") +
  facet_wrap(~model_f, ncol=3L, scales="free_y") +
  labs(title="Residual Distributions — All Models",
       subtitle="Residual = Predicted – Actual | Centred at 0 = unbiased",
       x="Residual (kg C/ha/season)", y="Density", caption="Test set | Monsoon only") +
  .theme_ml()
.save_png(p_resid, "S13_F3_residual_distributions.png", w=14, h=8)

# ── F4a: LSTM training curves ─────────────────────────────────────────────────
p_lstm_loss <- ggplot(
  melt(data.table(epoch=seq_along(lstm_fit$train_hist),
                  Train=lstm_fit$train_hist, Val=lstm_fit$val_hist),
       id.vars="epoch", variable.name="Split", value.name="Loss"),
  aes(x=epoch, y=Loss, colour=Split)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=c(Train="#C62828",Val="#1565C0"), name="") +
  labs(title="LSTM Training Curves",
       subtitle=sprintf("hidden=%d | layers=%d | lr=%.0e | batch=%d | pat=%d",
                        LSTM_HIDDEN,LSTM_LAYERS,LSTM_LR,LSTM_BATCH,LSTM_PATIENCE),
       x="Epoch", y="MSE Loss (scaled log-target)") +
  .theme_ml()
.save_png(p_lstm_loss, "S13_F4a_LSTM_training_curves.png", w=10, h=5)

# ── F4b: Transformer training curves ──────────────────────────────────────────
p_tf_loss <- ggplot(
  melt(data.table(epoch=seq_along(tf_fit$train_hist),
                  Train=tf_fit$train_hist, Val=tf_fit$val_hist),
       id.vars="epoch", variable.name="Split", value.name="Loss"),
  aes(x=epoch, y=Loss, colour=Split)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=c(Train="#6A1B9A",Val="#AB47BC"), name="") +
  labs(title="Transformer Training Curves",
       subtitle=sprintf("d_model=%d | heads=%d | layers=%d | lr=%.0e | batch=%d",
                        TF_D_MODEL,TF_NHEAD,TF_LAYERS,TF_LR,TF_BATCH),
       x="Epoch", y="MSE Loss (scaled log-target)") +
  .theme_ml()
.save_png(p_tf_loss, "S13_F4b_Transformer_training_curves.png", w=10, h=5)

# ── F5: Feature importance ────────────────────────────────────────────────────
n_imp     <- min(15L, length(rf_imp))
rf_imp_dt <- data.table(Feature=names(rf_imp)[1:n_imp],
                        Importance=as.numeric(rf_imp)[1:n_imp],
                        Model="Random Forest")
xgb_imp_dt <- xgb_imp_raw[seq_len(min(15L,nrow(xgb_imp_raw))),
                          .(Feature, Importance=Gain, Model="Gradient Boosting")]
imp_all <- rbind(rf_imp_dt, xgb_imp_dt)
imp_all[, Feature := reorder(Feature, Importance, FUN=mean)]

p_imp <- ggplot(imp_all, aes(x=Feature, y=Importance, fill=Model)) +
  geom_col(position="dodge", alpha=0.85, width=0.7) +
  scale_fill_manual(values=c("Random Forest"="#1565C0","Gradient Boosting"="#2E7D32")) +
  coord_flip() +
  labs(title="Feature Importance — RF & XGBoost (Top 15)",
       subtitle="RF: impurity | XGBoost: gain | Flat features: {sum/mean/max/sd/q75/last} × seq + static",
       x=NULL, y="Importance") +
  .theme_ml(base=10)
.save_png(p_imp, "S13_F5_feature_importance.png", w=13, h=8)

# ── F6: Residuals by scenario ─────────────────────────────────────────────────
p_scene <- ggplot(all_preds_plot,
                  aes(x=scene_id, y=residual, fill=model_f, colour=model_f)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=0.5, linetype="dashed") +
  geom_boxplot(alpha=0.40, outlier.size=0.5, linewidth=0.45) +
  scale_fill_manual(values=MODEL_COLS, name="Model") +
  scale_colour_manual(values=MODEL_COLS, name="Model") +
  labs(title="Residuals by Scenario (Test Set)",
       subtitle="Zero = perfect. Spread shows scenario-specific bias.",
       x="Scenario ID", y="Residual (kg C/ha/season)", caption="Monsoon spatial holdout") +
  .theme_ml() +
  theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="bottom")
.save_png(p_scene, "S13_F6_residuals_by_scenario.png", w=16, h=7)

# ── F7: Calibration ───────────────────────────────────────────────────────────
calib_dt <- all_preds_plot[, {
  brks <- quantile(predicted, probs=seq(0,1,0.1), na.rm=TRUE)
  bin  <- as.integer(cut(predicted, breaks=brks, include.lowest=TRUE))
  tmp  <- data.table(predicted=predicted, actual=actual, bin=bin)
  tmp[!is.na(bin), .(mean_pred=mean(predicted), mean_actual=mean(actual), n=.N), by=bin]
}, by=model_f]

p_calib <- ggplot(calib_dt, aes(x=mean_actual, y=mean_pred, colour=model_f)) +
  geom_abline(slope=1, intercept=0, colour="grey40", linewidth=0.7, linetype="dashed") +
  geom_line(alpha=0.60, linewidth=0.7) +
  geom_point(aes(size=n), alpha=0.75) +
  scale_colour_manual(values=MODEL_COLS, name="Model") +
  scale_size_continuous(name="N", range=c(2,6)) +
  labs(title="Calibration — Mean Predicted vs Actual by Decile",
       subtitle="Dashed = perfect calibration",
       x="Mean Actual CH4 (kg C/ha)", y="Mean Predicted CH4 (kg C/ha)") +
  .theme_ml()
.save_png(p_calib, "S13_F7_calibration_decile.png", w=10, h=7)


# =============================================================================
# 19. EXCEL WORKBOOK (5 sheets)
# =============================================================================
cat("\n[19] Writing Excel workbook ...\n")

hparam_dt <- data.table(
  Model = c(rep("LSTM",7), rep("Transformer",7), rep("Random Forest",3),
            rep("Gradient Boosting",7), "KNN", rep("All",6)),
  Parameter = c(
    "hidden","n_layers","dropout","epochs_max","epochs_actual","lr","batch",
    "d_model","n_heads","n_layers","ff_dim","dropout","epochs_max","epochs_actual",
    "n_trees","mtry","min_node_size",
    "eta","max_depth","subsample","colsample_bytree","n_rounds_cv","n_rounds_final","min_child_weight",
    "k_neighbours",
    "seq_len","n_seq_feat","n_static_feat","n_flat_feat","test_cell_frac","train_years"
  ),
  Value = c(
    LSTM_HIDDEN,LSTM_LAYERS,LSTM_DROPOUT,LSTM_EPOCHS,lstm_fit$n_epochs,LSTM_LR,LSTM_BATCH,
    TF_D_MODEL,TF_NHEAD,TF_LAYERS,TF_FF_DIM,TF_DROPOUT,TF_EPOCHS,tf_fit$n_epochs,
    500,max(1L,floor(sqrt(n_feat_flat))),5,
    0.05,6,0.80,0.80,600,best_rounds,5,
    best_k,
    SEQ_LEN,n_seq_feat,n_static_feat,n_feat_flat,
    TEST_CELL_FRAC,paste0(min(TRAIN_YEARS),"-",max(TRAIN_YEARS))
  )
)

preds_export <- copy(all_preds_dt)
preds_export[, `:=`(actual=round(actual,4), predicted=round(predicted,4),
                    residual=round(residual,4))]

resid_summary <- all_preds_dt[, .(
  mean_res=round(mean(residual,na.rm=TRUE),3),
  sd_res  =round(sd(residual,  na.rm=TRUE),3),
  p5      =round(quantile(residual,0.05,na.rm=TRUE),3),
  p25     =round(quantile(residual,0.25,na.rm=TRUE),3),
  median  =round(median(residual,  na.rm=TRUE),3),
  p75     =round(quantile(residual,0.75,na.rm=TRUE),3),
  p95     =round(quantile(residual,0.95,na.rm=TRUE),3),
  pct_within_10pct=round(mean(abs(residual/(actual+1e-8))<0.10)*100,1)
), by=model]

wb <- createWorkbook()
addWorksheet(wb,"1_Overall_Metrics");  writeData(wb,1,results_all)
addWorksheet(wb,"2_Test_Metrics");     writeData(wb,2,test_results)
addWorksheet(wb,"3_Predictions");      writeData(wb,3,preds_export)
addWorksheet(wb,"4_Residual_Summary"); writeData(wb,4,resid_summary)
addWorksheet(wb,"5_Hyperparameters");  writeData(wb,5,hparam_dt)
conditionalFormatting(wb,2, cols=which(names(test_results)=="R2"),
                      rows=2:(nrow(test_results)+1),
                      style=c("#FCE4D6","#F4B8A4","#E26B5A"),
                      rule=c(0.5,0.8,0.95), type="colourScale")

xlsx_path <- file.path(OUT_DIR,"S13_model_results.xlsx")
saveWorkbook(wb, xlsx_path, overwrite=TRUE)
message(sprintf("  Saved: %s", basename(xlsx_path)))


# =============================================================================
# 20. SAVE SCALING OBJECTS & ARCH SPECS
# =============================================================================
cat("\n[20] Saving metadata ...\n")

saveRDS(list(
  y_log_mean=y_log_mean, y_log_sd=y_log_sd, back_transform=back_transform,
  seq_features=avail_seq, seq_means=seq_means, seq_sds=seq_sds,
  static_features=avail_static, static_means=static_means, static_sds=static_sds,
  flat_means=flat_means, flat_sds=flat_sds, flat_feature_names=colnames(X_flat_train),
  test_cells=test_cells, train_cells=train_cells,
  idx_train=idx_train, idx_val=idx_val, idx_test=idx_test,
  TRAIN_YEARS=TRAIN_YEARS, VAL_YEARS=VAL_YEARS,
  MONSOON_DOY=MONSOON_DOY, SEQ_LEN=SEQ_LEN
), file.path(MODEL_DIR,"S13_scaling_meta.rds"))

saveRDS(list(n_seq=n_seq_feat, n_static=n_static_feat,
             hidden=LSTM_HIDDEN, n_layers=LSTM_LAYERS, dropout=LSTM_DROPOUT),
        file.path(MODEL_DIR,"S13_lstm_arch.rds"))

saveRDS(list(n_seq=n_seq_feat, n_static=n_static_feat,
             d_model=TF_D_MODEL, nhead=TF_NHEAD, n_layers=TF_LAYERS,
             ff_dim=TF_FF_DIM, dropout=TF_DROPOUT, pe_matrix=pe_mat, seq_len=SEQ_LEN),
        file.path(MODEL_DIR,"S13_transformer_arch.rds"))

message("  Saved: S13_scaling_meta.rds | S13_lstm_arch.rds | S13_transformer_arch.rds")


# =============================================================================
# 21. CONSOLE SUMMARY
# =============================================================================
cat("\n", strrep("=", 72), "\n")
cat("  SCRIPT 13 v1.1 — COMPLETE\n")
cat(strrep("=", 72), "\n\n")

cat(sprintf("  Season   : monsoon only (Day %d–%d = %d days)\n",
            MONSOON_DOY[1], MONSOON_DOY[2], SEQ_LEN))
cat(sprintf("  Samples  : %s  (%d cells × %d scenes × %d years)\n",
            format(N_samples, big.mark=","),
            uniqueN(sample_agg$cell_id), uniqueN(sample_agg$scene_id),
            uniqueN(sample_agg$cal_year)))
cat(sprintf("  Split    : train=%d | val=%d | test=%d (%d cells, %.0f%% holdout)\n",
            length(idx_train), length(idx_val), length(idx_test),
            length(test_cells), TEST_CELL_FRAC*100))
cat(sprintf("  Features : %d seq × %d time-steps | %d static | %d flat\n",
            n_seq_feat, SEQ_LEN, n_static_feat, n_feat_flat))

cat("\n  TEST SET RESULTS (sorted by R²):\n")
for (i in seq_len(nrow(test_results))) {
  r <- test_results[i]
  cat(sprintf("  %-24s  R²=%6.4f  RMSE=%7.2f  MAE=%6.2f  NRMSE=%5.1f%%\n",
              r$Model, r$R2, r$RMSE, r$MAE, r$NRMSE_pct))
}
cat(sprintf("\n  LSTM       : %.1f s | %d epochs\n", elapsed_lstm, lstm_fit$n_epochs))
cat(sprintf("  Transformer: %.1f s | %d epochs\n",  elapsed_tf,   tf_fit$n_epochs))

cat("\n  MODELS:\n")
for (f in c("lstm_monsoon_state.pt","transformer_monsoon_state.pt",
            "rf_monsoon.rds","xgb_monsoon.bin","knn_monsoon.rds",
            "S13_scaling_meta.rds","S13_lstm_arch.rds","S13_transformer_arch.rds"))
  cat(sprintf("    %s/%s\n", MODEL_DIR, f))

cat("\n  FIGURES (PNG):\n")
for (f in c("S13_F1_performance_comparison.png","S13_F2_predicted_vs_actual.png",
            "S13_F3_residual_distributions.png","S13_F4a_LSTM_training_curves.png",
            "S13_F4b_Transformer_training_curves.png","S13_F5_feature_importance.png",
            "S13_F6_residuals_by_scenario.png","S13_F7_calibration_decile.png"))
  cat(sprintf("    %s/%s\n", FIG_DIR, f))

cat("\n  TABLES:\n")
for (f in c("S13_metrics_all.csv","S13_metrics_test.csv",
            "S13_predictions_test.csv","S13_model_results.xlsx"))
  cat(sprintf("    %s/%s\n", OUT_DIR, f))

cat("\n  Run: Rscript 13_monsoon_sequential_kgml.R\n")
cat(strrep("=", 72), "\n\n")

# ---- 2_dndc_make_climate_files_fix.R ----
##### 2_make_climate_files.R - DNDC Climate Files #####
# Generates DNDC-format climate files (one per cell per year)
# Format: Station header + daily rows (Julian_day Tmax Tmin Precip)

suppressPackageStartupMessages(library(terra))

# ---- Paths ----
path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
setwd(path)
power_dir <- file.path("data", "raw", "weather", "power")
out_dir <- file.path("dndc_inputs", "climate_files")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load cells table ----
cells <- readRDS("data/cells.rds")
years <- 2005:2024

message(sprintf("\n=== DNDC Climate File Generator ==="))
message(sprintf("Cells: %d | Years: %d | Total files: %d", nrow(cells), length(years), nrow(cells) * length(years)))

# ---- Main loop: one year at a time (memory efficient) ----
t0 <- Sys.time()
for (yr in years){
  ty <- Sys.time()
  n_days <- ifelse((yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0), 366, 365)
  
  # Load 3 climate variables for this year
  tmax <- rast(file.path(power_dir, "tmax-2005_2024-92x101x8.5x29.nc"))
  tmin <- rast(file.path(power_dir, "tmin-2005_2024-92x101x8.5x29.nc"))
  prec <- rast(file.path(power_dir, "prec-2005_2024-92x101x8.5x29.nc"))
  
  # Subset to this year's layers (handle leap years)
  year_start <- which(format(time(tmax), "%Y") == as.character(yr))[1]
  year_layers <- year_start:(year_start + n_days - 1)
  
  tmax_yr <- tmax[[year_layers]]
  tmin_yr <- tmin[[year_layers]]
  prec_yr <- prec[[year_layers]]
  
  # Extract for all cells at once
  tmax_vals <- as.matrix(tmax_yr[cells$cell])
  tmin_vals <- as.matrix(tmin_yr[cells$cell])
  prec_vals <- as.matrix(prec_yr[cells$cell])
  
  # Write file for each cell
  for (i in 1:nrow(cells)){
    cell_id <- cells$cell[i]
    fname <- file.path(out_dir, sprintf("Cell_%d_%d", cell_id, yr))
    
    # DNDC format: header + daily data (Julian, Tmax, Tmin, Precip)
    writeLines(sprintf("Cell_%d", cell_id), fname)
    write.table(
      data.frame(
        day = 1:n_days,
        tmax = round(tmax_vals[i, ], 2),
        tmin = round(tmin_vals[i, ], 2),
        prec = round(prec_vals[i, ], 2)
      ),
      fname, append = TRUE, quote = FALSE, sep = "  ", row.names = FALSE, col.names = FALSE
    )
  }
  
  message(sprintf("Year %d: %d files (%.1f sec)", yr, nrow(cells), difftime(Sys.time(), ty, units = "secs")))
  rm(tmax, tmin, prec, tmax_yr, tmin_yr, prec_yr, tmax_vals, tmin_vals, prec_vals); gc(verbose = FALSE)
}

# ---- Summary ----
dt <- difftime(Sys.time(), t0, units = "mins")
total_files <- nrow(cells) * length(years)
message(sprintf("\n=== COMPLETE ==="))
message(sprintf("Files created: %d", total_files))
message(sprintf("Total time: %.1f minutes (%.0f files/min)", dt, total_files / as.numeric(dt)))
message(sprintf("Output: %s", out_dir))

# ---- Validation sample ----
sample_file <- list.files(out_dir, pattern = "_2005$", full.names = TRUE)[1]
if (!is.na(sample_file)){
  message("\nSample file check:")
  message(sprintf("  File: %s", basename(sample_file)))
  sample_lines <- readLines(sample_file, n = 5)
  message("  First 5 lines:")
  for (ln in sample_lines) message(sprintf("    %s", ln))
}

message("\n✓ Climate files ready for DNDC simulation")
##### END OF: 2_make_climate_files.R #####

# ---- 3_dndc_make_input_files_MULTI_fixed.R ----
##### 3_dndc_make_input_files_MULTI_fixed.R #####
suppressPackageStartupMessages(library(dplyr))

path <- "C:\\DNDC\\simulate_dndc"
setwd(path)

template_dir <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\input_files\\template"
output_dir <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\input_files"
climate_base <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\climate_files"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- readRDS("C:\\DNDC\\simulate_dndc\\data\\cells.rds")
years <- 2005:2024
n_years <- length(years)
scenarios <- sprintf("MMR_L1_%03d", 0:13)

message(sprintf("Cells: %d | Scenarios: %d | Total files: %d", nrow(cells), length(scenarios), nrow(cells) * length(scenarios)))

get_soil <- function(row, col) {
  val <- row[[col]]
  if (is.null(val) || length(val) == 0 || is.na(val)) NA_real_ else as.numeric(val)
}

calc_fc_wp <- function(clay_frac, sand_frac, soc_frac) {
  clay_pct <- clay_frac * 100
  sand_pct <- sand_frac * 100
  soc_pct <- ifelse(is.na(soc_frac), 2, soc_frac * 100)
  fc <- 0.2576 - 0.002 * sand_pct + 0.0036 * clay_pct + 0.0299 * soc_pct
  wp <- 0.026 + 0.005 * clay_pct + 0.0158 * soc_pct
  list(fc = max(0.15, min(0.55, fc)), wp = max(0.05, min(0.25, wp)))
}

set_line <- function(lines, pattern, value_str) {
  idx <- grep(pattern, lines, perl = TRUE)
  if (length(idx) == 0) return(lines)
  key <- sub("^(_{2,6}\\S+).*", "\\1", lines[idx[1]])
  lines[idx[1]] <- formatC(key, width = -51, flag = "-") |> paste0(value_str)
  lines
}

generate_dnd <- function(cell_row, years, template) {
  id <- cell_row$cell
  lines <- template

  lines <- set_line(lines, "^__Site_name\\s", sprintf("Myanmar_Rice_Cell_%d", id))
  lines <- set_line(lines, "^__Simulated_years\\s", sprintf("%d", n_years))
  lines <- set_line(lines, "^__Latitude\\s", sprintf("%.4f", cell_row$y))

  cf_idx <- grep("^__Climate_files\\s", lines)
  mode_idx <- grep("^__Climate_file_mode\\s", lines)
  old_rows <- which(grepl("^[0-9]+\\s+", lines) & seq_along(lines) > cf_idx & seq_along(lines) < mode_idx)
  if (length(old_rows)) lines <- lines[-old_rows]

  new_cf <- sapply(seq_along(years), function(i)
    sprintf("%-2d                                       %s\\Cell_%d_%d", i, climate_base, id, years[i]))

  lines <- set_line(lines, "^__Climate_files\\s", sprintf("%d", n_years))
  lines <- append(lines, new_cf, after = grep("^__Climate_files\\s", lines))
  lines <- set_line(lines, "^__Total_years\\s", sprintf("%d", n_years))

  bd <- get_soil(cell_row, "bdod_0-5cm")
  ph <- get_soil(cell_row, "pH_0-5cm")
  clay <- get_soil(cell_row, "clay_0-5cm")
  sand <- get_soil(cell_row, "sand_0-5cm")
  soc <- get_soil(cell_row, "soc_0-5cm")
  poros <- get_soil(cell_row, "poros_0-5cm")
  tex <- get_soil(cell_row, "texture_0-5cm")

  if (!is.na(bd)) lines <- set_line(lines, "^__Bulk_density\\s", sprintf("%.4f", bd))
  if (!is.na(ph)) lines <- set_line(lines, "^__pH\\s", sprintf("%.4f", ph))
  if (!is.na(clay)) lines <- set_line(lines, "^__Clay_fraction\\s", sprintf("%.4f", clay))
  if (!is.na(poros)) lines <- set_line(lines, "^__Porosity\\s", sprintf("%.4f", poros))
  if (!is.na(soc)) lines <- set_line(lines, "^__Top_layer_SOC\\s", sprintf("%.4f", soc))
  if (!is.na(tex)) lines <- set_line(lines, "^__Soil_texture_ID\\s", sprintf("%d", as.integer(tex)))

  if (!is.na(clay) && !is.na(sand)) {
    hw <- calc_fc_wp(clay, sand, soc)
    lines <- set_line(lines, "^__Field_capacity\\s", sprintf("%.4f", hw$fc))
    lines <- set_line(lines, "^__Wilting_point\\s", sprintf("%.4f", hw$wp))
  }

  lines
}

message("Generating input files...")
dnd_files <- character(nrow(cells) * length(scenarios))
file_idx <- 1

for (s in seq_along(scenarios)) {
  scen <- scenarios[s]
  scen_code <- sub(".*_", "", scen)
  template_file <- file.path(template_dir, paste0(scen, ".dnd"))

  if (!file.exists(template_file)) {
    message(sprintf("Warning: Template not found: %s", basename(template_file)))
    next
  }

  template <- readLines(template_file)

  for (i in seq_len(nrow(cells))) {
    lines <- generate_dnd(cells[i, ], years, template)
    out <- file.path(output_dir, sprintf("Cell_%d_S%s.dnd", cells$cell[i], scen_code))
    writeLines(lines, out)
    dnd_files[file_idx] <- out
    file_idx <- file_idx + 1
  }

  message(sprintf("  Scenario S%s (%s): %d files", scen_code, scen, nrow(cells)))
}

dnd_files <- dnd_files[1:(file_idx - 1)]

batch_file <- file.path(output_dir, "myanmar_rice_batch_all.txt")
batch_lines <- c(as.character(length(dnd_files)), normalizePath(dnd_files, winslash = "\\", mustWork = FALSE))
writeLines(batch_lines, batch_file)

s <- readLines(dnd_files[1])
message(sprintf("\n=== Complete: %d files, batch -> %s", length(dnd_files), basename(batch_file)))
message(sprintf("Sample file: %s", basename(dnd_files[1])))
message(grep("^__Site_name", s, value = TRUE)[1])
message(grep("^__Climate_files", s, value = TRUE)[1])

##### END OF: 3_dndc_make_input_files_MULTI_fixed.R #####

# ---- 4_dndc_make_batch.R ----
##### 4_dndc_make_batch_MULTI.R #####
# Creates DNDC batch file listing all multi-scenario cell .dnd input files

# ---- Paths ----
input_dir  <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\input_files"
batch_file <- file.path(input_dir, "myanmar_rice_batch_all.txt")

# ---- Find and sort all scenario .dnd files ----
# Expected pattern from script 3:
# Cell_<cellid>_S<scenario>.dnd
dnd_files <- list.files(
  input_dir,
  pattern = "^Cell_[0-9]+_S[0-9]{3}\\.dnd$",
  full.names = TRUE
)

if (length(dnd_files) == 0) stop("No multi-scenario .dnd files found in: ", input_dir)

# ---- Extract cell ID and scenario code for sorting ----
bn <- basename(dnd_files)
cell_nums <- as.integer(sub("^Cell_([0-9]+)_S[0-9]{3}\\.dnd$", "\\1", bn))
scen_nums <- as.integer(sub("^Cell_[0-9]+_S([0-9]{3})\\.dnd$", "\\1", bn))

ord <- order(scen_nums, cell_nums)
dnd_files <- dnd_files[ord]
cell_nums <- cell_nums[ord]
scen_nums <- scen_nums[ord]

# ---- Write batch file ----
batch_lines <- c(
  as.character(length(dnd_files)),
  normalizePath(dnd_files, winslash = "\\", mustWork = FALSE)
)

writeLines(batch_lines, batch_file)

# ---- Summary ----
cat(sprintf("=== DNDC Batch File Created ===\n"))
cat(sprintf("Output    : %s\n", batch_file))
cat(sprintf("Files     : %d\n", length(dnd_files)))
cat(sprintf("Scenarios : %d  (S%03d to S%03d)\n",
            length(unique(scen_nums)), min(scen_nums), max(scen_nums)))
cat(sprintf("Cells     : %d  (Cell_%d to Cell_%d)\n",
            length(unique(cell_nums)), min(cell_nums), max(cell_nums)))

cat(sprintf("\nPreview (first 6 lines):\n"))
cat(paste(readLines(batch_file, n = 6), collapse = "\n"), "\n")

cat(sprintf("\nTo run: Open DNDC -> Tools -> Run Batch -> select %s\n",
            basename(batch_file)))

##### END OF: 4_dndc_make_batch_MULTI.R #####

# ---- 5_dndc_extract_results_fix_v2.R ----
###############################################################################
# Script 5 — Extract and assemble DNDC batch outputs
#
# Inputs:
#   - C:/DNDC/simulate_dndc/dndc_inputs/input_files/myanmar_rice_batch_all.txt
#   - C:/DNDC/Result/Record/Batch/Case*
#
# Outputs:
#   - case_scene_index.csv
#   - annual_allcases.rds / annual_allcases.csv
#   - daily_SoilC_allcases.rds
#   - daily_SoilClimate_allcases.rds
#   - daily_SoilWater_allcases.rds
#   - daily_FieldCrop_allcases.rds
#   - daily_panel_allcases.rds
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(future)
  library(future.apply)
})

# =============================================================================
# 0. CONFIG
# =============================================================================

ROOT       <- "C:/DNDC/simulate_dndc"
INPUT_DIR  <- file.path(ROOT, "dndc_inputs", "input_files")
OUTPUT_DIR <- file.path(ROOT, "dndc_outputs")

BATCH_TXT  <- file.path(INPUT_DIR, "myanmar_rice_batch_all.txt")
BATCH_DIR  <- "C:/DNDC/Result/Record/Batch"

N_YEARS    <- 20L
START_YEAR <- 2005L
N_WORKERS  <- 4L

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(BATCH_TXT)) stop("Missing batch file: ", BATCH_TXT)
if (!dir.exists(BATCH_DIR)) stop("Missing DNDC batch output folder: ", BATCH_DIR)

# =============================================================================
# 1. CASE LOOKUP FROM BATCH FILE
# =============================================================================

raw <- trimws(gsub("\r", "", readLines(BATCH_TXT, warn = FALSE)))
n_declared <- suppressWarnings(as.integer(raw[1]))
dnd_lines <- raw[-1]
dnd_lines <- dnd_lines[nchar(dnd_lines) > 5]

LOOKUP <- data.table(case_id = seq_along(dnd_lines), dnd_path = dnd_lines)
LOOKUP[, cell_id  := as.integer(str_extract(dnd_path, "(?<=Cell_)\\d+(?=_S)"))]
LOOKUP[, scene_id := str_extract(dnd_path, "S\\d{3}(?=\\.dnd)")]
LOOKUP <- LOOKUP[!is.na(cell_id) & !is.na(scene_id)]

fwrite(
  LOOKUP[, .(case_id, cell_id, scene_id)],
  file.path(OUTPUT_DIR, "case_scene_index.csv")
)

message(sprintf(
  "Batch lookup: %d declared | %d parsed | %d cells | %d scenes",
  n_declared, nrow(LOOKUP), uniqueN(LOOKUP$cell_id), uniqueN(LOOKUP$scene_id)
))

# =============================================================================
# 2. MATCH DNDC OUTPUT FOLDERS
# =============================================================================

dirs <- list.dirs(BATCH_DIR, full.names = TRUE, recursive = FALSE)
dirs <- dirs[grepl("^Case\\d+-Myanmar_Rice_Cell_\\d+$", basename(dirs))]

if (!length(dirs)) {
  stop("No DNDC Case folders found in: ", BATCH_DIR,
       "\nExpected folder pattern: Case{N}-Myanmar_Rice_Cell_{M}")
}

case_index <- data.table(
  folder  = dirs,
  case_id = as.integer(str_extract(basename(dirs), "(?<=Case)\\d+"))
)

case_index <- merge(
  case_index,
  LOOKUP[, .(case_id, cell_id, scene_id)],
  by = "case_id",
  all.x = TRUE
)

case_index <- case_index[!is.na(scene_id)]
setorder(case_index, case_id)

message(sprintf(
  "Matched folders: %d | cells: %d | scenes: %d",
  nrow(case_index), uniqueN(case_index$cell_id), uniqueN(case_index$scene_id)
))

# =============================================================================
# 3. COLUMN SPECS
# =============================================================================

ANNUAL_COLS <- c(
  "Year",
  "GrainC1","LeafStemC1","RootC1",
  "GrainC2","LeafStemC2","RootC2",
  "GrainC3","LeafStemC3","RootC3",
  "SOC_0_10cm","SOC_0_20cm","SOC_0_30cm",
  "Ini_SOC","End_SOC","dSOC",
  "LitterC_input","RootC_input","ManureC_input",
  "Soil_CO2","CH4",
  "Ini_SON","Ini_SIN","End_SON","End_SIN","dSN",
  "Atmo_N_input","Fertilizer_N_input","Manure_N_input",
  "Litter_N_input","N_fixation",
  "Crop_N_uptake","N_leach","N_runoff",
  "N2O_flux","NO_flux","N2_flux","NH3_flux","Exch_NH4",
  "PET","Transpiration","Evaporation","WaterLeach","Runoff",
  "Irrigation","Precipitation",
  "MeanT","WindSpeed",
  "ColdStress","WaterStress","N_Stress",
  "Cut_CropC"
)

KEEP <- list(
  SoilC = c(
    "Day", "SOC", "dSOC", "DOC", "Microbe", "Humads", "Humus",
    "NPP", "NEE", "Photosynthesis", "Soil-heterotrophic-respiration",
    "CH4-DOC", "CH4-prod.", "CH4-oxid.", "CH4-flux", "CH4-pool",
    "Litter-C", "Manure-C", "DOC-leach"
  ),
  SoilWater = c(
    "Day", "Ponding", "Precipitation", "Irrigation",
    "Evaporation", "Transpiration", "Leaching", "Runoff",
    "IniSoilWater", "EndSoilWater"
  ),
  FieldCrop = c(
    "Day", "LeafC", "StemC", "RootC", "GrainC",
    "TDD", "GrowthIndex", "Water_stress", "N_stress", "LAI",
    "TotalCropN", "DailyCropGrowth", "DayGrainGrowth"
  )
)

SOILCLIM_SRC <- c(
  "Day",
  "1","5","10","20",
  "1.1","5.1","10.1","20.1",
  "1.3","10.3","20.3",
  "(mm).2","(mm).3",
  "1.5"
)

SOILCLIM_DST <- c(
  "Day",
  "SoilT_1cm","SoilT_5cm","SoilT_10cm","SoilT_20cm",
  "SoilM_1cm","SoilM_5cm","SoilM_10cm","SoilM_20cm",
  "SoilEh_1cm","SoilEh_10cm","SoilEh_20cm",
  "SoilWater_mm","DeepWater_mm",
  "Soil_pH_1cm"
)

# =============================================================================
# 4. READER FUNCTIONS
# =============================================================================

drop_bad <- function(dt) {
  dt[suppressWarnings(!is.na(as.numeric(as.character(dt[[1]]))))]
}

tag_rows <- function(dt, case_id, cell_id, scene_id, yr) {
  dt[, Day := as.integer(Day)]
  dt[, `:=`(
    sim_year = yr,
    cal_year = yr + START_YEAR - 1L,
    case_id  = case_id,
    cell_id  = cell_id,
    scene_id = scene_id
  )]
  dt
}

read_annual <- function(path, case_id, cell_id, scene_id) {
  dt <- tryCatch(fread(path, header = FALSE, fill = TRUE, showProgress = FALSE),
                 error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  
  n_named <- min(ncol(dt), length(ANNUAL_COLS))
  new_names <- ANNUAL_COLS[seq_len(n_named)]
  if (ncol(dt) > n_named) {
    new_names <- c(new_names, paste0("V", seq.int(n_named + 1L, ncol(dt))))
  }
  setnames(dt, new_names)
  
  dt[, `:=`(
    sim_year = as.integer(Year),
    cal_year = as.integer(Year) + START_YEAR - 1L,
    case_id  = case_id,
    cell_id  = cell_id,
    scene_id = scene_id
  )]
  dt
}

read_skip1 <- function(path, keep_cols, case_id, cell_id, scene_id, yr) {
  dt <- tryCatch(fread(path, skip = 1L, header = TRUE, fill = TRUE,
                       showProgress = FALSE),
                 error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  
  dt <- drop_bad(dt)
  cols <- keep_cols[keep_cols %in% names(dt)]
  if (!"Day" %in% cols) return(NULL)
  
  tag_rows(dt[, ..cols], case_id, cell_id, scene_id, yr)
}

read_soilclimate <- function(path, case_id, cell_id, scene_id, yr) {
  lns <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(lns)) return(NULL)
  
  h <- which(grepl("^Day,", lns))[1]
  if (is.na(h)) return(NULL)
  
  dt <- tryCatch(
    fread(text = paste(lns[h:length(lns)], collapse = "\n"),
          header = TRUE, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || !nrow(dt)) return(NULL)
  
  dt <- drop_bad(dt)
  present <- SOILCLIM_SRC %in% names(dt)
  if (!any(present)) return(NULL)
  
  out <- dt[, SOILCLIM_SRC[present], with = FALSE]
  setnames(out, SOILCLIM_DST[present])
  tag_rows(out, case_id, cell_id, scene_id, yr)
}

read_fieldcrop <- function(path, case_id, cell_id, scene_id, yr) {
  hdr <- tryCatch(readLines(path, n = 3L, warn = FALSE), error = function(e) NULL)
  if (is.null(hdr) || length(hdr) < 3L) return(NULL)
  
  r1 <- trimws(strsplit(hdr[2L], ",")[[1]])
  r2 <- trimws(strsplit(hdr[3L], ",")[[1]])
  n  <- max(length(r1), length(r2))
  r1 <- c(r1, rep("", n - length(r1)))
  r2 <- c(r2, rep("", n - length(r2)))
  out_names <- ifelse(r2 != "", r2, r1)
  
  dt <- tryCatch(fread(path, skip = 4L, header = FALSE, fill = TRUE,
                       showProgress = FALSE),
                 error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  
  setnames(dt, seq_len(min(ncol(dt), length(out_names))),
           out_names[seq_len(min(ncol(dt), length(out_names)))])
  
  dt <- drop_bad(dt)
  cols <- KEEP$FieldCrop[KEEP$FieldCrop %in% names(dt)]
  if (!"Day" %in% cols) return(NULL)
  
  tag_rows(dt[, ..cols], case_id, cell_id, scene_id, yr)
}

read_years <- function(folder, prefix, reader, case_id, cell_id, scene_id) {
  out <- vector("list", N_YEARS)
  
  for (yr in seq_len(N_YEARS)) {
    p <- file.path(folder, sprintf("%s_%d.csv", prefix, yr))
    if (!file.exists(p)) next
    out[[yr]] <- tryCatch(reader(p, case_id, cell_id, scene_id, yr),
                          error = function(e) NULL)
  }
  
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(NULL)
  rbindlist(out, fill = TRUE, use.names = TRUE)
}

extract_case <- function(folder, case_id, cell_id, scene_id) {
  ann_path <- file.path(folder, "Multi_year_summary.csv")
  
  list(
    Annual = if (file.exists(ann_path)) {
      read_annual(ann_path, case_id, cell_id, scene_id)
    } else NULL,
    
    SoilC = read_years(
      folder, "Day_SoilC",
      function(p, ci, cid, sid, yr) read_skip1(p, KEEP$SoilC, ci, cid, sid, yr),
      case_id, cell_id, scene_id
    ),
    
    SoilClimate = read_years(
      folder, "Day_SoilClimate",
      read_soilclimate,
      case_id, cell_id, scene_id
    ),
    
    SoilWater = read_years(
      folder, "Day_SoilWater",
      function(p, ci, cid, sid, yr) read_skip1(p, KEEP$SoilWater, ci, cid, sid, yr),
      case_id, cell_id, scene_id
    ),
    
    FieldCrop = read_years(
      folder, "Day_FieldCrop",
      read_fieldcrop,
      case_id, cell_id, scene_id
    )
  )
}

# =============================================================================
# 5. RUN EXTRACTION
# =============================================================================

plan(multisession, workers = N_WORKERS)
t0 <- proc.time()

res <- future_mapply(
  FUN = extract_case,
  folder = case_index$folder,
  case_id = case_index$case_id,
  cell_id = case_index$cell_id,
  scene_id = case_index$scene_id,
  SIMPLIFY = FALSE,
  future.seed = TRUE
)

plan(sequential)
elapsed_min <- (proc.time() - t0)[["elapsed"]] / 60

# =============================================================================
# 6. SAVE COMPONENT FILES
# =============================================================================

TYPES <- c("Annual", "SoilC", "SoilClimate", "SoilWater", "FieldCrop")
SORT_KEYS <- c("scene_id", "case_id", "cell_id", "sim_year", "cal_year", "Day")

inventory <- data.table(
  file_type = TYPES,
  rows = 0L,
  cols = 0L,
  size_mb = 0,
  status = "missing"
)

for (ft in TYPES) {
  parts <- Filter(function(x) !is.null(x) && nrow(x), lapply(res, `[[`, ft))
  if (!length(parts)) next
  
  dt <- rbindlist(parts, fill = TRUE, use.names = TRUE)
  rm(parts); gc()
  
  for (j in intersect(c("case_id", "cell_id", "sim_year", "cal_year", "Year", "Day"), names(dt))) {
    set(dt, j = j, value = as.integer(dt[[j]]))
  }
  
  setkeyv(dt, intersect(SORT_KEYS, names(dt)))
  
  fn <- if (ft == "Annual") "annual_allcases.rds" else sprintf("daily_%s_allcases.rds", ft)
  fp <- file.path(OUTPUT_DIR, fn)
  
  saveRDS(dt, fp)
  if (ft == "Annual") fwrite(dt, file.path(OUTPUT_DIR, "annual_allcases.csv"))
  
  inventory[file_type == ft, `:=`(
    rows = nrow(dt),
    cols = ncol(dt),
    size_mb = round(file.info(fp)$size / 1024^2, 1),
    status = "saved"
  )]
  
  rm(dt); gc()
}

# =============================================================================
# 7. BUILD DAILY PANEL
# =============================================================================

build_daily_panel <- function(output_dir) {
  panel_path <- file.path(output_dir, "daily_panel_allcases.rds")
  ftypes <- c("SoilC", "SoilClimate", "SoilWater", "FieldCrop")
  files <- setNames(file.path(output_dir, sprintf("daily_%s_allcases.rds", ftypes)), ftypes)
  
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing daily component files:\n", paste(missing, collapse = "\n"))
  
  join_keys <- c("case_id", "cell_id", "scene_id", "sim_year", "cal_year", "Day")
  
  comp <- lapply(files, function(f) {
    dt <- readRDS(f)
    for (j in intersect(c("case_id", "cell_id", "sim_year", "cal_year", "Day"), names(dt))) {
      set(dt, j = j, value = as.integer(dt[[j]]))
    }
    setkeyv(dt, intersect(join_keys, names(dt)))
    dt
  })
  
  scenes <- sort(unique(comp[[which.max(sapply(comp, nrow))]]$scene_id))
  chunks <- vector("list", length(scenes))
  
  for (i in seq_along(scenes)) {
    sid <- scenes[i]
    parts <- lapply(comp, function(dt) dt[scene_id == sid])
    parts <- Filter(function(x) !is.null(x) && nrow(x), parts)
    
    chunks[[i]] <- Reduce(
      function(a, b) {
        by_cols <- intersect(join_keys, intersect(names(a), names(b)))
        merge(a, b, by = by_cols, all = TRUE, sort = FALSE)
      },
      parts[-1],
      init = parts[[1]]
    )
  }
  
  rm(comp); gc()
  
  panel <- rbindlist(chunks, fill = TRUE, use.names = TRUE)
  setkeyv(panel, intersect(c("scene_id", join_keys), names(panel)))
  
  saveRDS(panel, panel_path)
  out <- data.table(
    file_type = "daily_panel",
    rows = nrow(panel),
    cols = ncol(panel),
    size_mb = round(file.info(panel_path)$size / 1024^2, 1),
    status = "saved"
  )
  
  rm(panel, chunks); gc()
  out
}

panel_inventory <- build_daily_panel(OUTPUT_DIR)
inventory <- rbind(inventory, panel_inventory, fill = TRUE)

# =============================================================================
# 8. SUMMARY
# =============================================================================

print(inventory)
message(sprintf(
  "Script 5 complete: %d matched cases | %d cells | %d scenes | %.1f min",
  nrow(case_index), uniqueN(case_index$cell_id), uniqueN(case_index$scene_id), elapsed_min
))

# ---- 8_dndc_enriched_daily_v2.R ----
###############################################################################
# Script 8 — Enrich daily DNDC panel with cell and scenario metadata
#
# Inputs:
#   - daily_panel_allcases.rds
#   - cells.rds
#   - scene_registry.csv
#
# Outputs:
#   - enriched_daily_panel.rds
#   - enriched_daily_panel.fst, if fst is installed
#   - enrichment_audit.txt
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

USE_FST <- requireNamespace("fst", quietly = TRUE)
if (USE_FST) library(fst)

# =============================================================================
# 0. CONFIG
# =============================================================================

ROOT       <- "C:/DNDC/simulate_dndc"
OUTPUT_DIR <- file.path(ROOT, "dndc_outputs")

PANEL_RDS <- file.path(OUTPUT_DIR, "daily_panel_allcases.rds")
CELLS_RDS <- file.path(ROOT, "dndc_inputs", "input_files", "processing", "cells.rds")
SCENE_CSV <- file.path(ROOT, "dndc_inputs", "input_files", "processing", "scene_registry.csv")

OUT_RDS   <- file.path(OUTPUT_DIR, "enriched_daily_panel.rds")
OUT_FST   <- file.path(OUTPUT_DIR, "enriched_daily_panel.fst")
AUDIT_TXT <- file.path(OUTPUT_DIR, "enrichment_audit.txt")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (f in c(PANEL_RDS, CELLS_RDS, SCENE_CSV)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

# =============================================================================
# 1. HELPERS
# =============================================================================

rename_if_present <- function(dt, old, new) {
  hit <- old %in% names(dt)
  if (any(hit)) setnames(dt, old[hit], new[hit])
  invisible(dt)
}

first_present <- function(dt, candidates) {
  x <- candidates[candidates %in% names(dt)]
  if (length(x)) x[1] else NA_character_
}

need_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(label, " missing columns: ", paste(miss, collapse = ", "))
}

rc_sum <- function(dt, pattern) {
  cols <- grep(pattern, names(dt), value = TRUE)
  if (!length(cols)) return(rep(0, nrow(dt)))
  rowSums(dt[, ..cols], na.rm = TRUE)
}

as_int0 <- function(x) as.integer(fifelse(is.na(x), 0L, as.integer(x)))
as_num0 <- function(x) fifelse(is.na(x), 0, as.numeric(x))

# =============================================================================
# 2. LOAD INPUTS
# =============================================================================

panel <- readRDS(PANEL_RDS)
setDT(panel)

cells <- readRDS(CELLS_RDS)
setDT(cells)

scene_raw <- fread(SCENE_CSV)

# =============================================================================
# 3. PREPARE CELL LOOKUP
# =============================================================================

cell_id_col <- first_present(
  cells,
  c("cell_id", "cellid", "cell", "id", "pixel_id", "grid_id", "cell_no")
)
if (is.na(cell_id_col)) stop("Cannot identify cell_id column in cells.rds")

setnames(cells, cell_id_col, "cell_id")

rename_if_present(
  cells,
  old = c("x", "y",
          "clay_0-5cm", "bdod_0-5cm", "pH_0-5cm", "soc_0-5cm",
          "poros_0-5cm", "sand_0-5cm", "silt_0-5cm", "elevation"),
  new = c("lon", "lat",
          "clay_fraction", "bulk_density", "soil_pH", "SOC_static",
          "porosity", "sand_fraction", "silt_fraction", "elevation_m")
)

cell_keep <- intersect(
  c("cell_id", "lon", "lat",
    "clay_fraction", "bulk_density", "soil_P", "soil_pH", "SOC_static",
    "porosity", "sand_fraction", "silt_fraction", "elevation_m",
    "texture_0-5cm", "field_capacity", "wilting_point",
    "hydro_conductivity", "watertable_depth", "soil_slope"),
  names(cells)
)

cells <- unique(cells[, ..cell_keep], by = "cell_id")

# =============================================================================
# 4. PREPARE SCENARIO LOOKUP
# =============================================================================

sid_col <- first_present(scene_raw, c("scene_id", "Scenario_ID", "scenario_id", "Scene_ID"))
if (is.na(sid_col)) stop("Cannot identify scene_id column in scene_registry.csv")

setnames(scene_raw, sid_col, "scene_id")
scene_raw[, scene_id := as.character(scene_id)]

SCENE_LABELS <- c(
  S000 = "Baseline",   S001 = "Injected fert", S002 = "Med-high fert",
  S003 = "Extreme fert", S004 = "Split fert", S005 = "Short flood",
  S006 = "AWD",        S008 = "No-till",       S009 = "Deep plow",
  S010 = "Continuous irrigation", S011 = "Rainfed",
  S012 = "Low residue", S013 = "High residue"
)

SCENE_CATS <- c(
  S000 = "baseline", S001 = "fertilizer", S002 = "fertilizer",
  S003 = "fertilizer", S004 = "fertilizer", S005 = "water",
  S006 = "water", S008 = "tillage", S009 = "tillage",
  S010 = "irrigation", S011 = "irrigation",
  S012 = "residue", S013 = "residue"
)

urea_total <- rc_sum(scene_raw, "^F[0-9]+_Urea$")
amm_total  <- rc_sum(scene_raw, "^F[0-9]+_Ammonium(_02)?$")
phos_total <- rc_sum(scene_raw, "^F[0-9]+_Phosphate$")
total_N    <- urea_total + amm_total

fert_method_code <- as_int0(scene_raw$F1_Fertilizing_method)
fert_method_label <- fcase(
  fert_method_code == 0L, "broadcast",
  fert_method_code == 1L, "banding",
  fert_method_code == 2L, "injection",
  default = "broadcast"
)

fl1_days <- pmax(
  0L,
  (as_int0(scene_raw$FL1_End_month) - as_int0(scene_raw$FL1_Start_month)) * 30L +
    (as_int0(scene_raw$FL1_End_day) - as_int0(scene_raw$FL1_Start_day))
)

fl2_start_month <- as_int0(scene_raw$FL1_End_month)

fl2_days <- pmax(
  0L,
  (as_int0(scene_raw$FL2_End_month) - fl2_start_month) * 30L +
    (as_int0(scene_raw$FL2_End_day) - as_int0(scene_raw$FL2_Start_day))
)

total_flood_days <- fl1_days + fl2_days

awd_flag <- as.integer(
  as_int0(scene_raw$FL1_Alter_wet_dry) == 1L |
    as_int0(scene_raw$FL2_Alter_wet_dry) == 1L
)

till_apps <- as_int0(scene_raw$CS1_Till_applications)
till_method_code <- as_int0(scene_raw$T1_Till_method)

till_method_label <- fcase(
  till_method_code == 0L, "none",
  till_method_code == 1L, "no_till",
  till_method_code == 2L, "ridge",
  till_method_code == 3L, "chisel",
  till_method_code == 4L, "disk",
  till_method_code == 5L, "moldboard",
  default = "none"
)

SCENES <- data.table(
  scene_id = scene_raw$scene_id,
  scene_label = fifelse(scene_raw$scene_id %in% names(SCENE_LABELS),
                        SCENE_LABELS[scene_raw$scene_id], scene_raw$scene_id),
  scene_category = fifelse(scene_raw$scene_id %in% names(SCENE_CATS),
                           SCENE_CATS[scene_raw$scene_id], "unknown"),
  Site_name = if ("Site_name" %in% names(scene_raw)) scene_raw$Site_name else NA_character_,
  
  fert_total_N_kgNha   = round(total_N, 2),
  fert_urea_kgNha      = round(urea_total, 2),
  fert_ammonium_kgNha  = round(amm_total, 2),
  fert_phosphate_kgPha = round(phos_total, 2),
  fert_n_applications  = as_int0(scene_raw$CS1_Fertilizer_applications),
  fert_method_code     = fert_method_code,
  fert_method_label    = fert_method_label,
  
  residue_summer_frac  = as.numeric(scene_raw$C1_Residue_left_in_field),
  residue_monsoon_frac = as.numeric(scene_raw$C2_Residue_left_in_field),
  
  FL1_flood_days = fl1_days,
  FL2_flood_days = fl2_days,
  total_flood_days = total_flood_days,
  awd_flag = awd_flag,
  irrigation_control_code = as_int0(scene_raw$CS1_Irrigation_control),
  
  till_applications = till_apps,
  till_method_code = till_method_code,
  till_method_label = till_method_label,
  
  max_yield_summer_kgha = as.numeric(scene_raw$C1_Maximum_yield),
  planting_month_summer = as_int0(scene_raw$C1_Planting_month),
  harvest_month_summer = as_int0(scene_raw$C1_Harvest_month),
  harvest_month_monsoon = as_int0(scene_raw$FL2_End_month) + 1L
)

SCENES <- unique(SCENES, by = "scene_id")

# =============================================================================
# 5. JOIN METADATA
# =============================================================================

need_cols(panel, c("cell_id", "scene_id"), "daily panel")

original_cols <- names(panel)

panel[, cell_id := as.integer(cell_id)]
panel[, scene_id := as.character(scene_id)]

panel <- merge(panel, cells, by = "cell_id", all.x = TRUE, sort = FALSE)
panel <- merge(panel, SCENES, by = "scene_id", all.x = TRUE, sort = FALSE)

setkeyv(panel, intersect(
  c("scene_id", "case_id", "cell_id", "sim_year", "cal_year", "Day"),
  names(panel)
))

# =============================================================================
# 6. SAVE OUTPUTS
# =============================================================================

saveRDS(panel, OUT_RDS)
sz_rds <- round(file.info(OUT_RDS)$size / 1024^2, 1)

sz_fst <- NA_real_
if (USE_FST) {
  fst::write_fst(panel, OUT_FST, compress = 50)
  sz_fst <- round(file.info(OUT_FST)$size / 1024^2, 1)
}

# =============================================================================
# 7. AUDIT
# =============================================================================

scene_cols <- names(SCENES)
cell_cols <- setdiff(names(cells), "cell_id")
all_cols <- names(panel)

audit <- data.table(
  metric = c(
    "rows", "columns_total", "unique_cells", "unique_scenes",
    "original_DNDC_columns", "cell_metadata_columns", "scene_metadata_columns",
    "missing_cell_metadata_cells", "missing_scene_metadata_scenes",
    "output_rds_MB", "output_fst_MB"
  ),
  value = as.character(c(
    nrow(panel),
    ncol(panel),
    uniqueN(panel$cell_id),
    uniqueN(panel$scene_id),
    sum(all_cols %in% original_cols),
    sum(all_cols %in% cell_cols),
    sum(all_cols %in% scene_cols),
    panel[is.na(lon) | is.na(lat), uniqueN(cell_id)],
    panel[is.na(scene_label), uniqueN(scene_id)],
    sz_rds,
    ifelse(is.na(sz_fst), "not_written", sz_fst)
  ))
)

scene_audit <- unique(
  panel[, intersect(
    c("scene_id", "scene_label", "scene_category",
      "fert_total_N_kgNha", "fert_urea_kgNha",
      "total_flood_days", "awd_flag",
      "till_applications", "till_method_label",
      "residue_summer_frac", "residue_monsoon_frac"),
    names(panel)
  ), with = FALSE],
  by = "scene_id"
)

sink(AUDIT_TXT)
cat("ENRICHED DAILY PANEL AUDIT\n")
cat(strrep("=", 72), "\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

cat("SUMMARY\n")
cat(strrep("-", 72), "\n")
print(audit)

cat("\nSCENE MANAGEMENT TABLE\n")
cat(strrep("-", 72), "\n")
print(scene_audit[order(scene_id)])

cat("\nCOLUMN INVENTORY\n")
cat(strrep("-", 72), "\n")
for (col in names(panel)) {
  src <- if (col %in% original_cols) "SIM" else if (col %in% scene_cols) "SCENE" else "CELL"
  cat(sprintf("%-6s %s\n", src, col))
}
sink()

# =============================================================================
# 8. SUMMARY
# =============================================================================

message(sprintf(
  "Script 8 complete: %s rows | %d columns | RDS %.1f MB | audit %s",
  format(nrow(panel), big.mark = ","), ncol(panel), sz_rds, basename(AUDIT_TXT)
))

if (USE_FST) {
  message(sprintf("FST written: %.1f MB -> %s", sz_fst, basename(OUT_FST)))
}

