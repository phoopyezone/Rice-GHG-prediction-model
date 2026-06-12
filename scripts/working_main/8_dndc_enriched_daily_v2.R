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