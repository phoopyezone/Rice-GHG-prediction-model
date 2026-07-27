###############################################################################
# Script 8 — Enrich Daily Panel with Cell & Scenario Metadata
# Myanmar Rice Paddy CH4 | KGML Pipeline
#
# WHAT THIS SCRIPT DOES
# Attaches two metadata tables to the 24.9M-row daily panel:
#
# (A) CELL lookup (cells.rds) joined by cell_id
#     Attaches: lat, lon, all soil properties (bulk_density, clay_fraction,
#               soil_pH, SOC_fraction, porosity, field_capacity, wilting_point,
#               hydro_conductivity, watertable_depth, soil_slope, etc.)
#
# (B) SCENE lookup (scene_registry.csv) joined by scene_id
#     The raw file has 132 DNDC parameter columns (F1..F6 fertilizer blocks,
#     FL1/FL2 flooding blocks, T1/T2 tillage, C1/C2 crop, CS1 system).
#     This script computes 24 clean KGML-ready management variables:
#       scene_label, scene_category,
#       fert_total_N_kgNha, fert_urea_kgNha, fert_ammonium_kgNha,
#       fert_phosphate_kgPha, fert_n_applications,
#       fert_method_code, fert_method_label,
#       residue_summer_frac, residue_monsoon_frac,
#       FL1_flood_days, FL2_flood_days, total_flood_days, awd_flag,
#       irrigation_control_code,
#       till_applications, till_method_code, till_method_label,
#       max_yield_summer_kgha,
#       planting_month_summer, harvest_month_summer, harvest_month_monsoon,
#       Site_name
#
# CONFIRMED STRUCTURE
#   cells.rds       : 263 rows, soil/spatial properties per grid cell
#   scene_registry  : 13 rows x 132 cols, Scenario_ID = S000..S013 (no S007)
#   daily panel     : 24.9M rows x 48 cols, keys = case_id/cell_id/scene_id
#
# MEMORY APPROACH
#   Both lookups are tiny. Data.table binary joins on sorted keys add new
#   columns in-place without duplicating the 24.9M-row panel. Each lookup
#   table is freed from RAM immediately after its join completes.
###############################################################################

library(data.table)

USE_FST <- requireNamespace("fst", quietly = TRUE)
if (USE_FST) library(fst)

# =============================================================================
# 0. PATHS  — only change these
# =============================================================================

PANEL_RDS  <- "C:/DNDC/simulate_dndc/dndc_outputs/daily_panel_allcases.rds"
CELLS_RDS  <- "C:/DNDC/simulate_dndc/dndc_inputs/input_files/processing/cells.rds"
SCENE_CSV  <- "C:/DNDC/simulate_dndc/dndc_inputs/input_files/processing/scene_registry.csv"
OUTPUT_DIR <- "C:/DNDC/simulate_dndc/dndc_outputs"

OUT_RDS   <- file.path(OUTPUT_DIR, "enriched_daily_panel.rds")
OUT_FST   <- file.path(OUTPUT_DIR, "enriched_daily_panel.fst")
AUDIT_TXT <- file.path(OUTPUT_DIR, "enrichment_audit.txt")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. CELL LOOKUP — cells.rds
# =============================================================================

message("Loading cells.rds ...")
cells_raw <- readRDS(CELLS_RDS)
if (!is.data.table(cells_raw)) cells_raw <- as.data.table(cells_raw)
message(sprintf("  %d rows x %d cols", nrow(cells_raw), ncol(cells_raw)))
message("  Columns: ", paste(names(cells_raw), collapse = ", "))

# Identify the cell_id column (may vary by Script 1 version)
CELL_ID_VARIANTS <- c("cell_id","cellid","id","cell","pixel_id","grid_id","cell_no")
cell_id_col <- NA_character_
for (v in CELL_ID_VARIANTS) {
  hit <- names(cells_raw)[tolower(names(cells_raw)) == v]
  if (length(hit)) { cell_id_col <- hit[1]; break }
}
if (is.na(cell_id_col)) {
  # Fall back to first integer column
  int_cols    <- names(cells_raw)[sapply(cells_raw, is.integer)]
  cell_id_col <- if (length(int_cols)) int_cols[1] else names(cells_raw)[1]
  message("  WARNING: defaulting cell_id to column '", cell_id_col,
          "'. Verify this is correct.")
}
message("  Cell ID column identified: '", cell_id_col, "'")
if (cell_id_col != "cell_id") setnames(cells_raw, cell_id_col, "cell_id")
cells_raw[, cell_id := as.integer(cell_id)]

# Rename known soil columns to canonical names regardless of Script 1 naming
# Keys = CANONICAL name, value = regex matching the actual column name in cells.rds
SOIL_RENAMES <- list(
  lat              = "^lat$|^latitude$|^Lat$|^Latitude$",
  lon              = "^lon$|^longitude$|^Lon$|^Longitude$",
  rice_area_ha     = "^rice_area_ha$|^spam_area$|^rice_area$|^area_ha$",
  bulk_density     = "^[Bb]ulk_density$|^bdod$",
  clay_fraction    = "^[Cc]lay_fraction$|^clay$",
  soil_pH          = "^pH$|^ph$|^soil_pH$|^phh2o$",
  SOC_fraction     = "^Top_layer_SOC$|^soc$|^soc_fraction$|^soc_top$",
  porosity         = "^[Pp]orosity$",
  field_capacity   = "^[Ff]ield_capacity$|^fc$",
  wilting_point    = "^[Ww]ilting_point$|^wp$",
  ksat_cmhr        = "^[Hh]ydro_conductivity$|^ksat$",
  watertable_depth = "^[Ww]atertable_depth$|^wtd$",
  soil_slope       = "^[Ss]oil_slope$|^slope$",
  soil_texture_id  = "^[Ss]oil_texture_[Ii][Dd]$|^texture_id$|^texture$",
  bypass_flow      = "^[Bb]ypass_flow$",
  soil_salinity    = "^[Ss]oil_salinity$",
  N_in_rainfall    = "^N_in_rainfall$|^n_dep$|^N_dep$",
  litter_fraction  = "^[Ll]itter_fraction$",
  humads_fraction  = "^[Hh]umads_fraction$",
  humus_fraction   = "^[Hh]umus_fraction$",
  humads_CN        = "^Humads_C.N$|^humads_cn$",
  humus_CN         = "^Humus_C.N$|^humus_cn$",
  SOC_profile_A    = "^SOC_profile_A$|^soc_profile_a$",
  SOC_profile_B    = "^SOC_profile_B$|^soc_profile_b$",
  init_NO3_ppm     = "^Initial_nitrate_ppm$|^init_no3$",
  init_NH4_ppm     = "^Initial_ammonium_ppm$|^init_nh4$"
)

# Apply renames where the actual column exists and differs from canonical
for (canon in names(SOIL_RENAMES)) {
  pattern <- SOIL_RENAMES[[canon]]
  hit     <- grep(pattern, names(cells_raw), value = TRUE)
  if (length(hit) && hit[1] != canon && hit[1] %in% names(cells_raw))
    setnames(cells_raw, hit[1], canon)
}

# Drop filler columns (named None/none/...)
drop_pat <- "^[Nn]one$|^[Nn]one_|^\\.none"
none_cols <- grep(drop_pat, names(cells_raw), value = TRUE)
if (length(none_cols)) cells_raw[, (none_cols) := NULL]

# Build final CELLS lookup (all columns except those that duplicate panel keys)
PANEL_KEYS <- c("case_id","scene_id","sim_year","cal_year","Day")
CELLS <- cells_raw[, setdiff(names(cells_raw), PANEL_KEYS), with = FALSE]

if (anyDuplicated(CELLS$cell_id)) {
  message("  Duplicate cell_ids detected — keeping first row per cell_id")
  CELLS <- CELLS[!duplicated(cell_id)]
}
setkeyv(CELLS, "cell_id")

message(sprintf("  CELLS lookup: %d cells, %d columns to attach",
                nrow(CELLS), ncol(CELLS) - 1L))
message("  Attached soil/spatial columns: ",
        paste(setdiff(names(CELLS), "cell_id"), collapse = ", "))

# =============================================================================
# 2. SCENE LOOKUP — scene_registry.csv
# =============================================================================

message("\nLoading scene_registry.csv ...")
scene_raw <- fread(SCENE_CSV, encoding = "UTF-8", na.strings = c("","NA","N/A"))

# Standardise the scenario ID column name
if ("Scenario_ID" %in% names(scene_raw))
  setnames(scene_raw, "Scenario_ID", "scene_id")
scene_raw[, scene_id := as.character(scene_id)]

message(sprintf("  %d rows x %d cols | Scenes: %s",
                nrow(scene_raw), ncol(scene_raw),
                paste(sort(scene_raw$scene_id), collapse=" ")))

# ── Scene descriptor lookup vectors ──────────────────────────────────────────
SCENE_LABELS <- c(
  S000="Baseline",      S001="CommInjection", S002="MedHighFert",
  S003="ExtremeFert",   S004="SplitApp",      S005="ShortFlood",
  S006="AWD",           S008="NoTill",        S009="DeepPlow",
  S010="ContIrrig",     S011="Rainfed",       S012="LowResidue",
  S013="HighResidue"
)
SCENE_CATS <- c(
  S000="baseline",
  S001="fertilizer", S002="fertilizer", S003="fertilizer", S004="fertilizer",
  S005="water",      S006="water",
  S008="tillage",    S009="tillage",
  S010="irrigation", S011="irrigation",
  S012="residue",    S013="residue"
)

# ── Helper: sum specific columns handling NAs ─────────────────────────────────
rc_sum <- function(dt, pattern) {
  cols <- grep(pattern, names(dt), perl = TRUE, value = TRUE)
  if (!length(cols)) return(rep(0.0, nrow(dt)))
  rowSums(dt[, cols, with = FALSE], na.rm = TRUE)
}

# ── Fertilizer totals ─────────────────────────────────────────────────────────
# Urea: F2_Urea through F6_Urea (F1 has no Urea application in this setup)
# Ammonium: F1_Ammonium, F1_Ammonium_02, F3_Ammonium, F3_Ammonium_02,
#            F4_Ammonium, F5_Ammonium, F6_Ammonium
#           (NOT Ammonium_bicarbonate — different fertilizer form)
# Phosphate: F1_Phosphate, F3_Phosphate, F4_Phosphate (not in F5/F6 here)

urea_total <- rc_sum(scene_raw, "^F[0-9]+_Urea$")
amm_total  <- rc_sum(scene_raw, "^F[0-9]+_Ammonium(_02)?$")
phos_total <- rc_sum(scene_raw, "^F[0-9]+_Phosphate$")
total_N    <- urea_total + amm_total

# ── Fertilizer method (from F1 application — primary application block) ───────
fert_method_code <- as.integer(
  fifelse(is.na(scene_raw$F1_Fertilizing_method), 0L,
          as.integer(scene_raw$F1_Fertilizing_method))
)
fert_method_label <- fcase(
  fert_method_code == 0L, "broadcast",
  fert_method_code == 1L, "banding",
  fert_method_code == 2L, "injection",
  default = "broadcast"
)
n_fert_apps <- as.integer(
  fifelse(is.na(scene_raw$CS1_Fertilizer_applications), 0L,
          as.integer(scene_raw$CS1_Fertilizer_applications))
)

# ── Residue fractions ─────────────────────────────────────────────────────────
residue_summer  <- scene_raw$C1_Residue_left_in_field   # summer (Crop 1)
residue_monsoon <- scene_raw$C2_Residue_left_in_field   # monsoon (Crop 2)

# ── Flood duration (days) ─────────────────────────────────────────────────────
# FL1: all four fields present (Start_month, Start_day, End_month, End_day)
# FL2: Start_month absent — use FL1_End_month as FL2 start month
fl1_days <- pmax(0L,
                 (as.integer(scene_raw$FL1_End_month)   - as.integer(scene_raw$FL1_Start_month)) * 30L +
                   (as.integer(scene_raw$FL1_End_day)     - as.integer(scene_raw$FL1_Start_day))
)

fl2_start_month <- as.integer(scene_raw$FL1_End_month)  # implicit from FL1 end
fl2_days <- pmax(0L,
                 (as.integer(scene_raw$FL2_End_month)   - fl2_start_month) * 30L +
                   (as.integer(scene_raw$FL2_End_day)     - as.integer(scene_raw$FL2_Start_day))
)

total_flood_days <- fl1_days + fl2_days

# ── AWD flag ──────────────────────────────────────────────────────────────────
awd_flag <- as.integer(
  (as.integer(scene_raw$FL1_Alter_wet_dry) == 1L) |
    (as.integer(scene_raw$FL2_Alter_wet_dry) == 1L)
)

# ── Tillage ───────────────────────────────────────────────────────────────────
till_apps <- as.integer(
  fifelse(is.na(scene_raw$CS1_Till_applications), 0L,
          as.integer(scene_raw$CS1_Till_applications))
)
till_method_code <- as.integer(
  fifelse(is.na(scene_raw$T1_Till_method), 0L,
          as.integer(scene_raw$T1_Till_method))
)
till_method_label <- fcase(
  till_method_code == 0L, "none",
  till_method_code == 1L, "no_till",
  till_method_code == 2L, "ridge",
  till_method_code == 3L, "chisel",
  till_method_code == 4L, "disk",
  till_method_code == 5L, "moldboard",
  default = "none"
)

# ── Irrigation control ────────────────────────────────────────────────────────
irrig_control <- as.integer(
  fifelse(is.na(scene_raw$CS1_Irrigation_control), 0L,
          as.integer(scene_raw$CS1_Irrigation_control))
)

# ── Crop parameters ───────────────────────────────────────────────────────────
max_yield_summer      <- scene_raw$C1_Maximum_yield
planting_month_summer <- as.integer(scene_raw$C1_Planting_month)
harvest_month_summer  <- as.integer(scene_raw$C1_Harvest_month)
# Monsoon harvest: FL2 ends Oct (month 10), harvest ~Nov (month 11)
harvest_month_monsoon <- as.integer(scene_raw$FL2_End_month) + 1L

# ── Assemble final SCENES lookup table ───────────────────────────────────────
SCENES <- data.table(
  scene_id              = scene_raw$scene_id,
  scene_label           = SCENE_LABELS[scene_raw$scene_id],
  scene_category        = SCENE_CATS[scene_raw$scene_id],
  Site_name             = scene_raw$Site_name,
  # Fertilizer
  fert_total_N_kgNha    = round(total_N,    2),
  fert_urea_kgNha       = round(urea_total, 2),
  fert_ammonium_kgNha   = round(amm_total,  2),
  fert_phosphate_kgPha  = round(phos_total, 2),
  fert_n_applications   = n_fert_apps,
  fert_method_code      = fert_method_code,
  fert_method_label     = fert_method_label,
  # Residue
  residue_summer_frac   = residue_summer,
  residue_monsoon_frac  = residue_monsoon,
  # Water / flooding
  FL1_flood_days        = fl1_days,
  FL2_flood_days        = fl2_days,
  total_flood_days      = total_flood_days,
  awd_flag              = awd_flag,
  irrigation_control_code = irrig_control,
  # Tillage
  till_applications     = till_apps,
  till_method_code      = till_method_code,
  till_method_label     = till_method_label,
  # Crop
  max_yield_summer_kgha = max_yield_summer,
  planting_month_summer = planting_month_summer,
  harvest_month_summer  = harvest_month_summer,
  harvest_month_monsoon = harvest_month_monsoon
)
setkeyv(SCENES, "scene_id")

message(sprintf("  SCENES lookup: %d scenarios, %d management columns computed",
                nrow(SCENES), ncol(SCENES) - 1L))
message("\n  Scene management summary (verification):")
print(SCENES[, .(scene_id, scene_label, fert_total_N_kgNha,
                 total_flood_days, awd_flag, till_method_label,
                 residue_summer_frac, irrigation_control_code)])

# =============================================================================
# 3. LOAD DAILY PANEL
# =============================================================================

message("\nLoading daily panel ...")
t0 <- proc.time()
panel <- readRDS(PANEL_RDS)
t_load <- round((proc.time()-t0)[["elapsed"]], 1)
message(sprintf("  %s rows x %d cols loaded in %.1f s",
                format(nrow(panel), big.mark=","), ncol(panel), t_load))
message("  Memory: ~", format(object.size(panel), units="GB"))

panel[, cell_id  := as.integer(cell_id)]
panel[, scene_id := as.character(scene_id)]

# Coverage check
n_cells_missing  <- length(setdiff(unique(panel$cell_id),  CELLS$cell_id))
n_scenes_missing <- length(setdiff(unique(panel$scene_id), SCENES$scene_id))
if (n_cells_missing  > 0) warning(n_cells_missing,  " cell_ids in panel not in cells.rds")
if (n_scenes_missing > 0) warning(n_scenes_missing, " scene_ids in panel not in scene_registry")

# =============================================================================
# 4. JOIN BOTH LOOKUPS
# =============================================================================

message("\nJoining cell metadata (soil/spatial) by cell_id ...")
t1 <- proc.time()
setkeyv(panel, "cell_id")
panel <- merge(panel, CELLS, by = "cell_id", all.x = TRUE, sort = FALSE)
rm(CELLS, cells_raw); gc()
message(sprintf("  Done %.1f s — panel now %d cols", (proc.time()-t1)[3], ncol(panel)))

message("Joining scene metadata (management) by scene_id ...")
t2 <- proc.time()
setkeyv(panel, "scene_id")
panel <- merge(panel, SCENES, by = "scene_id", all.x = TRUE, sort = FALSE)
rm(SCENES, scene_raw); gc()
message(sprintf("  Done %.1f s — panel now %d cols", (proc.time()-t2)[3], ncol(panel)))

# Restore sort order
message("Sorting ...")
sk <- c("scene_id","case_id","cell_id","sim_year","cal_year","Day")
setkeyv(panel, sk[sk %in% names(panel)])

message(sprintf("Final enriched panel: %s rows x %d cols",
                format(nrow(panel), big.mark=","), ncol(panel)))

# Quick spot-check
chk_cols <- c("lat","lon","bulk_density","clay_fraction","soil_pH",
              "SOC_fraction","fert_total_N_kgNha","total_flood_days",
              "awd_flag","till_method_label","residue_summer_frac",
              "scene_label","scene_category")
chk_cols <- chk_cols[chk_cols %in% names(panel)]
message("\nSpot-check (row 1, new columns):")
print(panel[1, chk_cols, with = FALSE])
message("\nNA counts in new columns:")
for (col in chk_cols) {
  n_na <- sum(is.na(panel[[col]]))
  message(sprintf("  %-30s %s NAs", col, format(n_na, big.mark=",")))
}

# =============================================================================
# 5. SAVE
# =============================================================================

message("\nSaving enriched_daily_panel.rds ...")
t3 <- proc.time()
saveRDS(panel, OUT_RDS, compress = "gzip")
t_rds  <- round((proc.time()-t3)[3], 1)
sz_rds <- round(file.info(OUT_RDS)$size / 1024^2, 1)
message(sprintf("  Saved: %.1f MB in %.1f s -> %s", sz_rds, t_rds, OUT_RDS))

if (USE_FST) {
  message("Saving enriched_daily_panel.fst (fast column reads) ...")
  t4 <- proc.time()
  write_fst(as.data.frame(panel), OUT_FST, compress = 50)
  t_fst  <- round((proc.time()-t4)[3], 1)
  sz_fst <- round(file.info(OUT_FST)$size / 1024^2, 1)
  message(sprintf("  Saved: %.1f MB in %.1f s -> %s", sz_fst, t_fst, OUT_FST))
}

# =============================================================================
# 6. AUDIT FILE
# =============================================================================

original_cols <- c("case_id","cell_id","scene_id","sim_year","cal_year","Day",
                   "SOC","dSOC","DOC","Microbe","Humads","Humus","NPP","NEE",
                   "Photosynthesis","Soil-heterotrophic-respiration",
                   "CH4-prod.","CH4-oxid.","CH4-flux","CH4-pool","CH4-DOC",
                   "Litter-C","Manure-C","DOC-leach",
                   "SoilT_1cm","SoilT_5cm","SoilT_10cm","SoilT_20cm",
                   "SoilM_1cm","SoilM_5cm","SoilM_10cm","SoilM_20cm",
                   "SoilEh_1cm","SoilEh_10cm","SoilEh_20cm",
                   "SoilWater_mm","DeepWater_mm","Soil_pH_1cm",
                   "Ponding","Precipitation","Irrigation","Evaporation",
                   "Transpiration","Leaching","Runoff","IniSoilWater","EndSoilWater",
                   "LeafC","StemC","RootC","GrainC","TDD","GrowthIndex",
                   "Water_stress","N_stress","LAI","TotalCropN","DailyCropGrowth","DayGrainGrowth")

scene_cols <- c("scene_label","scene_category","Site_name",
                "fert_total_N_kgNha","fert_urea_kgNha","fert_ammonium_kgNha",
                "fert_phosphate_kgPha","fert_n_applications","fert_method_code",
                "fert_method_label","residue_summer_frac","residue_monsoon_frac",
                "FL1_flood_days","FL2_flood_days","total_flood_days","awd_flag",
                "irrigation_control_code","till_applications","till_method_code",
                "till_method_label","max_yield_summer_kgha","planting_month_summer",
                "harvest_month_summer","harvest_month_monsoon")

all_cols  <- names(panel)
cell_cols <- setdiff(all_cols, c(original_cols, scene_cols))

sink(AUDIT_TXT)
cat("ENRICHED DAILY PANEL AUDIT\n")
cat(strrep("=", 72), "\n")
cat(sprintf("Generated       : %s\n", Sys.time()))
cat(sprintf("Rows            : %s\n", format(nrow(panel), big.mark=",")))
cat(sprintf("Columns (total) : %d\n", length(all_cols)))
cat(sprintf("  SIM  (original DNDC outputs) : %d\n",
            sum(all_cols %in% original_cols)))
cat(sprintf("  CELL (from cells.rds)        : %d\n", length(cell_cols)))
cat(sprintf("  SCENE (from scene_registry)  : %d\n",
            sum(all_cols %in% scene_cols)))
cat(sprintf("Output RDS : %s  (%.1f MB)\n", basename(OUT_RDS), sz_rds))
if (USE_FST) cat(sprintf("Output FST : %s  (%.1f MB)\n", basename(OUT_FST), sz_fst))
cat("\n")

cat("SCENE MANAGEMENT TABLE (derived variables per scenario)\n")
cat(strrep("-", 72), "\n")
cat(sprintf("%-6s  %-14s  %-8s  %-8s  %-5s  %5s  %5s  %-8s  %-6s  %-6s\n",
            "Scene","Label","N_total","Urea","Flood","AWD","Till","TillMeth","Res_S","Res_M"))
for (r in seq_len(nrow(panel[, .N, by=scene_id]))) {
  sid   <- sort(unique(panel$scene_id))[r]
  first <- panel[scene_id == sid][1L]
  cat(sprintf("%-6s  %-14s  %8.1f  %8.1f  %5d  %5d  %5d  %-8s  %6.2f  %6.2f\n",
              sid,
              ifelse("scene_label" %in% names(first), first$scene_label, "?"),
              ifelse("fert_total_N_kgNha" %in% names(first), first$fert_total_N_kgNha, NA_real_),
              ifelse("fert_urea_kgNha"    %in% names(first), first$fert_urea_kgNha,    NA_real_),
              ifelse("total_flood_days"   %in% names(first), first$total_flood_days,    0L),
              ifelse("awd_flag"           %in% names(first), first$awd_flag,            0L),
              ifelse("till_applications"  %in% names(first), first$till_applications,   0L),
              ifelse("till_method_label"  %in% names(first), first$till_method_label,  "?"),
              ifelse("residue_summer_frac" %in% names(first), first$residue_summer_frac, NA_real_),
              ifelse("residue_monsoon_frac" %in% names(first), first$residue_monsoon_frac, NA_real_)
  ))
}

cat("\nFULL COLUMN INVENTORY\n")
cat(strrep("-", 72), "\n")
for (col in all_cols) {
  src <- if (col %in% original_cols) "SIM  " else
    if (col %in% scene_cols)    "SCENE" else "CELL "
  x   <- panel[[col]]
  nas <- sum(is.na(x))
  if (is.numeric(x))
    cat(sprintf("  [%s] %-35s num  NA=%-8d  range=[%.4g, %.4g]\n",
                src, col, nas, min(x,na.rm=TRUE), max(x,na.rm=TRUE)))
  else
    cat(sprintf("  [%s] %-35s chr  NA=%-8d  ex: %s\n",
                src, col, nas,
                paste(head(unique(x[!is.na(x)]), 3), collapse=",")))
}
sink()
message("Audit written: ", AUDIT_TXT)

# =============================================================================
# 7. FINAL SUMMARY
# =============================================================================

message("\n", strrep("=", 70))
message("ENRICHMENT COMPLETE")
message(strrep("=", 70))
message(sprintf("  Rows          : %s", format(nrow(panel), big.mark=",")))
message(sprintf("  Columns total : %d", ncol(panel)))
message(sprintf("    SIM (DNDC outputs)   : %d", sum(all_cols %in% original_cols)))
message(sprintf("    CELL (soil/spatial)  : %d", length(cell_cols)))
message(sprintf("    SCENE (management)   : %d", sum(all_cols %in% scene_cols)))
message(sprintf("  Output RDS    : %.1f MB  ->  %s", sz_rds, basename(OUT_RDS)))
if (USE_FST)
  message(sprintf("  Output FST    : %.1f MB  ->  %s  (fast column reads)", sz_fst, basename(OUT_FST)))
message(sprintf("  Audit         : %s", basename(AUDIT_TXT)))
message("")
message("To use the enriched panel:")
message("  panel <- readRDS('", OUT_RDS, "')")
if (USE_FST) {
  message("  # Or fast column-selective read:")
  message("  cols <- c('scene_id','cell_id','sim_year','cal_year','Day',")
  message("            'CH4-flux','lat','lon','clay_fraction','soil_pH',")
  message("            'fert_total_N_kgNha','total_flood_days','awd_flag',")
  message("            'scene_label','scene_category','till_method_label')")
  message("  sub <- fst::read_fst('", OUT_FST, "', columns=cols)")
}