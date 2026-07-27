###############################################################################
# Script 5 — DNDC Batch Output Extraction
# Myanmar Rice Paddy CH4 | KGML Pipeline
#
# WHAT THIS SCRIPT DOES
# ─────────────────────
# 1. Reads the batch .txt file to build an exact case_id → cell_id + scene_id
#    lookup table (no formula — parsed directly from the .dnd filenames).
# 2. Scans the DNDC batch output folder for Case{N}-Myanmar_Rice_Cell_{M}
#    directories and matches them to the lookup.
# 3. Extracts four daily output file types per case × year:
#      Day_SoilC        CH4 flux/prod/oxid, SOC pools, NEE, NPP
#      Day_SoilClimate  soil T, moisture, Eh, pH profiles
#      Day_SoilWater    ponding, precip, ET, runoff
#      Day_FieldCrop    biomass, LAI, water/N stress
# 4. Extracts the annual Multi_year_summary.csv per case.
# 5. Saves ONE RDS per file type (all cases combined):
#      annual_allcases.rds / .csv
#      daily_SoilC_allcases.rds
#      daily_SoilClimate_allcases.rds
#      daily_SoilWater_allcases.rds
#      daily_FieldCrop_allcases.rds
#    Plus:
#      case_scene_index.csv   authoritative case_id/cell_id/scene_id lookup
#
# IDENTIFIERS ON EVERY ROW
# ─────────────────────────
#   case_id   integer   from folder "Case{N}-Myanmar_Rice_Cell_{M}"
#   cell_id   integer   spatial grid cell
#   scene_id  character "S000" ... "S013" parsed from batch file
#   sim_year  integer   1–20 (DNDC internal year index)
#   cal_year  integer   2005–2024 (sim_year + START_YEAR - 1)
#
# CONFIRMED BATCH STRUCTURE (myanmar_rice_batch_all_rerun.txt)
#   3,419 cases | 263 cells | 13 scenes (S000–S013, no S007)
#   Ordering: scenes outer loop, cells inner loop (263 cells per scene block)
#
# MEMORY APPROACH
# ───────────────
# Each file type is combined, saved, and freed before the next one loads.
# No cross-file join is attempted here. Script 6 joins per-case during
# seasonal aggregation (~7,300 rows at a time — no OOM risk).
###############################################################################

library(data.table)
library(stringr)
library(future)
library(future.apply)

# =============================================================================
# 0. CONFIGURATION  — edit paths and N_WORKERS only
# =============================================================================

BATCH_TXT  <- "C:/DNDC/simulate_dndc/dndc_inputs/input_files/template/myanmar_rice_batch_all_rerun.txt"
BATCH_DIR  <- "C:/DNDC/Result/Record/Batch"
OUTPUT_DIR <- "C:/DNDC/simulate_dndc/dndc_outputs"

# Path verification block — auto-finds the batch file if the path above is wrong.
# Check the correct path from the error message and update BATCH_TXT above.
if (!file.exists(BATCH_TXT)) {
  candidates <- c(
    BATCH_TXT,
    "C:/DNDC/simulate_dndc/myanmar_rice_batch_all_rerun.txt",
    "C:/DNDC/simulate_dndc/dndc_inputs/myanmar_rice_batch_all_rerun.txt",
    "C:/DNDC/simulate_dndc/dndc_inputs/input_files/myanmar_rice_batch_all_rerun.txt",
    "C:/DNDC/simulate_dndc/dndc_inputs/input_files/template/myanmar_rice_batch_all_rerun.txt"
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) {
    message("BATCH_TXT path updated automatically to: ", found[1])
    BATCH_TXT <- found[1]
  } else {
    hits <- tryCatch(
      list.files("C:/DNDC", pattern = "myanmar_rice_batch_all_rerun\\.txt",
                 recursive = TRUE, full.names = TRUE),
      error = function(e) character(0)
    )
    if (length(hits) > 0)
      stop("Batch file found at: ", hits[1],
           "\nUpdate BATCH_TXT at the top of this script and re-run.")
    stop("Cannot locate myanmar_rice_batch_all_rerun.txt.\n",
         "Set BATCH_TXT to the correct full path and re-run.")
  }
}
if (!dir.exists(BATCH_DIR))
  stop("BATCH_DIR not found: ", BATCH_DIR)

message("Paths confirmed:")
message("  BATCH_TXT  : ", BATCH_TXT)
message("  BATCH_DIR  : ", BATCH_DIR)
message("  OUTPUT_DIR : ", OUTPUT_DIR)

N_YEARS    <- 20L    # simulated years per case
START_YEAR <- 2005L  # first calendar year of simulation
N_WORKERS  <- 4L     # parallel workers; set to 2 if RAM < 16 GB

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. BUILD EXACT CASE -> CELL + SCENE LOOKUP FROM BATCH FILE
# =============================================================================

message("Reading batch file: ", BATCH_TXT)

raw_lines <- tryCatch(
  readLines(BATCH_TXT, warn = FALSE, encoding = "UTF-8"),
  error = function(e) readLines(BATCH_TXT, warn = FALSE)
)
raw_lines <- trimws(gsub("\r", "", raw_lines))

# Line 1 = total count of cases; lines 2+ = one .dnd path per case
n_declared <- suppressWarnings(as.integer(raw_lines[1]))
dnd_lines  <- raw_lines[raw_lines != "" & !grepl("^\\d+$", raw_lines)][-1]
# Simpler: just take all non-empty lines after line 1
dnd_lines  <- raw_lines[-1]
dnd_lines  <- dnd_lines[nchar(dnd_lines) > 5]

# Parse Cell_{cell_id}_{scene_id}.dnd from each path
LOOKUP <- data.table(
  case_id  = seq_along(dnd_lines),
  dnd_path = dnd_lines
)
LOOKUP[, cell_id  := as.integer(str_extract(dnd_path, "(?<=Cell_)\\d+(?=_S)"))]
LOOKUP[, scene_id := str_extract(dnd_path, "S\\d{3}(?=\\.dnd)")]

# Remove any lines that did not parse correctly
n_bad <- LOOKUP[is.na(cell_id) | is.na(scene_id), .N]
if (n_bad > 0) message(sprintf("WARNING: %d lines could not be parsed — skipping", n_bad))
LOOKUP <- LOOKUP[!is.na(cell_id) & !is.na(scene_id)]

message(sprintf(
  "Batch file: %d declared | %d parsed | %d cells | %d scenes",
  n_declared, nrow(LOOKUP), uniqueN(LOOKUP$cell_id), uniqueN(LOOKUP$scene_id)
))
message("Records per scene:")
print(LOOKUP[, .N, by = scene_id][order(scene_id)])

fwrite(LOOKUP[, .(case_id, cell_id, scene_id)],
       file.path(OUTPUT_DIR, "case_scene_index.csv"))
message("Saved: case_scene_index.csv\n")

# =============================================================================
# 2. MATCH OUTPUT FOLDERS TO LOOKUP
# =============================================================================

message("Scanning: ", BATCH_DIR)

all_dirs   <- list.dirs(BATCH_DIR, full.names = TRUE, recursive = FALSE)
FOLD_PAT   <- "^Case(\\d+)-Myanmar_Rice_Cell_(\\d+)$"
valid_dirs <- all_dirs[grepl(FOLD_PAT, basename(all_dirs))]

if (length(valid_dirs) == 0)
  stop("No matching folders found in BATCH_DIR: ", BATCH_DIR,
       "\nExpected pattern: Case{N}-Myanmar_Rice_Cell_{M}")

folder_dt <- data.table(
  folder  = valid_dirs,
  case_id = as.integer(str_extract(basename(valid_dirs), "(?<=Case)\\d+"))
)

# Join cell_id and scene_id from the batch file lookup (ground truth)
case_index <- merge(folder_dt,
                    LOOKUP[, .(case_id, cell_id, scene_id)],
                    by = "case_id", all.x = TRUE)

n_unmatched <- case_index[is.na(scene_id), .N]
if (n_unmatched > 0)
  message(sprintf(
    "WARNING: %d folders have no batch-file match and will be skipped", n_unmatched
  ))

case_index <- case_index[!is.na(scene_id)]
setorder(case_index, case_id)

message(sprintf(
  "Matched: %d folders | %d cells | %d scenes\n",
  nrow(case_index), uniqueN(case_index$cell_id), uniqueN(case_index$scene_id)
))

# =============================================================================
# 3. COLUMN DEFINITIONS
# =============================================================================

# Annual Multi_year_summary.csv — no header row; 52 positional columns per
# DNDC 9.5 manual. Columns beyond 52 are named V53, V54, …
ANNUAL_COLS <- c(
  "Year",
  "GrainC1",  "LeafStemC1", "RootC1",
  "GrainC2",  "LeafStemC2", "RootC2",
  "GrainC3",  "LeafStemC3", "RootC3",
  "SOC_0_10cm", "SOC_0_20cm", "SOC_0_30cm",
  "Ini_SOC",  "End_SOC",    "dSOC",
  "LitterC_input", "RootC_input", "ManureC_input",
  "Soil_CO2", "CH4",
  "Ini_SON",  "Ini_SIN",    "End_SON",  "End_SIN", "dSN",
  "Atmo_N_input", "Fertilizer_N_input", "Manure_N_input",
  "Litter_N_input", "N_fixation",
  "Crop_N_uptake", "N_leach", "N_runoff",
  "N2O_flux", "NO_flux",    "N2_flux",  "NH3_flux", "Exch_NH4",
  "PET",      "Transpiration", "Evaporation", "WaterLeach", "Runoff",
  "Irrigation", "Precipitation",
  "MeanT",    "WindSpeed",
  "ColdStress", "WaterStress", "N_Stress",
  "Cut_CropC"
)

# Day_SoilC — skip=1, standard CSV header on row 1
KEEP_SoilC <- c(
  "Day",
  "SOC", "dSOC", "DOC", "Microbe", "Humads", "Humus",
  "NPP", "NEE", "Photosynthesis", "Soil-heterotrophic-respiration",
  "CH4-DOC", "CH4-prod.", "CH4-oxid.", "CH4-flux", "CH4-pool",
  "Litter-C", "Manure-C", "DOC-leach"
)

# Day_SoilClimate — position-based (fread deduplicates repeated depth headers)
SOILCLIM_SRC <- c(
  "Day",
  "1", "5", "10", "20",           # soil temperature (C) at 1,5,10,20 cm
  "1.1", "5.1", "10.1", "20.1",   # volumetric water content
  "1.3", "10.3", "20.3",           # Eh redox potential (mV)
  "(mm).2", "(mm).3",              # SoilWater_mm, DeepWater_mm
  "1.5"                            # soil pH at 1 cm
)
SOILCLIM_DST <- c(
  "Day",
  "SoilT_1cm",  "SoilT_5cm",  "SoilT_10cm",  "SoilT_20cm",
  "SoilM_1cm",  "SoilM_5cm",  "SoilM_10cm",  "SoilM_20cm",
  "SoilEh_1cm", "SoilEh_10cm","SoilEh_20cm",
  "SoilWater_mm", "DeepWater_mm",
  "Soil_pH_1cm"
)

# Day_SoilWater — skip=1, standard CSV header on row 1
KEEP_SoilWater <- c(
  "Day",
  "Ponding", "Precipitation", "Irrigation",
  "Evaporation", "Transpiration", "Leaching", "Runoff",
  "IniSoilWater", "EndSoilWater"
)

# Day_FieldCrop — two-row merged header; data starts at row 5 (skip=4)
KEEP_FieldCrop <- c(
  "Day",
  "LeafC", "StemC", "RootC", "GrainC",
  "TDD", "GrowthIndex",
  "Water_stress", "N_stress", "LAI",
  "TotalCropN", "DailyCropGrowth", "DayGrainGrowth"
)

# =============================================================================
# 4. FILE READER FUNCTIONS
# =============================================================================

# Drop non-numeric sentinel rows that appear in some DNDC CSVs
drop_non_numeric <- function(dt) {
  dt[suppressWarnings(!is.na(as.numeric(as.character(dt[[1]]))))]
}

# Tag rows with identifiers and year columns
tag <- function(dt, case_id, cell_id, scene_id, yr) {
  dt[, Day      := as.integer(Day)]
  dt[, sim_year := yr]
  dt[, cal_year := yr + START_YEAR - 1L]
  dt[, case_id  := case_id]
  dt[, cell_id  := cell_id]
  dt[, scene_id := scene_id]
  dt
}

# ── Annual summary ────────────────────────────────────────────────────────────
read_annual <- function(path, case_id, cell_id, scene_id) {
  dt <- tryCatch(
    fread(path, header = FALSE, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  n_named <- min(ncol(dt), length(ANNUAL_COLS))
  extras  <- if (ncol(dt) > n_named)
    paste0("V", seq(n_named + 1L, ncol(dt))) else character(0L)
  setnames(dt, c(ANNUAL_COLS[seq_len(n_named)], extras))
  dt[, sim_year := as.integer(Year)]
  dt[, cal_year := sim_year + START_YEAR - 1L]
  dt[, case_id  := case_id]
  dt[, cell_id  := cell_id]
  dt[, scene_id := scene_id]
  dt
}

# ── Day_SoilC and Day_SoilWater  (skip=1, clean header) ──────────────────────
read_skip1 <- function(path, keep_cols, case_id, cell_id, scene_id, yr) {
  dt <- tryCatch(
    fread(path, skip = 1L, header = TRUE, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  dt <- drop_non_numeric(dt)
  if (nrow(dt) == 0L) return(NULL)
  wanted <- keep_cols[keep_cols %in% names(dt)]
  if (!"Day" %in% wanted) return(NULL)
  tag(dt[, wanted, with = FALSE], case_id, cell_id, scene_id, yr)
}

# ── Day_SoilClimate  (position-based extraction) ──────────────────────────────
read_soilclimate <- function(path, case_id, cell_id, scene_id, yr) {
  lns <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(lns) || length(lns) == 0L) return(NULL)
  hdr_i <- which(grepl("^Day,", lns))[1L]
  if (is.na(hdr_i)) return(NULL)
  data_lns <- lns[(hdr_i + 1L):length(lns)]
  data_lns <- data_lns[nchar(trimws(data_lns)) > 0L]
  if (length(data_lns) == 0L) return(NULL)
  dt <- tryCatch(
    fread(text   = paste(c(lns[hdr_i], data_lns), collapse = "\n"),
          header = TRUE, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  dt <- drop_non_numeric(dt)
  if (nrow(dt) == 0L) return(NULL)
  present <- SOILCLIM_SRC %in% names(dt)
  if (!any(present)) return(NULL)
  out <- tryCatch(dt[, SOILCLIM_SRC[present], with = FALSE],
                  error = function(e) NULL)
  if (is.null(out) || !"Day" %in% SOILCLIM_DST[present]) return(NULL)
  setnames(out, SOILCLIM_DST[present])
  tag(out, case_id, cell_id, scene_id, yr)
}

# ── Day_FieldCrop  (merged two-row header; skip=4) ────────────────────────────
read_fieldcrop <- function(path, case_id, cell_id, scene_id, yr) {
  lns <- tryCatch(readLines(path, n = 5L), error = function(e) NULL)
  if (is.null(lns) || length(lns) < 3L) return(NULL)
  r1  <- trimws(strsplit(lns[2L], ",")[[1L]])
  r2  <- trimws(strsplit(lns[3L], ",")[[1L]])
  n   <- max(length(r1), length(r2))
  r1  <- c(r1, rep("", n - length(r1)))
  r2  <- c(r2, rep("", n - length(r2)))
  hdr <- ifelse(r2 != "", r2, r1)
  dt  <- tryCatch(
    fread(path, skip = 4L, header = FALSE, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  setnames(dt,
           seq_len(min(ncol(dt), length(hdr))),
           hdr[seq_len(min(ncol(dt), length(hdr)))])
  dt <- drop_non_numeric(dt)
  if (nrow(dt) == 0L) return(NULL)
  wanted <- KEEP_FieldCrop[KEEP_FieldCrop %in% names(dt)]
  if (!"Day" %in% wanted) return(NULL)
  tag(dt[, wanted, with = FALSE], case_id, cell_id, scene_id, yr)
}

# =============================================================================
# 5. PER-CASE EXTRACTION  (runs inside each parallel worker)
# =============================================================================

extract_one_case <- function(folder, case_id, cell_id, scene_id, n_years) {
  
  # Read all years for one file type and combine
  all_years <- function(reader_fn, prefix) {
    out <- vector("list", n_years)
    for (yr in seq_len(n_years)) {
      p <- file.path(folder, sprintf("%s_%d.csv", prefix, yr))
      if (!file.exists(p)) next
      out[[yr]] <- tryCatch(
        reader_fn(p, case_id, cell_id, scene_id, yr),
        error = function(e) NULL
      )
    }
    out <- Filter(Negate(is.null), out)
    if (length(out) == 0L) return(NULL)
    rbindlist(out, fill = TRUE, use.names = TRUE)
  }
  
  ann_path <- file.path(folder, "Multi_year_summary.csv")
  
  list(
    Annual = if (file.exists(ann_path))
      tryCatch(read_annual(ann_path, case_id, cell_id, scene_id),
               error = function(e) NULL)
    else NULL,
    
    SoilC = all_years(
      function(p, ci, cid, sid, yr) read_skip1(p, KEEP_SoilC, ci, cid, sid, yr),
      "Day_SoilC"
    ),
    
    SoilClimate = all_years(read_soilclimate, "Day_SoilClimate"),
    
    SoilWater = all_years(
      function(p, ci, cid, sid, yr) read_skip1(p, KEEP_SoilWater, ci, cid, sid, yr),
      "Day_SoilWater"
    ),
    
    FieldCrop = all_years(read_fieldcrop, "Day_FieldCrop")
  )
}

# =============================================================================
# 6. PARALLEL EXTRACTION
# =============================================================================

plan(multisession, workers = N_WORKERS)
message(sprintf("Extracting %d cases on %d workers ...", nrow(case_index), N_WORKERS))
t0 <- proc.time()

all_results <- future_mapply(
  FUN         = extract_one_case,
  folder      = case_index$folder,
  case_id     = case_index$case_id,
  cell_id     = case_index$cell_id,
  scene_id    = case_index$scene_id,
  MoreArgs    = list(n_years = N_YEARS),
  SIMPLIFY    = FALSE,
  future.seed = TRUE
)

plan(sequential)
elapsed <- (proc.time() - t0)[["elapsed"]]
message(sprintf("Done in %.1f minutes.\n", elapsed / 60))

# =============================================================================
# 7. COMBINE AND SAVE — one file type at a time
# =============================================================================

FILE_TYPES <- c("Annual", "SoilC", "SoilClimate", "SoilWater", "FieldCrop")

SORT_KEYS <- list(
  Annual      = c("scene_id", "case_id", "cell_id", "sim_year"),
  SoilC       = c("scene_id", "case_id", "cell_id", "sim_year", "Day"),
  SoilClimate = c("scene_id", "case_id", "cell_id", "sim_year", "Day"),
  SoilWater   = c("scene_id", "case_id", "cell_id", "sim_year", "Day"),
  FieldCrop   = c("scene_id", "case_id", "cell_id", "sim_year", "Day")
)

message("Combining and saving (one file type at a time) ...")
inventory <- data.table(
  file_type = FILE_TYPES, rows = 0L, cols = 0L,
  size_mb = 0.0, status = "empty"
)

for (ft in FILE_TYPES) {
  
  message(sprintf("  %-14s ", ft), appendLF = FALSE)
  
  parts <- Filter(Negate(is.null), lapply(all_results, `[[`, ft))
  
  if (length(parts) == 0L) {
    message("NO DATA")
    inventory[file_type == ft, status := "missing"]
    next
  }
  
  combined <- rbindlist(parts, fill = TRUE, use.names = TRUE)
  rm(parts); gc()
  
  # Enforce integer types on key columns
  for (col in c("case_id","cell_id","sim_year","cal_year","Year","Day")) {
    if (col %in% names(combined))
      set(combined, j = col, value = as.integer(combined[[col]]))
  }
  
  # Sort
  sk <- SORT_KEYS[[ft]][SORT_KEYS[[ft]] %in% names(combined)]
  if (length(sk)) setkeyv(combined, sk)
  
  # Save
  fname <- if (ft == "Annual") "annual_allcases.rds" else
    sprintf("daily_%s_allcases.rds", ft)
  fpath <- file.path(OUTPUT_DIR, fname)
  saveRDS(combined, fpath)
  if (ft == "Annual") fwrite(combined, sub("\\.rds$", ".csv", fpath))
  
  sz <- round(file.info(fpath)$size / 1024^2, 1)
  message(sprintf("%s rows x %d cols | %.1f MB -> %s",
                  format(nrow(combined), big.mark=","), ncol(combined), sz, fname))
  
  inventory[file_type == ft,
            `:=`(rows    = nrow(combined),
                 cols    = ncol(combined),
                 size_mb = sz,
                 status  = "saved")]
  rm(combined); gc()
}


# =============================================================================
# 8. BUILD daily_panel_allcases.rds  (scene-by-scene merge to control RAM)
# =============================================================================
# Strategy: load all 4 daily component RDS files, then loop over the 13 scenes.
# Each scene chunk is ~263 cells x 20 years x 365 days = ~1.9M rows x ~57 cols
# = ~0.9 GB, well within a 16 GB machine.  Chunks are appended sequentially.
#
# Key columns on every row: case_id, cell_id, scene_id, sim_year, cal_year, Day
# Data columns: all SoilC + SoilClimate + SoilWater + FieldCrop columns
# =============================================================================

message("\nBuilding daily_panel_allcases.rds (scene-by-scene merge) ...")

PANEL_PATH   <- file.path(OUTPUT_DIR, "daily_panel_allcases.rds")
JOIN_KEYS    <- c("case_id", "cell_id", "scene_id", "sim_year", "cal_year", "Day")
DAILY_FTYPES <- c("SoilC", "SoilClimate", "SoilWater", "FieldCrop")

# Check all component files exist
comp_files <- setNames(
  file.path(OUTPUT_DIR, sprintf("daily_%s_allcases.rds", DAILY_FTYPES)),
  DAILY_FTYPES
)
missing_comp <- comp_files[!file.exists(comp_files)]
if (length(missing_comp) > 0)
  stop("Missing component files needed for panel merge:\n",
       paste(" ", missing_comp, collapse = "\n"),
       "\nRun the extraction section above first.")

# Load all 4 components into named list; set keys for fast binary join
message("  Loading components ...")
comp <- lapply(names(comp_files), function(ft) {
  message(sprintf("    %s ...", ft), appendLF = FALSE)
  dt <- readRDS(comp_files[[ft]])
  # Ensure key columns are integer
  for (col in c("case_id","cell_id","sim_year","cal_year","Day"))
    if (col %in% names(dt)) set(dt, j=col, value=as.integer(dt[[col]]))
  setkeyv(dt, JOIN_KEYS[JOIN_KEYS %in% names(dt)])
  message(sprintf(" %s rows", format(nrow(dt), big.mark=",")))
  dt
})
names(comp) <- names(comp_files)

# All unique scene_ids (from the largest component)
ref_dt    <- comp[[which.max(sapply(comp, nrow))]]
all_scenes <- sort(unique(ref_dt$scene_id))
message(sprintf("  Scenes to merge: %s", paste(all_scenes, collapse=" ")))

# Scene-by-scene merge — each chunk is safe in RAM
panel_chunks <- vector("list", length(all_scenes))

for (s_i in seq_along(all_scenes)) {
  
  sid <- all_scenes[s_i]
  message(sprintf("  Scene %s (%d/%d) ...", sid, s_i, length(all_scenes)),
          appendLF = FALSE)
  
  # Slice this scene from every component
  scene_parts <- lapply(comp, function(dt) dt[scene_id == sid])
  scene_parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, scene_parts)
  
  if (length(scene_parts) == 0) {
    message(" no data, skipping")
    next
  }
  
  # Merge the 4 components on their common key columns
  scene_panel <- Reduce(
    function(a, b) {
      by_cols <- intersect(JOIN_KEYS, intersect(names(a), names(b)))
      merge(a, b, by = by_cols, all = TRUE, sort = FALSE)
    },
    scene_parts[-1],
    init = scene_parts[[1]]
  )
  
  # Sort within scene
  sk <- JOIN_KEYS[JOIN_KEYS %in% names(scene_panel)]
  setkeyv(scene_panel, sk)
  
  panel_chunks[[s_i]] <- scene_panel
  message(sprintf(" %s rows x %d cols", format(nrow(scene_panel), big.mark=","),
                  ncol(scene_panel)))
  
  rm(scene_parts, scene_panel); gc()
}

# Free component RAM before final rbind
rm(comp); gc()

# Combine all scene chunks
message("  Combining all scenes ...")
panel <- rbindlist(Filter(Negate(is.null), panel_chunks),
                   fill = TRUE, use.names = TRUE)
rm(panel_chunks); gc()

# Final sort: scene_id -> case_id -> cell_id -> sim_year -> Day
final_sk <- c("scene_id","case_id","cell_id","sim_year","cal_year","Day")
setkeyv(panel, final_sk[final_sk %in% names(panel)])

message(sprintf("  Panel: %s rows x %d cols",
                format(nrow(panel), big.mark=","), ncol(panel)))
message(sprintf("  Identifiers present: case_id=%s | cell_id=%s | scene_id=%s",
                "case_id" %in% names(panel),
                "cell_id" %in% names(panel),
                "scene_id" %in% names(panel)))

saveRDS(panel, PANEL_PATH)
panel_mb <- round(file.info(PANEL_PATH)$size / 1024^2, 1)
message(sprintf("  Saved: daily_panel_allcases.rds | %.1f MB", panel_mb))

# Quick CH4 sanity check on the panel
if ("CH4-flux" %in% names(panel)) {
  v <- panel[["CH4-flux"]]
  v <- v[!is.na(v) & v > 0]
  message(sprintf("  CH4-flux (>0): n=%s | mean=%.4f | p95=%.4f kg C/ha/d",
                  format(length(v), big.mark=","), mean(v), quantile(v, 0.95)))
}

rm(panel); gc()

# =============================================================================
# 9. VALIDATION REPORT
# =============================================================================

message("\n", strrep("=", 68))
message("EXTRACTION SUMMARY")
message(strrep("=", 68))
message(sprintf("  Cases matched  : %d", nrow(case_index)))
message(sprintf("  Unique cells   : %d", uniqueN(case_index$cell_id)))
message(sprintf("  Unique scenes  : %s",
                paste(sort(unique(case_index$scene_id)), collapse=" ")))
message(sprintf("  Elapsed        : %.1f min", elapsed / 60))
message("")
print(inventory)

# CH4 spot-check
sc_path <- file.path(OUTPUT_DIR, "daily_SoilC_allcases.rds")
if (file.exists(sc_path)) {
  sc <- readRDS(sc_path)
  message("\nCH4 check (daily_SoilC_allcases):")
  message(sprintf("  Rows / cases / cells / scenes : %s / %d / %d / %d",
                  format(nrow(sc),big.mark=","),
                  uniqueN(sc$case_id), uniqueN(sc$cell_id), uniqueN(sc$scene_id)))
  if ("CH4-flux" %in% names(sc)) {
    v <- sc[["CH4-flux"]]
    v <- v[!is.na(v) & v > 0]
    message(sprintf("  CH4-flux >0: n=%s | mean=%.4f | p95=%.4f kg C/ha/d",
                    format(length(v),big.mark=","), mean(v), quantile(v,0.95)))
    # Per-scene annual sum normalised by cell-years
    sc[, cy_key := paste(cell_id, cal_year)]
    ch4_s <- sc[, .(
      n_cell_yrs = uniqueN(cy_key),
      ch4_sum_per_cell_yr = round(sum(`CH4-flux`, na.rm=TRUE) / uniqueN(cy_key), 2)
    ), by = scene_id][order(scene_id)]
    message("\n  Annual CH4 total (kg C/ha) per cell-year by scene:")
    print(ch4_s)
  }
  rm(sc); gc()
}

# Annual check
an_path <- file.path(OUTPUT_DIR, "annual_allcases.rds")
if (file.exists(an_path)) {
  an <- readRDS(an_path)
  message("\nAnnual summary check:")
  message(sprintf("  Rows / scenes : %s / %d",
                  format(nrow(an),big.mark=","), uniqueN(an$scene_id)))
  if (all(c("CH4","GrainC1") %in% names(an))) {
    s <- an[, .(CH4_mean   = round(mean(CH4,     na.rm=TRUE), 2),
                Grain1_mean= round(mean(GrainC1,  na.rm=TRUE), 1)),
            by = scene_id][order(scene_id)]
    message("  Annual CH4 (kg C/ha/yr) and GrainC1 (kg C/ha/yr) by scene:")
    print(s)
  }
  rm(an); gc()
}

message("\nScript 5 complete. Pass these files to Script 6:")
fns <- c("case_scene_index.csv", "annual_allcases.rds",
         "daily_SoilC_allcases.rds", "daily_SoilClimate_allcases.rds",
         "daily_SoilWater_allcases.rds", "daily_FieldCrop_allcases.rds")
for (fn in fns) message(sprintf("  %s/%s", OUTPUT_DIR, fn))