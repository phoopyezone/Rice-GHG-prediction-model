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