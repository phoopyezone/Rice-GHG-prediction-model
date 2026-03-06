##### 5_dndc_extract_outputs.R #####
# Extracts daily DNDC outputs for all 207 Myanmar rice cells
# Combines 5 daily CSVs per year per cell into one long-format panel dataset
# Outpu F: CH4-flux, yield (GrainC), NEE, N2O-flux
# Input features:  soil climate, water, crop state, management
# Output: one RDS per file type + one merged daily panel RDS
#
# File sources per cell-year:
#   Day_SoilC_{yr}.csv       -> CH4 pathway + carbon stocks/fluxes
#   Day_SoilClimate_{yr}.csv -> soil temp, moisture, Eh, water table
#   Day_SoilWater_{yr}.csv   -> ponding, precipitation, ET, runoff
#   Day_FieldCrop_{yr}.csv   -> crop biomass, LAI, GrainC (yield proxy)
#   Day_FieldManage_{yr}.csv -> fertiliser, irrigation, management events
#   Multi_year_summary.csv   -> annual grain yield + annual CH4 (for labels)
##### 5_dndc_extract_outputs.R #####
# Extracts daily DNDC outputs for all 207 Myanmar rice cells
# Reads 5 daily file types x 20 years per cell, merges into one panel
###############################################################################


library(data.table)
library(stringr)

# ── Paths ──────────────────────────────────────────────────────────────────
batch_dir  <- "C:/DNDC/Result/Record/Batch"
output_dir <- "C:/DNDC/simulate_dndc/dndc_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── SoilClimate: target column names after position-based extraction ────────
# Row structure: row0=title, row1=group labels, row2=depth header, row3=blank
# fread deduplicates repeated depth cols:
#   Temp    = 1, 5, 10, 20
#   Moisture= 1.1, 5.1, 10.1, 20.1
#   O2      = 1.2, 10.2, 20.2  (skipped - not needed)
#   Eh      = 1.3, 10.3, 20.3
#   Ice_wfps= 1.4, 10.4, 20.4  (skipped)
#   SoilWater/DeepWater = (mm).2, (mm).3
#   Soil_pH = 1.5, 10.5, 20.5
soilclimate_src  <- c("Day","1","5","10","20",
                      "1.1","5.1","10.1","20.1",
                      "1.3","10.3","20.3",
                      "(mm).2","(mm).3","1.5")
soilclimate_dest <- c("Day",
                      "SoilT_1cm","SoilT_5cm","SoilT_10cm","SoilT_20cm",
                      "SoilM_1cm","SoilM_5cm","SoilM_10cm","SoilM_20cm",
                      "SoilEh_1cm","SoilEh_10cm","SoilEh_20cm",
                      "SoilWater_mm","DeepWater_mm","Soil_pH_1cm")

# ── Columns to keep per file type ──────────────────────────────────────────
# Names come directly from the actual file headers (verified from uploads)
keep_cols <- list(
  # Day_SoilC: skip=1, header on row1
  # Actual header: Day,Very labile litter,Labile litter,...,CH4-prod.,CH4-oxid.,CH4-flux,...
  SoilC = c("Day",
            "SOC","dSOC","DOC","Humads","Humus","Microbe",
            "CH4-prod.","CH4-oxid.","CH4-flux","CH4-pool",
            "NPP","NEE","Soil-heterotrophic-respiration","Photosynthesis",
            "Litter-C","Manure-C","DOC-leach"),
  
  # Day_SoilClimate: position-extracted, renamed via soilclimate_dest
  SoilClimate = soilclimate_dest,
  
  # Day_SoilWater: skip=1, header on row1
  # Actual header: Day,IniSoilWater,IniDeepWater,IniWaterPool,...,Ponding,...
  SoilWater = c("Day",
                "Ponding","Precipitation","Irrigation",
                "Evaporation","Transpiration","Leaching","Runoff",
                "IniSoilWater","EndSoilWater"),
  
  # Day_FieldCrop: rows 1+2 merged (row2 fills blanks in row1)
  # Merged header has: Day,Temperature,...,LeafC,StemC,RootC,GrainC,CropID,TDD,GrowthIndex,...
  FieldCrop = c("Day",
                "LeafC","StemC","RootC","GrainC",
                "TDD","GrowthIndex",
                "Water_stress","N_stress","LAI",
                "TotalCropN","DailyCropGrowth","DayGrainGrowth"),
  
  # Day_FieldManage: row0=title, row1=blank, row2=header, row3=blank
  # Actual header: Day,Irrigation (mm),Fertilizer (kgN/ha),Fertilizer (kgP/ha),...
  FieldManage = c("Day","Irrigation (mm)","Fertilizer (kgN/ha)","Manure (kgN/ha)")
)

# ── Multi_year_summary column names (from DNDC manual + actual file) ────────
# Manual-defined order (52 columns). If file has more, extras kept as-is.
# NOTE: if your file has 72 columns, run names(fread(path, nrows=0)) once
# and paste below to fill in positions 53-72.
annual_col_names <- c(
  "Year",
  # Crop biomass (3 crops x 3 components)
  "GrainC1","LeafStemC1","RootC1",
  "GrainC2","LeafStemC2","RootC2",
  "GrainC3","LeafStemC3","RootC3",
  # SOC at 3 depths + annual dynamics
  "SOC_0_10cm","SOC_0_20cm","SOC_0_30cm",
  "Ini_SOC","End_SOC","dSOC",
  # Carbon inputs/outputs
  "LitterC_input","RootC_input","ManureC_input",
  "Soil_CO2","CH4",
  # Nitrogen pools
  "Ini_SON","Ini_SIN","End_SON","End_SIN","dSN",
  # Nitrogen inputs
  "Atmo_N_input","Fertilizer_N_input","Manure_N_input","Litter_N_input","N_fixation",
  # Nitrogen outputs
  "Crop_N_uptake","N_leach","N_runoff",
  "N2O_flux","NO_flux","N2_flux","NH3_flux",
  "Exch_NH4",
  # Water budget
  "PET","Transpiration","Evaporation","WaterLeach","Runoff",
  "Irrigation","Precipitation",
  # Climate + stress
  "MeanT","WindSpeed",
  "ColdStress","WaterStress","N_Stress",
  "Cut_CropC"   # col 52
  # Positions 53-72: unknown - upload Multi_year_summary.csv to identify
)

# ── Custom readers ──────────────────────────────────────────────────────────
read_soilC <- function(path, cell_id, yr) {
  dt <- tryCatch(fread(path, skip=1, header=TRUE, fill=TRUE, showProgress=FALSE),
                 error=function(e) NULL)
  if (is.null(dt) || nrow(dt)==0) return(NULL)
  dt <- dt[suppressWarnings(!is.na(as.numeric(as.character(dt[[1]]))))]
  if (nrow(dt)==0) return(NULL)
  dt[, intersect(keep_cols$SoilC, names(dt)), with=FALSE][, `:=`(cell_id=cell_id, year=yr)]
}

read_soilclimate <- function(path, cell_id, yr) {
  lns <- tryCatch(readLines(path, warn=FALSE), error=function(e) NULL)
  if (is.null(lns)) return(NULL)
  hdr_idx <- which(grepl("^Day,", lns))[1]
  if (is.na(hdr_idx)) return(NULL)
  data_lns <- lns[(hdr_idx+1):length(lns)]
  data_lns <- data_lns[nchar(trimws(data_lns)) > 0]
  if (length(data_lns)==0) return(NULL)
  dt <- tryCatch(fread(text=paste(c(lns[hdr_idx], data_lns), collapse="\n"),
                       header=TRUE, fill=TRUE, showProgress=FALSE), error=function(e) NULL)
  if (is.null(dt) || nrow(dt)==0) return(NULL)
  dt <- dt[suppressWarnings(!is.na(as.numeric(as.character(dt[[1]]))))]
  if (nrow(dt)==0) return(NULL)
  present <- soilclimate_src %in% names(dt)
  if (!any(present)) return(NULL)
  out <- tryCatch(dt[, soilclimate_src[present], with=FALSE], error=function(e) NULL)
  if (is.null(out)) return(NULL)
  setnames(out, soilclimate_dest[present])
  out[, `:=`(cell_id=cell_id, year=yr)]
}

read_soilwater <- function(path, cell_id, yr) {
  dt <- tryCatch(fread(path, skip=1, header=TRUE, fill=TRUE, showProgress=FALSE),
                 error=function(e) NULL)
  if (is.null(dt) || nrow(dt)==0) return(NULL)
  dt <- dt[suppressWarnings(!is.na(as.numeric(as.character(dt[[1]]))))]
  if (nrow(dt)==0) return(NULL)
  dt[, intersect(keep_cols$SoilWater, names(dt)), with=FALSE][, `:=`(cell_id=cell_id, year=yr)]
}

read_fieldcrop <- function(path, cell_id, yr) {
  lns <- tryCatch(readLines(path, n=5), error=function(e) NULL)
  if (is.null(lns) || length(lns) < 3) return(NULL)
  r1  <- trimws(strsplit(lns[2], ",")[[1]])
  r2  <- trimws(strsplit(lns[3], ",")[[1]])
  n   <- max(length(r1), length(r2))
  r1  <- c(r1, rep("", n - length(r1)))
  r2  <- c(r2, rep("", n - length(r2)))
  hdr <- ifelse(r2 != "", r2, r1)
  dt  <- tryCatch(fread(path, skip=4, header=FALSE, fill=TRUE, showProgress=FALSE),
                  error=function(e) NULL)
  if (is.null(dt) || nrow(dt)==0) return(NULL)
  setnames(dt, seq_len(min(ncol(dt), length(hdr))), hdr[seq_len(min(ncol(dt), length(hdr)))])
  dt <- dt[suppressWarnings(!is.na(as.numeric(as.character(dt[[1]]))))]
  if (nrow(dt)==0) return(NULL)
  dt[, intersect(keep_cols$FieldCrop, names(dt)), with=FALSE][, `:=`(cell_id=cell_id, year=yr)]
}

read_fieldmanage <- function(path, cell_id, yr) {
  lns <- tryCatch(readLines(path, warn=FALSE), error=function(e) NULL)
  if (is.null(lns)) return(NULL)
  hdr_idx <- which(grepl("^Day,", lns))[1]
  if (is.na(hdr_idx)) return(NULL)
  hdr      <- trimws(strsplit(lns[hdr_idx], ",")[[1]])
  data_lns <- lns[(hdr_idx+1):length(lns)]
  data_lns <- data_lns[nchar(trimws(data_lns)) > 0]
  if (length(data_lns)==0) return(NULL)
  dt <- tryCatch(fread(text=paste(data_lns, collapse="\n"), header=FALSE,
                       fill=TRUE, showProgress=FALSE), error=function(e) NULL)
  if (is.null(dt) || nrow(dt)==0) return(NULL)
  setnames(dt, seq_len(min(ncol(dt), length(hdr))), hdr[seq_len(min(ncol(dt), length(hdr)))])
  dt <- dt[suppressWarnings(!is.na(as.numeric(as.character(dt[[1]]))))]
  if (nrow(dt)==0) return(NULL)
  dt[, intersect(keep_cols$FieldManage, names(dt)), with=FALSE][, `:=`(cell_id=cell_id, year=yr)]
}

read_annual <- function(path, cell_id) {
  dt <- tryCatch(fread(path, header=FALSE, fill=TRUE, showProgress=FALSE),
                 error=function(e) NULL)
  if (is.null(dt) || nrow(dt)==0) return(NULL)
  # Assign names up to however many columns exist
  n_cols   <- ncol(dt)
  n_named  <- min(n_cols, length(annual_col_names))
  col_names <- c(annual_col_names[seq_len(n_named)],
                 if (n_cols > n_named) paste0("V", seq(n_named+1, n_cols)) else character(0))
  setnames(dt, col_names)
  dt[, cell_id := cell_id]
  dt
}

readers <- list(SoilC=read_soilC, SoilClimate=read_soilclimate, SoilWater=read_soilwater,
                FieldCrop=read_fieldcrop, FieldManage=read_fieldmanage)

# ── Discover cell folders ───────────────────────────────────────────────────
cell_folders <- list.dirs(batch_dir, full.names=TRUE, recursive=FALSE)
cell_folders <- cell_folders[grepl("Myanmar_Rice_Cell_", basename(cell_folders))]
cat(sprintf("Found %d cell folders\n", length(cell_folders)))

# ── Main loop ───────────────────────────────────────────────────────────────
results <- list(SoilC=list(), SoilClimate=list(), SoilWater=list(),
                FieldCrop=list(), FieldManage=list(), Annual=list())

for (i in seq_along(cell_folders)) {
  folder  <- cell_folders[i]
  cell_id <- as.integer(str_extract(basename(folder), "(?<=Cell_)\\d+"))
  if (i %% 25 == 1) cat(sprintf("  [%d/%d] cell %d\n", i, length(cell_folders), cell_id))
  
  results$Annual[[i]] <- read_annual(file.path(folder, "Multi_year_summary.csv"), cell_id)
  
  for (ftype in names(readers)) {
    yr_list <- lapply(1:20, function(yr)
      tryCatch(readers[[ftype]](file.path(folder, sprintf("Day_%s_%d.csv", ftype, yr)), cell_id, yr),
               error=function(e) NULL))
    results[[ftype]][[i]] <- rbindlist(yr_list, fill=TRUE, use.names=TRUE)
  }
}

# ── Combine + clean ─────────────────────────────────────────────────────────
cat("Combining all cells...\n")
combined <- lapply(results, function(x) rbindlist(x, fill=TRUE, use.names=TRUE))

# Rename FieldManage cols to clean names
setnames(combined$FieldManage,
         c("Irrigation (mm)","Fertilizer (kgN/ha)","Manure (kgN/ha)"),
         c("Irrigation_mm","Fertilizer_kgN_ha","Manure_kgN_ha"), skip_absent=TRUE)

# ── Merge into daily panel ───────────────────────────────────────────────────
cat("Merging into daily panel...\n")
keys  <- c("cell_id", "year", "Day")
panel <- Reduce(function(a,b) merge(a, b, by=keys, all=TRUE), combined[names(readers)])
setorder(panel, cell_id, year, Day)
cat(sprintf("Panel: %s rows x %d cols\n", format(nrow(panel), big.mark=","), ncol(panel)))

# ── Save ─────────────────────────────────────────────────────────────────────
saveRDS(panel,           file.path(output_dir, "daily_panel_allcells.rds"))
saveRDS(combined$Annual, file.path(output_dir, "annual_summary_allcells.rds"))
fwrite(combined$Annual,  file.path(output_dir, "annual_summary_allcells.csv"))
cat("Done. Output:", output_dir, "\n")