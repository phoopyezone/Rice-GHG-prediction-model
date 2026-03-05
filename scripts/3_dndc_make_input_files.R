##### 3_dndc_make_input_files.R #####
# Generates DNDC input files (.dnd) - one per cell, covering 20 years
# Template: C:\DNDC\simulate_dndc\dndc_inputs\input files\template\Cell_93.dnd
# Output:   C:\DNDC\simulate_dndc\dndc_inputs\input files\
# Climate section: updated per cell with C:\DNDC\... absolute paths
# Soil section: pulled from cells.rds by cell ID
# Crop / management sections: unchanged from template

suppressPackageStartupMessages(library(dplyr))

# ---- Paths ----
path          <- "C:\\DNDC\\simulate_dndc"
setwd(path)

template_file <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\input_files\\template\\Cell_93.dnd"
output_dir    <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\input_files"
climate_base  <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\climate_files"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cells   <- readRDS("C:\\DNDC\\simulate_dndc\\data\\cells.rds")
years   <- 2005:2024
n_years <- length(years)
template <- readLines(template_file)

message(sprintf("Cells: %d | Years: %d -> %d .dnd files", nrow(cells), n_years, nrow(cells)))

# ---- Soil mapping from cells.rds columns ----
# Known column names from script 1: bdod_0-5cm, clay_0-5cm, sand_0-5cm,
# silt_0-5cm, soc_0-5cm, pH_0-5cm, poros_0-5cm, texture_0-5cm

get_soil <- function(row, col) {
  val <- row[[col]]
  if (is.null(val) || length(val) == 0 || is.na(val)) NA_real_ else as.numeric(val)
}

# ---- Field capacity / wilting point pedotransfer (Saxton & Rawls) ----
calc_fc_wp <- function(clay_frac, sand_frac, soc_frac) {
  clay_pct <- clay_frac * 100; sand_pct <- sand_frac * 100
  soc_pct  <- ifelse(is.na(soc_frac), 2, soc_frac * 100)
  fc <- 0.2576 - 0.002 * sand_pct + 0.0036 * clay_pct + 0.0299 * soc_pct
  wp <- 0.026  + 0.005 * clay_pct + 0.0158 * soc_pct
  list(fc = max(0.15, min(0.55, fc)),
       wp = max(0.05, min(0.25, wp)))
}

# ---- Line replacement helper (fixed-width to match template style) ----
set_line <- function(lines, pattern, value_str) {
  idx <- grep(pattern, lines, perl = TRUE)
  if (length(idx) == 0) return(lines)
  key  <- sub("^(_{2,6}\\S+).*", "\\1", lines[idx[1]])
  lines[idx[1]] <- formatC(key, width = -51, flag = "-") |>
    paste0(value_str)
  lines
}

# ---- Main generator ----
generate_dnd <- function(cell_row, years, template) {
  id    <- cell_row$cell
  lines <- template
  
  # 1. SITE INFO
  lines <- set_line(lines, "^__Site_name\\s",
                    sprintf("Myanmar_Rice_Cell_%d", id))
  lines <- set_line(lines, "^__Simulated_years\\s",
                    sprintf("%d", n_years))
  lines <- set_line(lines, "^__Latitude\\s",
                    sprintf("%.4f", cell_row$y))
  
  # 2. CLIMATE FILES - replace block between __Climate_files and __Climate_file_mode
  cf_idx   <- grep("^__Climate_files\\s", lines)
  mode_idx <- grep("^__Climate_file_mode\\s", lines)
  
  # Remove old file-index lines (numeric entries between the two markers)
  old_rows <- which(grepl("^[0-9]+\\s+", lines) &
                      seq_along(lines) > cf_idx &
                      seq_along(lines) < mode_idx)
  if (length(old_rows)) lines <- lines[-old_rows]
  
  # Recalculate mode_idx after deletion, then insert new climate lines
  mode_idx2 <- grep("^__Climate_file_mode\\s", lines)
  new_cf <- sapply(seq_along(years), function(i)
    sprintf("%-2d                                       %s\\Cell_%d_%d",
            i, climate_base, id, years[i]))
  
  lines <- set_line(lines, "^__Climate_files\\s", sprintf("%d", n_years))
  lines <- append(lines, new_cf, after = grep("^__Climate_files\\s", lines))
  
  # 3. CROP SYSTEM: keep total years consistent
  lines <- set_line(lines, "^__Total_years\\s", sprintf("%d", n_years))
  
  # 4. SOIL PARAMETERS
  bd   <- get_soil(cell_row, "bdod_0-5cm")
  ph   <- get_soil(cell_row, "pH_0-5cm")
  clay <- get_soil(cell_row, "clay_0-5cm")
  sand <- get_soil(cell_row, "sand_0-5cm")
  soc  <- get_soil(cell_row, "soc_0-5cm")
  poros <- get_soil(cell_row, "poros_0-5cm")
  tex  <- get_soil(cell_row, "texture_0-5cm")
  
  if (!is.na(bd))    lines <- set_line(lines, "^__Bulk_density\\s",    sprintf("%.4f", bd))
  if (!is.na(ph))    lines <- set_line(lines, "^__pH\\s",              sprintf("%.4f", ph))
  if (!is.na(clay))  lines <- set_line(lines, "^__Clay_fraction\\s",   sprintf("%.4f", clay))
  if (!is.na(poros)) lines <- set_line(lines, "^__Porosity\\s",        sprintf("%.4f", poros))
  if (!is.na(soc))   lines <- set_line(lines, "^__Top_layer_SOC\\s",   sprintf("%.4f", soc))
  if (!is.na(tex))   lines <- set_line(lines, "^__Soil_texture_ID\\s", sprintf("%d", as.integer(tex)))
  
  if (!is.na(clay) && !is.na(sand)) {
    hw <- calc_fc_wp(clay, sand, soc)
    lines <- set_line(lines, "^__Field_capacity\\s",  sprintf("%.4f", hw$fc))
    lines <- set_line(lines, "^__Wilting_point\\s",   sprintf("%.4f", hw$wp))
  }
  
  lines
}

# ---- Generate all .dnd files ----
message("Generating input files...")
dnd_files <- character(nrow(cells))

for (i in seq_len(nrow(cells))) {
  lines <- generate_dnd(cells[i, ], years, template)
  out   <- file.path(output_dir, sprintf("Cell_%d.dnd", cells$cell[i]))
  writeLines(lines, out)
  dnd_files[i] <- out
  if (i %% 5 == 0 || i == nrow(cells))
    message(sprintf("  %d/%d", i, nrow(cells)))
}

# ---- Batch list ----
batch_file  <- file.path(output_dir, "myanmar_rice_batch.txt")
batch_lines <- c(as.character(nrow(cells)),
                 normalizePath(dnd_files, winslash = "\\", mustWork = FALSE))
writeLines(batch_lines, batch_file)

# ---- Summary / quick validation ----
s <- readLines(dnd_files[1])
message(sprintf("\n=== Done: %d files, batch -> %s", nrow(cells), basename(batch_file)))
message(grep("^__Site_name",     s, value = TRUE)[1])
message(grep("^__Climate_files", s, value = TRUE)[1])
message(grep("^[0-9]+\\s+.*Cell_", s, value = TRUE)[1])
message(grep("^[0-9]+\\s+.*Cell_", s, value = TRUE) |> tail(1))
message(grep("^__(Bulk_density|pH|Clay_fraction|Top_layer_SOC)", s, value = TRUE) |> paste(collapse = "\n"))

##### END OF: 3_dndc_make_input_files.R #####