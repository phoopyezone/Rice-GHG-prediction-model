##### 3_make_dndc_input_files.R #####
# Generates DNDC input files (.dnd) - one per cell, covering 20 years
# Each .dnd references all 20 climate files (one per year)
# Creates batch list file for DNDC batch mode

suppressPackageStartupMessages(library(dplyr))

path <- "G:/My Drive/Research/simulation/main_dndc/simulate_dndc"
setwd(path)

template_file <- file.path("dndc_inputs", "input files", "template", "mm_dndc_input_template.txt")
output_dir <- file.path("dndc_inputs", "input files")
climate_rel_path <- "dndc_inputs/climate_files"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cells <- readRDS("data/cells.rds")
years <- 2005:2024
n_years <- length(years)

message(sprintf("\n=== DNDC Input File Generator ==="))
message(sprintf("Cells: %d | Years: %d", nrow(cells), n_years))
message(sprintf("Output: %d .dnd files + 1 batch file", nrow(cells)))

# Check column names and print for reference
message("\nAvailable soil columns:")
soil_cols <- grep("bdod|clay|sand|silt|soc|pH|poros", names(cells), value = TRUE, ignore.case = TRUE)
message(paste(soil_cols, collapse = ", "))

# ---- Read template ----
template <- readLines(template_file)

# ---- Helper: Safe column access ----
get_col <- function(df, patterns) {
  for (pat in patterns) {
    matches <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) return(df[[matches[1]]])
  }
  return(NA)
}

# ---- Main generation function ----
generate_dnd_file <- function(cell_row, years, template) {
  cell_id <- cell_row$cell
  lines <- template
  
  # 1) UPDATE SITE INFO
  lines <- gsub("^(__Site_name\\s+).*",
                sprintf("\\1Myanmar_Rice_Cell_%d", cell_id), lines)
  
  lines <- gsub("^(__Simulated_years\\s+).*",
                sprintf("\\1%d", n_years), lines)
  
  lines <- gsub("^(__Latitude\\s+).*",
                sprintf("\\1%.4f", cell_row$y), lines)
  
  # 2) UPDATE CLIMATE FILES SECTION
  climate_count_idx <- grep("^__Climate_files\\s+", lines)
  old_file_idx <- grep("^[0-9]+\\s+.*climate", lines, ignore.case = TRUE)
  
  lines[climate_count_idx] <- sprintf("__Climate_files%46s%d", "", n_years)
  lines <- lines[-old_file_idx]
  
  climate_lines <- sapply(1:n_years, function(i) {
    sprintf("%-2d                                       %s/Cell_%d_%d",
            i, climate_rel_path, cell_id, years[i])
  })
  lines <- append(lines, climate_lines, after = climate_count_idx)
  
  # 3) UPDATE CROPPING SYSTEM
  lines <- gsub("^(__Total_years\\s+).*",
                sprintf("\\1%d", n_years), lines)
  
  # 4) UPDATE SOIL PARAMETERS - Use flexible column matching
  # Bulk density
  bd <- get_col(cell_row, c("bdod_0-5cm", "bdod_0.5cm", "bdod"))
  if (!is.na(bd)) {
    lines <- gsub("^(__Bulk_density\\s+).*",
                  sprintf("\\1%.4f", bd), lines)
  }
  
  # pH
  ph <- get_col(cell_row, c("pH_0-5cm", "pH_0.5cm", "phh2o_0-5cm", "pH"))
  if (!is.na(ph)) {
    lines <- gsub("^(__pH\\s+).*",
                  sprintf("\\1%.4f", ph), lines)
  }
  
  # Clay
  clay <- get_col(cell_row, c("clay_0-5cm", "clay_0.5cm", "clay"))
  if (!is.na(clay)) {
    lines <- gsub("^(__Clay_fraction\\s+).*",
                  sprintf("\\1%.4f", clay), lines)
  }
  
  # Porosity
  poros <- get_col(cell_row, c("poros_0-5cm", "poros_0.5cm", "porosity"))
  if (!is.na(poros)) {
    lines <- gsub("^(__Porosity\\s+).*",
                  sprintf("\\1%.4f", poros), lines)
  }
  
  # SOC
  soc <- get_col(cell_row, c("soc_0-5cm", "soc_0.5cm", "soc"))
  if (!is.na(soc)) {
    lines <- gsub("^(__Top_layer_SOC\\s+).*",
                  sprintf("\\1%.4f", soc), lines)
  }
  
  # Sand (for calculations)
  sand <- get_col(cell_row, c("sand_0-5cm", "sand_0.5cm", "sand"))
  
  # Calculate field capacity and wilting point
  if (!is.na(clay) && !is.na(sand)) {
    clay_pct <- clay * 100
    sand_pct <- sand * 100
    soc_pct <- ifelse(is.na(soc), 2, soc * 100)
    
    fc <- 0.2576 - 0.002 * sand_pct + 0.0036 * clay_pct + 0.0299 * soc_pct
    wp <- 0.026 + 0.005 * clay_pct + 0.0158 * soc_pct
    
    fc <- max(0.15, min(0.55, fc))
    wp <- max(0.05, min(0.25, wp))
    
    lines <- gsub("^(__Field_capacity\\s+).*",
                  sprintf("\\1%.4f", fc), lines)
    lines <- gsub("^(__Wilting_point\\s+).*",
                  sprintf("\\1%.4f", wp), lines)
  }
  
  # Soil texture ID from clay content
  if (!is.na(clay)) {
    clay_pct <- clay * 100
    if (clay_pct < 7) {
      texture_id <- 2
    } else if (clay_pct < 15) {
      texture_id <- 3
    } else if (clay_pct < 27) {
      texture_id <- 4
    } else if (clay_pct < 40) {
      texture_id <- 8
    } else {
      texture_id <- 12
    }
    lines <- gsub("^(__Soil_texture_ID\\s+).*",
                  sprintf("\\1%d", texture_id), lines)
  }
  
  lines
}

# ---- Generate all .dnd files ----
message("\nGenerating input files...")
dnd_files <- character(nrow(cells))

for (i in 1:nrow(cells)) {
  lines <- generate_dnd_file(cells[i, ], years, template)
  
  output_file <- file.path(output_dir, sprintf("Cell_%d.dnd", cells$cell[i]))
  writeLines(lines, output_file)
  dnd_files[i] <- output_file
  
  if (i %% 10 == 0 || i == nrow(cells)) {
    message(sprintf("  Progress: %d/%d", i, nrow(cells)))
  }
}

# ---- Create batch list file ----
message("\nCreating batch list file...")
batch_file <- file.path(output_dir, "myanmar_rice_batch.txt")

batch_lines <- c(
  as.character(nrow(cells)),
  normalizePath(dnd_files, winslash = "\\", mustWork = FALSE)
)

writeLines(batch_lines, batch_file)

# ---- Summary ----
message("\n" <- paste(rep("=", 60), collapse = ""))
message("=== SUMMARY ===")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("Input files (.dnd): %d", nrow(cells)))
message(sprintf("Years per file: %d (2005-2024)", n_years))
message(sprintf("Total simulations: %d cell-years", nrow(cells) * n_years))
message(sprintf("\nOutput: %s", output_dir))
message(sprintf("Batch file: %s", basename(batch_file)))

# ---- Validation ----
message("\n=== VALIDATION ===")

sample_file <- dnd_files[1]
sample_lines <- readLines(sample_file)

site_line <- grep("^__Site_name", sample_lines, value = TRUE)[1]
years_line <- grep("^__Simulated_years", sample_lines, value = TRUE)[1]
climate_count_line <- grep("^__Climate_files", sample_lines, value = TRUE)[1]

message(sprintf("\nSample: %s", basename(sample_file)))
message(sprintf("  %s", site_line))
message(sprintf("  %s", years_line))
message(sprintf("  %s", climate_count_line))

climate_refs <- grep("^[0-9]+\\s+.*Cell_", sample_lines, value = TRUE)
message(sprintf("\nClimate references: %d", length(climate_refs)))
if (length(climate_refs) >= 2) {
  message(sprintf("  First: %s", climate_refs[1]))
  message(sprintf("  Last: %s", tail(climate_refs, 1)))
}

soil_lines <- grep("^__(Bulk_density|pH|Clay_fraction|Porosity|Top_layer_SOC)", 
                   sample_lines, value = TRUE)
message("\nSoil parameters updated:")
for (line in head(soil_lines, 5)) message(sprintf("  %s", line))

file_sizes <- file.info(dnd_files)$size / 1024
message(sprintf("\nFile sizes: %.1f - %.1f KB", min(file_sizes), max(file_sizes)))

batch_preview <- readLines(batch_file, n = 3)
message("\nBatch file:")
message(sprintf("  Count: %s", batch_preview[1]))
message(sprintf("  First: %s", basename(batch_preview[2])))

message("\nb READY FOR DNDC BATCH SIMULATION")
message(sprintf("\nTo run: Load %s in DNDC batch mode", basename(batch_file)))

##### END OF: 3_make_dndc_input_files.R #####
