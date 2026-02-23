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
  tmax <- rast(file.path(power_dir, "tmax-2005_2024-91.5x101.5x8x29.nc"))
  tmin <- rast(file.path(power_dir, "tmin-2005_2024-91.5x101.5x8x29.nc"))
  prec <- rast(file.path(power_dir, "prec-2005_2024-91.5x101.5x8x29.nc"))
  
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