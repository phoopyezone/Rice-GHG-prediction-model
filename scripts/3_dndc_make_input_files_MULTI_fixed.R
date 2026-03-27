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
