##### SIMPLE_make_batch.R #####
# Minimal script to create DNDC batch file

setwd("G:/My Drive/Research/simulation/main_dndc/simulate_dndc")

# Find all .dnd files
dnd_files <- list.files("dndc_inputs/input files", 
                        pattern = "\\.dnd$", 
                        full.names = TRUE)

# Sort by cell number
cell_nums <- as.integer(gsub(".*Cell_(\\d+)\\.dnd", "\\1", dnd_files))
dnd_files <- dnd_files[order(cell_nums)]

# Create batch file (Windows format)
batch_content <- c(
  length(dnd_files),  # First line: count
  normalizePath(dnd_files, winslash = "\\")  # Full paths with backslashes
)

# Write
batch_file <- "dndc_inputs/input files/myanmar_rice_batch.txt"
writeLines(as.character(batch_content), batch_file)

# Report
cat(sprintf("\n✓ Created: %s\n", batch_file))
cat(sprintf("  Files: %d\n", length(dnd_files)))
cat(sprintf("  Range: Cell_%d to Cell_%d\n\n", min(cell_nums), max(cell_nums)))

# Preview
cat("Preview (first 3 lines):\n")
cat(paste(readLines(batch_file, n = 3), collapse = "\n"), "\n")

##### END #####