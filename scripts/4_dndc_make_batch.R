##### 4_dndc_make_batch.R #####
# Creates DNDC batch file listing all cell .dnd input files

# ---- Paths ----
input_dir  <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\input_files"
batch_file <- file.path(input_dir, "myanmar_rice_batch.txt")

# ---- Find and sort all .dnd files ----
dnd_files <- list.files(input_dir, pattern = "\\.dnd$", full.names = TRUE)

if (length(dnd_files) == 0) stop("No .dnd files found in: ", input_dir)

# Sort numerically by cell ID
cell_nums <- as.integer(gsub(".*Cell_(\\d+)\\.dnd$", "\\1", basename(dnd_files)))
dnd_files <- dnd_files[order(cell_nums)]
cell_nums <- sort(cell_nums)

# ---- Write batch file ----
batch_lines <- c(
  as.character(length(dnd_files)),
  normalizePath(dnd_files, winslash = "\\", mustWork = FALSE)
)

writeLines(batch_lines, batch_file)

# ---- Summary ----
cat(sprintf("=== DNDC Batch File Created ===\n"))
cat(sprintf("Output : %s\n", batch_file))
cat(sprintf("Cells  : %d  (Cell_%d to Cell_%d)\n",
            length(dnd_files), min(cell_nums), max(cell_nums)))
cat(sprintf("\nPreview (first 4 lines):\n"))
cat(paste(readLines(batch_file, n = 4), collapse = "\n"), "\n")
cat(sprintf("\nTo run: Open DNDC -> Tools -> Run Batch -> select myanmar_rice_batch.txt\n"))

##### END OF: 4_dndc_make_batch.R #####