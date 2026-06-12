##### 4_dndc_make_batch_MULTI.R #####
# Creates DNDC batch file listing all multi-scenario cell .dnd input files

# ---- Paths ----
input_dir  <- "C:\\DNDC\\simulate_dndc\\dndc_inputs\\input_files"
batch_file <- file.path(input_dir, "myanmar_rice_batch_all.txt")

# ---- Find and sort all scenario .dnd files ----
# Expected pattern from script 3:
# Cell_<cellid>_S<scenario>.dnd
dnd_files <- list.files(
  input_dir,
  pattern = "^Cell_[0-9]+_S[0-9]{3}\\.dnd$",
  full.names = TRUE
)

if (length(dnd_files) == 0) stop("No multi-scenario .dnd files found in: ", input_dir)

# ---- Extract cell ID and scenario code for sorting ----
bn <- basename(dnd_files)
cell_nums <- as.integer(sub("^Cell_([0-9]+)_S[0-9]{3}\\.dnd$", "\\1", bn))
scen_nums <- as.integer(sub("^Cell_[0-9]+_S([0-9]{3})\\.dnd$", "\\1", bn))

ord <- order(scen_nums, cell_nums)
dnd_files <- dnd_files[ord]
cell_nums <- cell_nums[ord]
scen_nums <- scen_nums[ord]

# ---- Write batch file ----
batch_lines <- c(
  as.character(length(dnd_files)),
  normalizePath(dnd_files, winslash = "\\", mustWork = FALSE)
)

writeLines(batch_lines, batch_file)

# ---- Summary ----
cat(sprintf("=== DNDC Batch File Created ===\n"))
cat(sprintf("Output    : %s\n", batch_file))
cat(sprintf("Files     : %d\n", length(dnd_files)))
cat(sprintf("Scenarios : %d  (S%03d to S%03d)\n",
            length(unique(scen_nums)), min(scen_nums), max(scen_nums)))
cat(sprintf("Cells     : %d  (Cell_%d to Cell_%d)\n",
            length(unique(cell_nums)), min(cell_nums), max(cell_nums)))

cat(sprintf("\nPreview (first 6 lines):\n"))
cat(paste(readLines(batch_file, n = 6), collapse = "\n"), "\n")

cat(sprintf("\nTo run: Open DNDC -> Tools -> Run Batch -> select %s\n",
            basename(batch_file)))

##### END OF: 4_dndc_make_batch_MULTI.R #####