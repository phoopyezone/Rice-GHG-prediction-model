# Compile all R scripts in working_main into a single R file

folder_path <- "G:/My Drive/Research/simulation/main_dndc/scripts/working_main"

r_files <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)

if (length(r_files) == 0) {
  stop("No .R files found in: ", folder_path)
}

cat("Found", length(r_files), "R script(s):\n")
cat(paste0("  ", basename(r_files), "\n"), sep = "")

lines <- character(0)
for (f in r_files) {
  lines <- c(
    lines,
    paste0("# ---- ", basename(f), " ----"),
    readLines(f, warn = FALSE),
    ""
  )
}

output_file <- file.path(folder_path, "compiled_scripts.R")
writeLines(lines, output_file)

cat("\nCompiled output saved to:\n ", output_file, "\n")
