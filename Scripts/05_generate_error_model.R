# Script to generate error models

# get wd
path <- getwd()

# Dataset-specific directories
output_dir  <- file.path(path, "Data", "Temp", test_data_name)
rds_dir     <- file.path(output_dir, "R_objects")
results_dir <- file.path(path, "Results", test_data_name)

# ensure directories exist
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# read in the lists of filtered reads
filtFs <- readRDS(file = file.path(rds_dir, "04_fnFs.filtN.rds"))
filtRs <- readRDS(file = file.path(rds_dir, "04_fnRs.filtN.rds"))
out    <- readRDS(file = file.path(rds_dir, "04_out.rds"))

# check that files exist (some may have been dropped)
exists <- file.exists(filtFs) & file.exists(filtRs)

if (!any(exists)) {
  stop("No valid filtered files found — check previous step.")
}

# run error learning only on valid files
errF <- learnErrors(filtFs[exists], multithread = TRUE, nbases = 2e8)
errR <- learnErrors(filtRs[exists], multithread = TRUE, nbases = 2e8)

# save error models
saveRDS(errF, file = file.path(rds_dir, "05_errF.rds"))
saveRDS(errR, file = file.path(rds_dir, "05_errR.rds"))

# write plot to file
pdf(
  file = file.path(results_dir, "05_error_rate_plots.pdf"),
  width = 10,
  height = 10
)