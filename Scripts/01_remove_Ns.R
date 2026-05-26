# Remove N's for DADA2 functionality using filter and trim

# get wd & data location
path <- getwd()

# Use dataset-specific output directory
output_dir <- file.path(path, "Data", "Temp", test_data_name)

# Define subdirectories
filtN_dir <- file.path(output_dir, "01_filtN")
rds_dir   <- file.path(output_dir, "R_objects")

# Create directories if they don't exist
dir.create(filtN_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# group raw files into "forward" reads and "reverse" reads
fnFs <- sort(manifest$absolute_forward_path)
fnRs <- sort(manifest$absolute_backward_path)

# set up filepaths for output files
fnFs.filtN <- file.path(filtN_dir, basename(fnFs))
fnRs.filtN <- file.path(filtN_dir, basename(fnRs))

# run filterAndTrim for all files
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0)

# write objects to pass to next script
saveRDS(fnFs.filtN, file = file.path(rds_dir, "01_fnFs.filtN.rds"))
saveRDS(fnRs.filtN, file = file.path(rds_dir, "01_fnRs.filtN.rds"))