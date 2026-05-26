# Script to filter the data to remove any poor quality reads

# get wd
path <- getwd()

# Dataset-specific directories
output_dir  <- file.path(path, "Data", "Temp", test_data_name)
rds_dir     <- file.path(output_dir, "R_objects")
filt_dir    <- file.path(output_dir, "04_filterAndTrim")
results_dir <- file.path(path, "Results", test_data_name)

# ensure directories exist
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# read in the lists of cutadapt files for trimming
fnFs.cut <- readRDS(file = file.path(rds_dir, "02_fnFs.cut.rds"))
fnRs.cut <- readRDS(file = file.path(rds_dir, "02_fnRs.cut.rds"))

# set up filepaths for output files
fnFs.filtN <- file.path(filt_dir, basename(fnFs.cut))
fnRs.filtN <- file.path(filt_dir, basename(fnRs.cut))

# run filterAndTrim based on user defined parameters
out <- filterAndTrim(
  fnFs.cut, fnFs.filtN,
  fnRs.cut, fnRs.filtN,
  maxN = 0,
  maxEE = maxEE,
  truncQ = truncQ,
  minLen = minLen,
  truncLen = truncLen,
  multithread = TRUE
)

# save outputs
saveRDS(fnFs.filtN, file = file.path(rds_dir, "04_fnFs.filtN.rds"))
write.table(fnFs.filtN,
            file = file.path(filt_dir, "04_fnFs.filtN.txt"),
            row.names = FALSE, col.names = FALSE)

saveRDS(fnRs.filtN, file = file.path(rds_dir, "04_fnRs.filtN.rds"))
write.table(fnRs.filtN,
            file = file.path(filt_dir, "04_fnRs.filtN.txt"),
            row.names = FALSE, col.names = FALSE)

saveRDS(out, file = file.path(rds_dir, "04_out.rds"))

# generate quality plots on the quality-trimmed data
pdf(
  file = file.path(results_dir, "04_post_trim_quality_plots.pdf"),
  width = 10,
  height = 10
)

# create new manifest for filtered samples
filtered_data_loc <- filt_dir
manifest_filtered <- make_manifest(filtered_data_loc)

write.csv(
  manifest_filtered,
  file = file.path(filt_dir, "04_manifest_filtered.csv"),
  row.names = FALSE
)

head(manifest_filtered)