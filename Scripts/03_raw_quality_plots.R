# Script to explore quality of sequences

# get wd
path <- getwd()

# Dataset-specific directories
output_dir <- file.path(path, "Data", "Temp", test_data_name)
rds_dir    <- file.path(output_dir, "R_objects")
results_dir <- file.path(path, "Results", test_data_name)

# ensure results directory exists
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

## read in the lists of cutadapt files
fnFs.cut <- readRDS(file = file.path(rds_dir, "02_fnFs.cut.rds"))
fnRs.cut <- readRDS(file = file.path(rds_dir, "02_fnRs.cut.rds"))

# extract sample names and write to R object
sample.names <- manifest$sampleID

print("Sample names:")
print(sample.names)

# save outputs
saveRDS(sample.names, file = file.path(rds_dir, "03_sample_names.rds"))

write.table(
  sample.names,
  file = file.path(output_dir, "03_sample_names.txt"),
  row.names = FALSE,
  col.names = FALSE
)

# write plot to file
pdf(
  file = file.path(results_dir, "03_pre_trim_quality_plots.pdf"),
  width = 10,
  height = 10
)