# Set-up script - packages and functions. Run for each individual script.

# Install packages
installed <- rownames(installed.packages())

ensure_cran <- function(pkgs, repos = "https://www.stats.bris.ac.uk/R/") {
  to_install <- setdiff(pkgs, installed)
  if (length(to_install)) {
    install.packages(to_install, dependencies = TRUE, repos = repos)
  }
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  ))
}

ensure_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  to_install <- setdiff(pkgs, installed)
  if (length(to_install)) {
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  }
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  ))
}

## CRAN packages
cran_pkgs <- c("devtools", "ggplot2", "optparse")
ensure_cran(cran_pkgs)

## Bioconductor packages
bioc_pkgs <- c("Biostrings", "ShortRead", "dada2")
ensure_bioc(bioc_pkgs)

# Define functions

## Retrieve all gzipped FASTQ files
get_files <- function(results_loc) {
  dir_abspath <- normalizePath(results_loc, mustWork = TRUE)
  file_paths <- list.files(
    path = dir_abspath,
    pattern = "\\.fastq\\.gz$",
    full.names = TRUE
  )
  return(file_paths)
}

## make manifest table
make_manifest <- function(path) {

  files <- list.files(
    path,
    pattern = "\\.fastq\\.gz$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No .fastq.gz files found in the directory.")
  }

  filenames <- basename(files)

  # Detect direction (R1 / R2) in both formats
  direction <- ifelse(
    grepl("(_R1_|\\.R1\\.)", filenames),
    "forward",
    ifelse(
      grepl("(_R2_|\\.R2\\.)", filenames),
      "reverse",
      NA
    )
  )

  if (any(is.na(direction))) {
    stop("Some files do not contain R1 or R2 identifiers.")
  }

  # Extract sampleID for both naming schemes
  sampleID <- ifelse(
    grepl("_R[12]_", filenames),
    sub("(_R[12]_.*)$", "", filenames),     # Illumina format
    sub("(\\.R[12]\\.fastq\\.gz)$", "", filenames)  # Dot format
  )

  df <- data.frame(
    sampleID = sampleID,
    direction = direction,
    full_path = normalizePath(files),
    stringsAsFactors = FALSE
  )

  # Collapse to one row per sample
  result <- reshape(
    df,
    idvar = "sampleID",
    timevar = "direction",
    direction = "wide"
  )

  # Enforce requested column names
  colnames(result) <- c(
    "sampleID",
    "absolute_forward_path",
    "absolute_backward_path"
  )

  return(result)
}

## Identify universal section of the file names (e.g. 001_R1.fastq.gz)
identify_file_extension<-function(input_directory){
  extensions <- list()
  n=1
  kill="0"
  while(kill=="0"){
    curr_char=NULL
    for (input_file in basename(list.files(input_directory))){
      if (is.null(curr_char)){
        curr_char=substr(input_file,(nchar(input_file)+1)-n,nchar(input_file))
      }
      else{
	result = substr(input_file,(nchar(input_file)+1)-n,nchar(input_file))
        if (result!=curr_char){
          kill="1"
        }else{previous=substr(input_file,(nchar(input_file)+1)-n,nchar(input_file))}
      }
    }
    n=n+1
  }
  extensions$R1<-previous
  extensions$R2<-substr(input_file,(nchar(input_file)+2)-n,nchar(input_file))
  return(extensions)
}

## Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][2]

## Calculate reads for tracking table
getN <- function(x) sum(getUniques(x))