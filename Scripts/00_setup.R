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
cran_pkgs <- c("Rcpp", "devtools", "ggplot2", "optparse","taxonomizr", "dplyr", "seqinr")
ensure_cran(cran_pkgs)

## Bioconductor packages
bioc_pkgs <- c("Biostrings", "ShortRead", "dada2","phyloseq")
ensure_bioc(bioc_pkgs)

# Define other functions

## Retrieve all gzipped FASTQ files
get_files <- function(results_loc) {
  dir_abspath <- normalizePath(results_loc, mustWork = TRUE)
  file_paths <- list.files(
    path = dir_abspath,
    pattern = "\\.(fastq|fq)\\.gz$",
    full.names = TRUE
  )
  return(file_paths)
}

## make manifest table
make_manifest <- function(path) {

  files <- list.files(
    path,
    pattern = "\\.(fastq|fq)\\.gz$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No .fastq.gz or .fq.gz files found in the directory.")
  }

  filenames <- basename(files)

  # Detect direction (R1 / R2 / _1 / _2)
  direction <- ifelse(
    grepl("(_R1_|\\.R1\\.|_1\\.(fastq|fq)\\.gz$)", filenames),
    "forward",
    ifelse(
      grepl("(_R2_|\\.R2\\.|_2\\.(fastq|fq)\\.gz$)", filenames),
      "reverse",
      NA
    )
  )

  if (any(is.na(direction))) {
    stop(
      paste(
        "Some files do not contain a valid forward/reverse identifier:",
        paste(filenames[is.na(direction)], collapse = ", ")
      )
    )
  }

  # Extract sampleID for different naming schemes
  sampleID <- ifelse(
    grepl("_R[12]_", filenames),
    sub("(_R[12]_.*)$", "", filenames),                    # Illumina format
    ifelse(
      grepl("\\.R[12]\\.", filenames),
      sub("(\\.R[12]\\.(fastq|fq)\\.gz)$", "", filenames), # Dot format
      sub("(_[12]\\.(fastq|fq)\\.gz)$", "", filenames)     # _1/_2 format
    )
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

## Function to identify universal section of the file names (e.g. 001_R1.fastq.gz) - not working for hard file names (made new above)
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