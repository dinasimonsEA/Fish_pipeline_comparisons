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
## Identify universal section of the file names (e.g. 001_R1.fastq.gz)
guess_file_extension<-function(input_directory){
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