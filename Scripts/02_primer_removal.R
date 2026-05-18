# Script to remove primers. Finds and counts all primer orientations and runs cutadapt. 

# get wd & data location
path<-getwd()
test_data_loc <- paste("Data/Raw/RingTrial_Sean")

# functions

## function to check for all the possible orientations of primer
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}


message("Step 1: Identifying primers")
Sys.sleep(1)

## function to count primer occurrences
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# finds orients of selected primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# redefine filtNs based on what is actually in the dirs (sometimes an input sample is completely removed by filterAndTrim)
filt_dir <- file.path(path, "Data/Temp/filtN")

fnFs.filtN <- sort(
  file.path(
    filt_dir,
    basename(manifest$absolute_forward_path)
  )[file.exists(file.path(
    filt_dir,
    basename(manifest$absolute_forward_path)
  ))]
)

fnRs.filtN <- sort(
  file.path(
    filt_dir,
    basename(manifest$absolute_backward_path)
  )[file.exists(file.path(
    filt_dir,
    basename(manifest$absolute_backward_path)
  ))]
)

# save in case different from the filtN lists generated from the raw inputs
saveRDS(fnFs.filtN, file = paste(path, "/Data/Temp/R_objects/02_fnFs.filtN.rds", sep=""))
saveRDS(fnRs.filtN, file = paste(path, "/Data/Temp/R_objects/02_fnRs.filtN.rds", sep=""))

# count primers pre trim
pre_trim_primer_counts <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

## check
print("Table for pre-trimmed primer counts:")
print(pre_trim_primer_counts)

## write tables
write.table(pre_trim_primer_counts, paste(path, "/Data/Temp/02_pre_trim_primer_counts.tsv", sep=""), col.names=NA, sep="\t")

# now run cutadapt
message("Step 2: Running cutadapt")
Sys.sleep(1)

## create directory for cutadapt
path.cut <- file.path(paste(path, "/Data/Temp/cutadapt", sep=""))
if(!dir.exists(path.cut)) dir.create(path.cut)

# specify options needed by cutadapt
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC) # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)

print("Cutadapt flags are:")
print(R1.flags)
print(R2.flags)

## redefine the cut dirs based on what is going in
fnFs.cut <- gsub("/filtN", "/cutadapt", fnFs.filtN)
fnRs.cut <- gsub("/filtN", "/cutadapt", fnRs.filtN)

## write objects to pass to next script
saveRDS(fnFs.cut, file = paste(path, "/Data/Temp/R_objects/02_fnFs.cut.rds", sep=""))
saveRDS(fnRs.cut, file = paste(path, "/Data/Temp/R_objects/02_fnRs.cut.rds", sep=""))

## check the order of files going in/out
files_match <- identical(
  sapply(strsplit(fnFs.cut,  "cutadapt/"), `[`, 2),
  sapply(strsplit(fnFs.filtN, "filtN/"),    `[`, 2)
)

if (!files_match) {
  message("List of input files to cutadapt does not match the proposed list of output files.")
  message("Quitting.")
  quit(status = 1)
} else {
  message("Input files and filtered output files match. Continuing.")
}

# decode arguments for cutadapt
all_args <- c("-m", minimum, "-n", copies)

# run cutadapt loop
for (i in seq_along(fnFs.filtN)) {
  system2(
    cutadapt_loc,
    args = c(
      R1.flags,
      R2.flags,
      "-o", fnFs.cut[i],
      "-p", fnRs.cut[i],
      fnFs.filtN[i],
      fnRs.filtN[i],
      "-j", "0",
      "--discard-untrimmed",
      all_args
    )
  )
}

# make post trim table to check
post_trim_primer_counts <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

## check
print("Table for post-trimmed primer counts:")
print(post_trim_primer_counts)

# write post trim table
write.table(post_trim_primer_counts, paste(path, "/Data/Temp/02_post_trim_primer_counts.tsv", sep=""), col.names=NA)   

