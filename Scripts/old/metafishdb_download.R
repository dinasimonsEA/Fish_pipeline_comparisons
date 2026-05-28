#download metfish db
# load packages (install if required)
library("tidyverse")
install.packages("ape")
library("ape")

# load REMOTE references and cleaning scripts (requires internet connection)
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")

# subset reference library table by metabarcode fragment (primer set) from the following options:
reflib.sub <- subset_references(df=reflib.cleaned, metabarcode="12s.riaz")

# [OPTIONAL] taxonomically dereplicate and filter on sequence length
# 'proplen=0.5' removes sequences shorter than 50% of median sequence length
# 'proplen=0' retains all sequences
reflib.sub <- derep_filter(df=reflib.sub, derep=TRUE, proplen=0.5)

# write out reference library in various formats to current working directory
# currently supported formats are: [1] sintax, [2] dada2 (taxonomy), [3] dada2 (species), [4] Qiime2 (fasta and taxonomy), and [5] plain dbid (GenBank or BOLD database identifiers)
write_references_fasta(df=reflib.sub)