library(dplyr)
library(readr)


#This script outputs a mean % methylated per CpG Island for one specified chromosome (for now). The chromosome can be specified in the bedMethyl intersect file. This script can also be run
# on its own with the correct inputs



args = commandArgs(TRUE)

arg1 <- args[1]
arg2 <- args[2]

setwd(arg1)

mod_kit_header <- c("chr", "start", "stop", "mod", "score", "strand", "Nvalid_cov", "per_mod", "Nmod", "Ncanon", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall", "cpg", "sample","Haplotype")

bedfile <- paste0(arg2, ".bed")

bed <- read_delim(bedfile)

names(bed) <- mod_kit_header #label columns

#calculate the average fraction of methylated reads at each CpG island for each haplotype
perIsland <- bed %>%
  dplyr::mutate(Haplotype=ifelse(Haplotype == 1, "hp1", "hp2")) %>%
  dplyr::group_by(chr,cpg,Haplotype,sample) %>%
  dplyr::reframe(start = first(start),stop =last(stop),mean_per_mod=mean(per_mod, na.rm =TRUE),std=sd(per_mod),mean_coverage = mean(Nvalid_cov))

outfile <- paste0(arg2, "perIsland.tsv")

write_tsv(perIsland, outfile)