library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

#Read in data


# List all the files with a specific file extension (e.g., .csv)
file_names <- list.files(pattern = "\\.bed$")

# Read and combine the files into a single data frame
combined_data <- do.call(rbind, lapply(file_names, function(x) read_tsv(x, col_names = FALSE)))

mod_kit_header <- c("chr", "start", "stop", "mod", "score", "strand", "Nvalid_cov", "per_mod", "Nmod", "Ncanon", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall", 
                    "dmr_chr", "dmr_start", "dmr_stop", "dmr", "inheritance", "sample","haplotype")


names(combined_data) <- mod_kit_header

formatted_pileup <- combined_data %>%
  tidyr::unite(place, dmr_chr, dmr_start, sep = ":") %>%
  tidyr::unite(dmr_locus, place, dmr_stop, sep = "-") %>%
  dplyr::mutate(haplotype = ifelse(haplotype == 1, "hp1", "hp2"))



readr::write_tsv(formatted_pileup, "first100_BWS_PWA_pileup_20240207.tsv")


condense <- formatted_pileup %>%
  dplyr::group_by(chr,dmr,haplotype,sample,dmr_locus,inheritance) %>%
  dplyr::reframe(start = first(start),stop =last(stop),mean_per_mod=mean(per_mod, na.rm =TRUE),std=sd(per_mod),mean_coverage = mean(Nvalid_cov))


readr::write_tsv(condense, "first100_BWS_PWA_pileup_condensed_20240207.tsv")
