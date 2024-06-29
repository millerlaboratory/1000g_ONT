library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)



#Sophia Gibson
#04/22/2024
#Final Repeat Plot


#Read in data, format



data <- read.delim2("first100_strchive_count.tsv")

metadata <- read.csv("STR-disease-loci.processed.csv")


metadata_format <- metadata %>%
  mutate(same = ifelse(reference_motif_reference_orientation == pathogenic_motif_reference_orientation, TRUE, FALSE)) %>% # nolint: line_length_linter.
  dplyr::select(chrom, start_hg38, stop_hg38, gene, type, pathogenic_min, pathogenic_max, same, Inheritance, id) %>% 
  tidyr::unite(place, chrom, start_hg38, sep=":") %>%
  tidyr::unite(vntr, place, stop_hg38, sep="-")


motif <- read.delim("/Volumes/172.23.192.137/users/sgibson/Projects/1KGP-ONT/1000g_preprint/all_simple_repeats/STR-disease-loci.processed_revised.bed", header=FALSE) %>%
  tidyr::unite(place, V1, V2, sep=":") %>%
  tidyr::unite(vntr, place, V3, sep="-") %>%
  dplyr::rename(motifs = V4) %>%
  dplyr::select(vntr, motifs)



plot_data <- left_join(data, metadata_format)


##############

#Figure 5A


#Part 1, small, simple alleles

small_panel_sample <- plot_data %>%
  dplyr::filter(gene != "DMD") %>% #not being plotted
  dplyr::filter(gene != "POLG") %>% #not being plotted
  dplyr::filter(gene != "ATXN10") %>% #too complex, need to further analyze
  dplyr::filter(same == TRUE) %>% #Removes complex repeats
  dplyr::group_by(id) %>%
  summarise(max=max(count)) %>%
  dplyr::filter(max <=20) %>%
  select(id)

small_plot_data <- left_join(small_panel_sample, plot_data)

small_plot_data$id_pathogenic_min <- paste(small_plot_data$id, small_plot_data$pathogenic_min, sep="_")




numeric_labeller <- function(variable){
  value <- as.character(variable)
  value <- sapply(strsplit(value, "_"), `[`, 3)
  return(value)
}

left <- ggplot(small_plot_data, aes(x=id, y=count)) +
  geom_violin(aes(fill=Inheritance)) +
  scale_fill_brewer(palette = "Accent") +
  coord_flip()+
  facet_wrap(~ id_pathogenic_min, ncol=2, scales = "free_y", strip.position = "right", labeller = labeller(id_pathogenic_min = numeric_labeller))+
  ylab("Number of Motifs")+
  theme_light()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        strip.text.y.right= element_text(angle = 0, vjust = 0.5, color = "black"),
        strip.background.y = element_rect(fill = "white", color = "grey"))





large_panel_sample <- plot_data %>%
  dplyr::filter(same == TRUE) %>%
  dplyr::group_by(id) %>%
  summarise(max=max(count)) %>%
  dplyr::filter(max > 20) %>%
  filter(id != "SCA27B_FGF14") %>% #Large, complex plot
  filter(id != "SCA10_ATXN10") %>% #Large, complex plot
  select(id)

large_plot_data <- left_join(large_panel_sample, plot_data)

large_plot_data$Inheritance <- factor(large_plot_data$Inheritance, levels = c("AD", "AD/AR", "AR", "XR", "XD"))
large_plot_data$id_pathogenic_min <- paste(large_plot_data$id, large_plot_data$pathogenic_min, sep="_")

right <- ggplot(large_plot_data, aes(x=id, y=count)) +
  geom_violin(aes(fill=Inheritance)) +
  scale_fill_brewer(palette = "Accent") +
  coord_flip()+
  facet_wrap(~ id_pathogenic_min, ncol=2, scales = "free_y", strip.position = "right", labeller = labeller(id_pathogenic_min = numeric_labeller))+
  ylab("Number of Motifs")+
  theme_light()+
  theme(
    axis.title.y=element_blank(),
    strip.text.y.right= element_text(angle = 0, vjust = 0.5, color = "black"),
    strip.background.y = element_rect(fill = "white", color = "grey"),
    axis.title.x = element_blank())



fig5A <- cowplot::plot_grid(left, right, ncol=2)

library(svglite)

svglite(file = "fig5A_simple_STR.svg", width = 15, height = 5)
print(fig5A)
dev.off()






#Part 2, large and/or complex alleles



complex_count_overall <- plot_data %>%
  dplyr::filter(same != TRUE)

ATXN10_FGF15 <- plot_data %>%
  dplyr::filter(id == "SCA27B_FGF14" | id == "SCA10_ATXN10")

complex_count_overall <- rbind(complex_count_overall, ATXN10_FGF15)


complex_count_overall$id_pathogenic_min <- paste(complex_count_overall$id, complex_count_overall$pathogenic_min, sep="_")




small_complex_sample <- complex_count_overall %>%
  dplyr::group_by(id) %>%
  summarise(max=max(count)) %>%
  dplyr::filter(max <=100) %>%
  select(id)

small_complex <- left_join(small_complex_sample, complex_count_overall)




numeric_labeller <- function(variable){
  value <- as.character(variable)
  value <- sapply(strsplit(value, "_"), `[`, 3)
  return(value)
}

top <- ggplot(small_complex, aes(x=id, y=count)) +
  geom_violin(aes(fill=Inheritance)) +
  scale_fill_brewer(palette = "Accent") +
  coord_flip()+
  facet_wrap(~ id_pathogenic_min, ncol=2, scales = "free_y", strip.position = "right", labeller = labeller(id_pathogenic_min = numeric_labeller))+
  ylab("Number of Motifs")+
  theme_light()+
  theme(legend.position = "none",axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        strip.text.y.right= element_text(angle = 0, vjust = 0.5, color = "black"),
        strip.background.y = element_rect(fill = "white", color = "grey"))






large_complex_sample <- complex_count_overall %>%
  dplyr::group_by(id) %>%
  summarise(max=max(count)) %>%
  dplyr::filter(max > 100) %>%
  select(id)

large_complex <- left_join(large_complex_sample, complex_count_overall)


bottom <- ggplot(large_complex, aes(x=id, y=count)) +
  geom_violin(aes(fill=Inheritance)) +
  scale_fill_brewer(palette = "Accent") +
  coord_flip()+
  facet_wrap(~ id_pathogenic_min, ncol=2, scales = "free_y", strip.position = "right", labeller = labeller(id_pathogenic_min = numeric_labeller))+
  ylab("Number of Motifs")+
  theme_light()+
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        strip.text.y.right= element_text(angle = 0, vjust = 0.5, color = "black"),
        strip.background.y = element_rect(fill = "white", color = "grey"))



fig5b <- cowplot::plot_grid(top, bottom, nrow=2, rel_heights = c(1, 2))

library(svglite)

svglite(file = "fig5_large_and_complex_STR_with_RFC1.svg", width = 15, height = 5)
print(fig5b)
dev.off()




#A select few of these loci are shown in the final figure


##############
#Figure 5B

motif_number <- motif %>%
  tidyr::separate_rows(motifs, sep = ",") %>%
  dplyr::group_by(vntr) %>%
  dplyr::mutate(anno=row_number()-1)

motif_number$anno <- as.character(motif_number$anno)

RFC1_alleles <- plot_data %>%
  dplyr::filter(id == "CANVAS_RFC1") %>%
  tidyr::separate_rows(anno, sep = ",") %>%
  group_by(sample,allele) %>%
  mutate(pos=row_number()) %>%
  unite(sample_allele, sample, allele, sep = "_")



RFC1_with_motif <- left_join(RFC1_alleles, motif_number)


RFC1_motif_counts <- RFC1_with_motif %>%
  group_by(sample_allele, motifs) %>%
  summarise(count = n(), .groups = 'drop') %>%
  spread(motifs, count, fill = 0) %>%
  dplyr::mutate(count = AAAAG + AAAGG + AAGAG + AGAGG +AAGGG)



heatmap_format <- tidyr::pivot_longer(RFC1_motif_counts, cols = starts_with("A"), names_to = "motif", values_to = "num") %>%
  dplyr::mutate(frac=num/count)


heatmap_format$sample_allele <- reorder(heatmap_format$sample_allele, heatmap_format$count)
heatmap_format$motif <- factor(heatmap_format$motif, levels = c("AAAAG", "AAAGG", "AAGAG", "AGAGG", "AAGGG"))


main_plot <- ggplot(heatmap_format, aes(x=motif,y=sample_allele, fill=frac)) +
  geom_tile(color='white') +
  scale_fill_gradient2(low="#D8DBD9", mid = "#2FB66C", high="#176AB3", midpoint = 0.5) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# Create the barplot
barplot <- ggplot(heatmap_format, aes(x=sample_allele, y=count)) +
  geom_bar(stat="identity") +
  theme_void() +
  coord_flip()

# Convert the barplot to a grob
barplot_grob <- ggplotGrob(barplot)

# Add the barplot to the main plot
heatmap <- main_plot + annotation_custom(grob = barplot_grob, xmin = 6, xmax = Inf, ymin = -Inf, ymax = Inf)



library(svglite)

svglite(file = "fig5_RFC1_heatmap_new_colors.svg", width = 5, height = 10)
print(heatmap)
dev.off()






#Fig 5C

example_samples <- RFC1_with_motif %>%
  dplyr::filter(sample_allele=="HG02409_hp2" | sample_allele=="HG02409_hp1" |
                  sample_allele=="HG00105_hp2" | sample_allele=="HG00105_hp1" |
                  sample_allele=="HG01122_hp1" | sample_allele=="HG01122_hp2" |
                  sample_allele=="HG00331_hp1" | sample_allele=="HG00331_hp2" |
                  sample_allele=="HG01862_hp1" | sample_allele=="HG01862_hp2")



example_samples$sample_allele <- factor(example_samples$sample_allele, levels = c("HG01862_hp2", "HG01862_hp1", "HG00331_hp2","HG00331_hp1" ,"HG01122_hp2",
                                                                                  "HG01122_hp1", "HG00105_hp1", "HG00105_hp2", "HG02409_hp1","HG02409_hp2"))

example_samples$motifs <- factor(example_samples$motifs, levels = c("AAAAG", "AAAGG", "AAGAG", "AGAGG", "AAGGG", "ACAGG"))


rcf1_large <- ggplot(example_samples, aes(x=pos,y=sample_allele, fill=motifs)) +
  geom_tile() +
  scale_fill_brewer(palette = "Accent")+
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        axis.title.y=element_blank())




library(svglite)

svglite(file = "fig5_RFC1_large_alleles.svg", width = 8, height = 2)
print(rfc1_large)
dev.off()







#Fig 5D

ATXN10_alleles <- plot_data %>%
  dplyr::filter(id == "SCA10_ATXN10") %>%
  tidyr::separate_rows(anno, sep = ",") %>%
  group_by(sample,allele) %>%
  mutate(pos=row_number()) %>%
  unite(sample_allele, sample, allele, sep = "_")



ATXN10_with_motif <- left_join(ATXN10_alleles, motif_number)







large_atxn10 <- ATXN10_with_motif %>%
  dplyr::filter(sample_allele == "HG01122_hp1" | sample_allele == "HG01122_hp2" | sample_allele == "HG02252_hp1" | sample_allele == "HG02252_hp2" |
                  sample_allele == "HG02345_hp1" | sample_allele == "HG02345_hp2")



large_atxn10$motifs <- factor(large_atxn10$motifs, levels = c("ATTCT", "ATTCC", "ATTTCT", "ATTCCT"))

large_atxn10$sample_allele <- factor(large_atxn10$sample_allele, levels = c("HG02345_hp2","HG02345_hp1","HG02252_hp2","HG02252_hp1",
                                                                            "HG01122_hp1", "HG01122_hp2"))


ATXN10_alleles <- ggplot(large_atxn10, aes(x=pos,y=sample_allele, fill=motifs)) +
  geom_tile() +
  scale_fill_brewer(palette = "Accent")+
  geom_vline(xintercept = 800) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        axis.title.y=element_blank())



library(svglite)

svglite(file = "fig5_ATXN10_alleles.svg", width = 8, height = 2)
print(ATXN10_alleles)
dev.off()



#Expansion Hunter results which are added in Illustrator


STRchive_EH_loci_id <- read_delim("/Volumes/172.23.192.137/users/sgibson/Projects/1KGP-ONT/1000g_preprint/EH_STR/STRchive_EH_loci_id.bed", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_names = FALSE, trim_ws = TRUE)


header <- c("chr", "start", "end", "format", "strchive_chr", "strchive_start", "strchive_end", "id")


names(STRchive_EH_loci_id) <- header


format <- STRchive_EH_loci_id %>%
  tidyr::separate_rows(format, sep = ",") %>%
  tidyr::separate(format, into = c("sample", "REPCN", "REPCI"), sep="=") %>%
  tidyr::drop_na() 

format$sample <- substr(format$sample, start = 1, stop = 7)


plot <- format %>%
  dplyr::mutate(hp = "hp1/hp2") %>%
  tidyr::separate_rows(REPCN,REPCI,hp, sep = "/")

RFC1 <- plot %>%
  dplyr::filter(id == "CANVAS_RFC1") %>%
  dplyr::filter(sample == "HG02409" | sample == "HG00105" | sample == "HG01122" | sample == "HG00331" | sample == "HG01862")


RFC1$sample <- factor(RFC1$sample, levels = c("HG01862", "HG00331", "HG01122", "HG00105", "HG02409"))

RFC1$REPCN <- as.numeric(RFC1$REPCN)


ci <- RFC1 %>%
  tidyr::separate(REPCI, into = c("low","high"), sep = "-")

ci$sample_allele <- paste(ci$sample, ci$hp, sep = "_")


ci$high <- as.numeric(ci$high)

ci$low <- as.numeric(ci$low)

ci$REPCN <- as.numeric(ci$REPCN)


RCF1_plot <- ggplot(ci, aes(x=sample, y=REPCN, fill=hp)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=low, ymax=high), position = position_dodge2(width = 0.2)) +
  scale_fill_brewer(palette="Accent") +
  coord_flip() +
  theme_light() +
  theme(axis.title = element_blank())

library(svglite)

svglite(file = "Fig5_RFC1_EH_CI.svg", width = 7, height = 4)
print(RFC1_plot)
dev.off()





ATXN10 <- plot %>%
  dplyr::filter(id == "SCA10_ATXN10") %>%
  dplyr::filter(sample == "HG01122" | sample == "HG02252" | sample == "HG02345")


ATXN10$sample <- factor(ATXN10$sample, levels = c("HG02345", "HG02252", "HG01122"))

ATXN10$REPCN <- as.numeric(ATXN10$REPCN)



ci_A <- ATXN10 %>%
  tidyr::separate(REPCI, into = c("low","high"), sep = "-")

ci_A$high <- as.numeric(ci_A$high)

ci_A$low <- as.numeric(ci_A$low)

ci_A$sample_allele <- paste(ci_A$sample, ci_A$hp, sep = "_")


ATXN10_plot <- ggplot(ci_A, aes(x=sample, y=REPCN, fill=hp)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=low, ymax=high), position = position_dodge2(width = 0.2)) +
  scale_fill_brewer(palette="Accent") +
  coord_flip() +
  theme_light() +
  theme(axis.title = element_blank())

library(svglite)

svglite(file = "Fig5_ATXN10_EH_CI.svg", width = 7, height = 3)
print(ATXN10_plot)
dev.off()


















