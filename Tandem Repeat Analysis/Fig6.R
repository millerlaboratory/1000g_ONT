library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(tidyr)



#Fig 6B
condense <- read.delim("/Volumes/172.23.192.137/users/sgibson/Projects/1KGP-ONT/1000g_preprint/methylation/pileup_files/first100_BWS_PWA_pileup_condensed_20240207.tsv")

samples <- condense %>%
  group_by(sample) %>%
  reframe(all="all")


BWS <- condense %>%
  dplyr::filter(chr == "chr11") %>%
  dplyr::mutate(assoc = "BWS/SRS")



BWS_subset <- BWS %>%
  dplyr::filter(dmr == "H19" | dmr == "KvDMR1") %>%
  dplyr::mutate(dmr=ifelse(dmr == "KvDMR1", "KCNQ10T1", "H19"))


sample_count <- BWS_subset %>%
  dplyr::group_by(sample) %>%
  reframe(n=n())


BWS <- ggplot(BWS_subset, aes(x=dmr, y=mean_per_mod)) +
  geom_point(shape=21, position = "jitter") +
  ylab("Average Fraction Methylated")+
  #geom_line(aes(group=sample_hp)) +
  facet_grid(cols=vars(assoc), scales = "free_x") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", colour = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=16))


library(svglite)

svglite(file = "fig5B_BWS.svg", width = 4, height = 5)
print(BWS)
dev.off()





PWA_subset <- condense %>%
  dplyr::filter(chr == "chr15") %>%
  dplyr::mutate(assoc = "PWS/AS") %>%
  dplyr::filter(dmr == "SNURF") %>%
  dplyr::mutate(dmr=ifelse(dmr == "SNURF", "SNRPN/SNURF", NA))




# Assuming your samples are named "sample1" and "sample2"
selected_samples <- c("GM19473", "HG00525")

# Create a new variable that specifies whether each row is one of the two samples
PWA_subset$selected_sample <- ifelse(PWA_subset$sample %in% selected_samples, PWA_subset$sample, "Other")




PWA <- ggplot(PWA_subset, aes(x=dmr, y=mean_per_mod, fill=selected_sample)) +
  geom_point(shape=21, position = "jitter") +
  ylab("Average Fraction Methylated")+
  #geom_line(aes(group=sample_hp)) +
  facet_grid(cols=vars(assoc), scales = "free_x") +
  scale_fill_manual(values=c("red", "blue", "white"))+
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", colour = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=16))


library(svglite)

svglite(file = "fig5C_PWA.svg", width = 4, height = 5)
print(PWA)
dev.off()







off_samples <- plot_data %>%
  dplyr::filter(assoc == "PWA") %>%
  dplyr::filter(mean_per_mod > 25) %>%
  dplyr::filter(mean_per_mod < 75)

both_hap <- left_join(off_samples, plot_data, by = join_by(sample, assoc), relationship = "many-to-many")


readr::write_tsv(both_hap, "first100_PWA_over25_under75.tsv")





#6A

#CHRX

FIRST_100_XX_CpGIsland <- read.delim("/Volumes/172.23.192.137/users/sgibson/Projects/1KGP-ONT/1000g_methylation/modkit_v0.1.11/FIRST_100_XX/FIRST_100_XX_CpGIsland.tsv")



first_100_predicted_skew <- read.csv("/Volumes/172.23.192.137/users/sgibson/reference/first_100_predicted_skew.txt", sep="")

skew <- left_join(first_100_predicted_skew, FIRST_100_XX_CpGIsland) %>%
dplyr::filter(mean_coverage >= 5)


reformat <- skew %>%
group_by(chr, start, stop, cpg, sample) %>%
reframe(frac=paste(mean_per_mod, collapse = ","), std=paste(paste(std, collapse = ","))) %>%
tidyr::separate(frac, into=c("hp1_frac", "hp2_frac"), sep = ",") %>%
tidyr::separate(std, into=c("hp1_std", "hp2_std"), sep = ",")



reformat$hp1_frac <- as.numeric(reformat$hp1_frac)
reformat$hp2_frac <- as.numeric(reformat$hp2_frac)
reformat$hp1_std <- as.numeric(reformat$hp1_std)
reformat$hp2_std <- as.numeric(reformat$hp2_std)



filter <- reformat %>%
dplyr::mutate("hp1_limit"=ifelse(hp1_frac > hp2_frac, hp1_frac-hp1_std, hp1_frac+hp1_std)) %>%
dplyr::mutate("hp2_limit"=ifelse(hp2_frac > hp1_frac, hp2_frac-hp2_std, hp2_frac+hp2_std)) %>%
dplyr::mutate("overlap_1"=ifelse(hp1_frac > hp2_frac & hp2_limit > hp1_limit, TRUE, FALSE)) %>%
dplyr::mutate("overlap_2"=ifelse(hp2_frac > hp1_frac & hp1_limit > hp2_limit, TRUE, FALSE))


#Reoving cpgs that are within 1SD of each other
no_overlap <- filter %>%
dplyr::filter(overlap_1==FALSE) %>%
dplyr::filter(overlap_2==FALSE) %>%
select(sample, cpg)

count <- no_overlap %>%
group_by(cpg) %>%
reframe(n=n())

#Keeping cpgs found in 75% of skew set
informative_cpg <- count %>%
 filter(n>=16)

cpgs <- informative_cpg %>%
 select(cpg)


#readr::write_tsv(cpgs, "~/Volumes/172.23.192.137/users/sgibson/reference/1000G_informative_cpgIslandsV0.2.tsv")




#FIRST_100_XX_CpGIsland <- read.delim("/Volumes/172.23.192.137/users/sgibson/1000g_methylation/modkit_v0.1.11/FIRST_100_XX/FIRST_100_XX_CpGIsland.tsv")


#cpgs <- read_tsv("/Volumes/172.23.192.137/users/sgibson/reference/1000G_informative_cpgIslandsV0.2.tsv")

x_plot_data <- left_join(cpgs, FIRST_100_XX_CpGIsland) %>%
  filter(mean_coverage >= 5)


#readr::write_tsv(x_plot_data, "/Volumes/172.23.192.137/users/sgibson/1000g_preprint/methylation/FIRST_100_XX_PLOT_DATA.tsv")


#Set cutoff for skew
re_format <- x_plot_data %>%
  dplyr::mutate(pos=(start+stop)/2) %>%
  dplyr::group_by(cpg,sample) %>%
  dplyr::reframe(frac=paste(mean_per_mod, collapse=","), pos=paste(pos, collapse=","), hap=paste(Haplotype, collapse = ",")) %>%
  tidyr::separate(frac, into = c("frac1", "frac2"), sep=",") %>%
  tidyr::drop_na()#drop cpgs where there aren't two haplotypes with coverage over 5


re_format$frac1 <- as.numeric(re_format$frac1)
re_format$frac2 <- as.numeric(re_format$frac2)


diff <- re_format %>%
  dplyr::mutate(diff = abs(frac1-frac2)) %>%
  group_by(sample) %>%
  dplyr::reframe(median_diff = median(diff))


ggplot(diff, aes(x=median_diff, stat="identity")) +
  geom_histogram(binwidth = 1)


skewed <- diff %>%
  dplyr::filter(median_diff >= 50)









#PLOT

sample_plot <- x_plot_data %>%
  dplyr::filter(sample == "HG01801" | sample == "HG01414") %>%
  dplyr::mutate(skew = ifelse(sample == "HG01801", "Predicted Skewed", "Predicted Random"))


plot <- ggplot(sample_plot, aes(x=start, y=mean_per_mod, fill=Haplotype, shape=Haplotype)) +
  geom_point(color="black", size=3) +
  scale_fill_manual(values=c("#D95F02", "#7570B3")) +
  scale_shape_manual(values=c(21, 24)) +
  facet_grid(rows = vars(skew))+
  xlab("CpG Island Start Position (bp)") +
  ylab("Average Fraction Methylated") +
  ylim(0, 100) +
  xlim(0,156040895) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", colour = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=16))


cytoband_color <- read.delim("/Volumes/172.23.192.137/users/sgibson/Projects/github_repo/RShinyDataBrowser/Data_Explorer_V1/Data/cytoband_color.tsv") %>%
  dplyr::filter(chr=="chrX")


rect <- ggplot(cytoband_color) +
  geom_rect(aes(xmin = start, xmax = stop, ymin = 0, ymax = 0.1, fill=color), color = "black", linewidth = 0.1, alpha = 0.5) +
  theme_void() +
  scale_fill_identity()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(10, 10, 30, 10)) 



combined_plot <- (plot / rect) + plot_layout(ncol = 1, heights = c(9, .5))

combined_plot

library(svglite)

svglite(file = "fig5_chrx_methylation_triangles.svg", width = 10, height = 5)
print(combined_plot)
dev.off()



