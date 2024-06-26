#script to plot QV scores from yak reports


library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)


#reading yak report
yak_flye <- read_tsv('./yak_qv/yak_flye_qv/QV_flye_stats.txt') %>% filter(FR == "QV")


directory_path <- "./flye_1000G_assembly/1000g_quast_flye_samples"
subdirectories <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)
subdirectories = subdirectories[grepl('*UL*', subdirectories)==FALSE] 


#adding subdirectories to yak_flye
yak_flye$sample_name = subdirectories

yak_flye$sample_name = substr(yak_flye$sample_name,48,65)

#filtering for UL


#yak_flye = yak_flye %>% mutate(sample_type = ifelse(grepl('*UL*', sample_name), 'UL', 'L'))
#yak_flye = yak_flye %>% filter(sample_type == "L")

#SHASTA
yak_shasta_hap1 <- read_tsv('./yak_qv/yak_shasta_hapdup/QV_hap1_stats.txt') %>% filter(FR == "QV")
yak_shasta_hap2 <- read_tsv('./yak_qv/yak_shasta_hapdup/QV_hap2_stats.txt') %>% filter(FR == "QV")

#making violin plot
combined_yak_stats = data.frame(yak_flye$sample_name, yak_flye$`0`, yak_shasta_hap1$`0`, yak_shasta_hap2$`0`)
colnames(combined_yak_stats) <- c("sample_name", "flye", "hap1", "hap2")
combined_yak_stats$shasta = (combined_yak_stats$hap1 + combined_yak_stats$hap2)/2

# Pivot the data
melted_data <- pivot_longer(combined_yak_stats, cols = c(flye, shasta), names_to = "Assembler", values_to = "Values")

##############################################################
#1000G Figure
##############################################################
# Create the violin plot
ggplot(melted_data, aes(x = Assembler, y = Values, fill = Assembler)) +
  geom_violin(scale = "width", width = 0.5, trim = FALSE) +
  labs(
    
    y = "QV score") +
  scale_fill_manual(values = c("orange", "grey")) +
  theme_bw()





