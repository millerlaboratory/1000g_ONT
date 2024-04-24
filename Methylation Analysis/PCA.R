library(dplyr)
library(tidyr)
library(readr)
library(stats)
library(data.table)
library(scales)
library(ggplot2)

cpg.df <- fread("/Volumes/172.23.192.137/users/sgibson/MeOW_files/MeOW-main/results/meow.reference.annotated.tsv")

cpg.avg.table <- cpg.df %>% subset(select=-c(chr, pos, gene)) %>% drop_na() %>% group_by(cpg) %>% summarize_all(mean)
cpg.avg.table <- data.frame(cpg.avg.table)
row.names(cpg.avg.table) <- cpg.avg.table$cpg
cpg.avg.table <- subset(cpg.avg.table, select=-c(cpg)) %>% drop_na()


#Outlier, need to get the cpg Islands to merge with the control dataset

GM19462_cpg_methyl_indexed <- read_delim("/n/users/sgibson/MeOW_files/MeOW-main/results/GM19462.cpg.methyl.indexed.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

meow_reference_index <- read_delim("/n/users/sgibson/MeOW_files/MeOW-main/results/meow.reference.index.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::rename(index = 0)

meow_index <- meow_reference_index %>%
  dplyr::mutate(index = `0`) %>%
  select(-`0`)

GM19462_index <- left_join(GM19462_cpg_methyl_indexed, meow_index)

test <- head(GM19462_index) %>%
  rename(chr_pos = chr.pos)

bed <- GM19462_index %>%
  rename(chr_pos = chr.pos) %>%
  tidyr::separate(chr_pos, into = c("chr", "start"), sep = "\\.")

bed$start <- as.numeric(bed$start)

bed_file <- bed %>%
  dplyr::mutate(stop = start+1) %>%
  select(chr, start, stop, GM19462) %>%
  tidyr::drop_na()
#options(scipen = 999)

write.table(bed_file, "1000g_preprint/methylation/GM19462_cpg_freq.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

cpgIslands_sorted_named_core <- read_delim("MeOW_files/MeOW-main/resources/cpgIslands.sorted.named.core.tsv", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE) %>%
  select(chr, start, stop, cpg)


write.table(cpgIslands_sorted_named_core, "1000g_preprint/methylation/cpgIslands.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


GM19462_cpg_freq_islands <- read_delim("1000g_preprint/methylation/GM19462_cpg_freq_islands.bed", 
                                       delim = "\t", escape_double = FALSE, 
                                       col_names = FALSE, trim_ws = TRUE) %>%
  select(X4, X8) %>%
  tidyr::drop_na()


GM19462_cpg_mean <- GM19462_cpg_freq_islands %>%
  dplyr::rename(`...1` = X8) %>%
  dplyr::group_by(`...1`) %>%
  dplyr::reframe(GM19462 = mean(X4))




cpg_avg_table <- read_csv("/Volumes/172.23.192.137/users/sgibson/1000g_preprint/methylation/cpg.avg.table.csv")


full_cpg_avg_table <- left_join(cpg_avg_table, GM19462_cpg_mean) %>%
  tidyr::drop_na()



df <- full_cpg_avg_table[, colSums(full_cpg_avg_table != 0) > 0]

for_pca <- df %>%
  select(-c(`...1`))

for_pca <- data.frame(for_pca)

row.names(for_pca) <- df$...1

samples <- colnames(for_pca)

d.pcx <- prcomp(t(for_pca), scale=FALSE, center = TRUE)


d.mvar <- sum(d.pcx$sdev^2)
PC1.label <- paste("PC1: ", percent(sum(d.pcx$sdev[1]^2)/d.mvar, accuracy=0.1))

PC2.label <- paste("PC2: ", percent(sum(d.pcx$sdev[2]^2)/d.mvar, accuracy=0.1))
PC3.label <- paste("PC3: ", percent(sum(d.pcx$sdev[3]^2)/d.mvar, accuracy=0.1))

print(c(PC1.label, PC2.label, PC3.label))


pcx.importance <- summary(d.pcx)$importance
#extract the percent variance like so:
percent.pcx <- pcx.importance[2,]

barplot(percent.pcx[1:20])

pca.points <- data.frame(d.pcx$x)


#S1_DEMOGRAPHICS <- read_delim("methylation/gene_follow_up/S1_DEMOGRAPHICS.tsv", 
                              #delim = "\t", escape_double = FALSE, 
                              #col_names = FALSE, trim_ws = TRUE, skip = 1) %>%
  #rename(samples=X1, sex=X2, pop_code=X3, pop_name=X4, super_pop=X5, super_name=X6)

sample.list <- data.frame("samples"=samples)


sample_check <- data.frame("samples"=samples, "pca" = TRUE)

S1_DEMOGRAPHICS_and_METADATA <- read_csv("/Volumes/172.23.192.137/users/sgibson/1000g_preprint/methylation/S1_DEMOGRAPHICS_and_METADATA.csv") %>%
  dplyr::rename(samples = `Sample Name`, Superpopulation = `Superpopulation code`) 


check <- left_join(S1_DEMOGRAPHICS_and_METADATA, sample_check) %>%
  select(samples, `Seq Lab`, pca)


metadata <- left_join(sample.list, S1_DEMOGRAPHICS_and_METADATA)

pca_plot <- cbind(pca.points, metadata)


m.g.colors <- c( "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                 "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
                 "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861" )


S19 <- ggplot2::ggplot(pca_plot, aes(x=PC1, y=PC2, color=Superpopulation, shape=Sex)) +
  geom_point(size=3, stroke=1, aes(fill=Superpopulation)) +
  scale_shape_manual(values=c(21,24)) +
  #ggtitle("Super Population") +
  scale_color_manual(values=m.g.colors) +
  scale_fill_manual(values=m.g.colors) +
  #geom_label(aes(label=samples))+
  ylab(PC2.label) + xlab(PC1.label) +
  annotate("text", x = 250, y = -175, label = "GM18864", vjust=1) +
  theme_light()


library(svglite)

svglite(file = "S20_PCA.svg", width = 7, height = 5)
print(S19)
dev.off()
