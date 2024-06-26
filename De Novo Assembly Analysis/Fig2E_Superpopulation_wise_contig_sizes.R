#1000g figure

#script to read in information on ancestries for first 100 samples from google sheets file
#script calculates contig sizes using bed file from assembly
library(dplyr)
library(ggplot2)
library(openxlsx)
library(viridis)
library(tidyverse)


sample_sheet <- read.xlsx('./1000G ONT Sample Sheet.xlsx', sheet = 'Gus_database_active', colNames = TRUE, )

col_names <- sample_sheet[1, ]
colnames(sample_sheet) <- col_names

sample_sheet <- sample_sheet[-1, ]

rm(col_names)

sample_sheet <- sample_sheet %>%
  select(unique(colnames(.)))

#filtering FIRST_100

sample_sheet = sample_sheet %>% filter(Freeze == 'FIRST_100') %>% select(`Sample Name`, `Superpopulation code`,`Superpopulation name`)
colnames(sample_sheet) <- c("sample_name", "superpop_code", "superpop_name")



#All ancestries
# Initialize AFR_contig_sizes vector
path_to_align_flye <- "./flye_1000G_assembly/align-flye-2.24-hg38"
#AFR_contig_sizes <- numeric()
#median_AFR_contig_sizes <- numeric()
median_contig_sizes = numeric()
total_assembly_size = numeric()
total_assembly_size_filtered = numeric()

# List all directories in the align_flye folder
align_flye_dirs <- list.dirs(path = path_to_align_flye, full.names = TRUE)
align_flye_dirs <- align_flye_dirs[grepl(".LSK110.", align_flye_dirs)]

# Iterate through each directory
for (dir_path in align_flye_dirs) {
  
  # Check if any of the AFR_samples is present in the directory name
  #if (any(grepl(paste0(AFR_samples, collapse = "|"), dir_path))) {
    
    # Construct the path to asm20_assembly.bed file
  bed_file_path <- file.path(dir_path, "asm20_assembly.bed")
    
    # Check if the bed file exists
  if (file.exists(bed_file_path)) {
    # Read the bed file
    bed_data <- read.table(bed_file_path, header = FALSE)
    
    # Compute the differences between column V3-V2 for all rows
    differences <- bed_data %>% group_by(V4) %>% 
      summarize(start = min(V2), stop = max(V3)) %>%
      mutate(size = stop-start) %>%
      pull(size)
    
    
  }
    
    # Append the differences to AFR_contig_sizes vector
    #AFR_contig_sizes <- c(AFR_contig_sizes, differences)
    total_assembly_size = c(total_assembly_size, sum(differences))
    differences = differences[differences >= 1000000]
    total_assembly_size_filtered = c(total_assembly_size_filtered, sum(differences))
    median_contig_sizes = c(median_contig_sizes, median(differences))
   
    
  }


median_flye_contig_sizes <- data.frame(str_sub(align_flye_dirs, 44,50), median_contig_sizes)

colnames(median_flye_contig_sizes) <- c("sample_name", "median_contig_sizes")
combined_sheet = left_join(sample_sheet, median_flye_contig_sizes, by = "sample_name") %>%
  mutate(superpop_code = factor(superpop_code))


#total_assembly size
total_flye_assembly_sizes <- data.frame(str_sub(align_flye_dirs, 44,50), total_assembly_size)
colnames(total_flye_assembly_sizes) <- c("sample_name", "assembly_sizes")
combined_sheet_assembly_size = left_join(sample_sheet, total_flye_assembly_sizes, by = "sample_name") %>%
  mutate(superpop_code = factor(superpop_code))

total_flye_assembly_sizes <- data.frame(str_sub(align_flye_dirs, 44,50), total_assembly_size)
colnames(total_flye_assembly_sizes) <- c("sample_name", "assembly_sizes")
combined_sheet_assembly_size = left_join(sample_sheet, total_flye_assembly_sizes, by = "sample_name") %>%
  mutate(superpop_code = factor(superpop_code))
#SHASTA-HAPDUP

#All ancestries
# Initialize AFR_contig_sizes vector
path_to_align_shasta <- "./shasta_1000g_assembly/"
#AFR_contig_sizes <- numeric()
#median_AFR_contig_sizes <- numeric()
mat_median_sizes = numeric()
pat_median_sizes = numeric()
shasta_median_sizes = numeric()

#assembly sizes
mat_assembly_sizes = numeric()
pat_assembly_sizes = numeric()
mat_assembly_sizes_filtered = numeric()
pat_assembly_sizes_filtered = numeric()
shasta_assembly_sizes = numeric()
shasta_assembly_sizes_filtered = numeric()
# List all directories in the align_flye folder
align_shasta_dirs <- list.files(path = path_to_align_shasta, full.names = TRUE)
align_shasta_dirs <- align_shasta_dirs[grepl(".LSK110.R9", align_shasta_dirs)]

# Iterate through each directory
for (dir_path in align_shasta_dirs) {
  
  # Check if any of the AFR_samples is present in the directory name
  #if (any(grepl(paste0(AFR_samples, collapse = "|"), dir_path))) {
  
  # Construct the path to asm20_assembly.bed file
  bed_file_path <- file.path(dir_path)
  
  # Check if the bed file exists
  if (file.exists(bed_file_path)) {
    # Read the bed file
    bed_data <- read.table(bed_file_path, header = FALSE)
    
    # Compute the differences between column V3-V2 for all rows
    differences <- bed_data %>% group_by(V4) %>% 
      summarize(start = min(V2), stop = max(V3)) %>%
      mutate(size = stop-start) %>%
      pull(size)
    
    
  }
  
  # Append the differences to AFR_contig_sizes vector
  #AFR_contig_sizes <- c(AFR_contig_sizes, differences)
  shasta_assembly_sizes = c(shasta_assembly_sizes, sum(differences))
  
  differences = differences[differences >= 1000000]
  shasta_assembly_sizes_filtered = c(shasta_assembly_sizes_filtered, sum(differences))
  shasta_median_sizes = c(shasta_median_sizes, median(differences))
  
  
  
}

median_shasta_contig_sizes = data.frame(str_sub(align_shasta_dirs, 26,59), shasta_median_sizes)
colnames(median_shasta_contig_sizes) <- c("sample_name", "median_contig_sizes")

median_shasta_contig_sizes = median_shasta_contig_sizes %>%
  mutate(hap = ifelse(grepl("mat", sample_name), "mat", "pat"))

median_shasta_contig_sizes$sample_name = str_sub(median_shasta_contig_sizes$sample_name, 1,7)
median_shasta_contig_sizes = median_shasta_contig_sizes %>% group_by(sample_name) %>% summarize(mean_hap_contig_sizes = mean(median_contig_sizes))

#assembly size
total_shasta_assembly_sizes = data.frame(str_sub(align_shasta_dirs, 26,59), shasta_assembly_sizes)
colnames(total_shasta_assembly_sizes) <- c("sample_name", "assembly_sizes")

total_shasta_assembly_sizes = total_shasta_assembly_sizes %>%
  mutate(hap = ifelse(grepl("mat", sample_name), "mat", "pat"))

total_shasta_assembly_sizes$sample_name = str_sub(total_shasta_assembly_sizes$sample_name, 1,7)
total_shasta_assembly_sizes = total_shasta_assembly_sizes %>% group_by(sample_name) %>% summarize(mean_assembly_size = mean(assembly_sizes))

#contig_sizes
combined_flye_sheet = left_join(sample_sheet, median_flye_contig_sizes, by = "sample_name") %>%
  mutate(superpop_code = factor(superpop_code))

combined_shasta_sheet = left_join(sample_sheet, median_shasta_contig_sizes, by = "sample_name") %>%
  mutate(superpop_code = factor(superpop_code))

combined_sheet = merge(combined_flye_sheet, combined_shasta_sheet) %>% select(-superpop_name)
colnames(combined_sheet) <- c("sample_name", "superpop_code", "flye", "shasta")

combined_sheet = combined_sheet %>% pivot_longer(flye:shasta, 
                                names_to = "assembler",
                                values_to = "contig_sizes")

combined_sheet$superpop_assembler = paste0(combined_sheet$superpop_code, "_",combined_sheet$assembler)

combined_sheet = combined_sheet %>% arrange(superpop_assembler)

color_palette <- viridis_pal()(length(unique(combined_sheet$superpop_assembler)))
#boxplot

range(combined_sheet$contig_sizes)
combined_sheet %>%
  ggplot(width = 0.5, mapping = aes(x = superpop_assembler, y = contig_sizes, fill = superpop_assembler)) + 
  geom_boxplot() +
  stat_summary(
    fun = "median",
    geom = "point",
    #width = 0.5,
    size = 0,
    position = position_dodge(0.75)
  ) +  # Add median lines
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +  # Remove legend
  xlab("Superpopulation") +
  ylab("Contig size") +
  scale_y_continuous(limits = c(12000000, 50000000), labels = function(x) paste0(x/1e6, "Mb")) +
  scale_fill_manual(values = scales::alpha(color_palette, alpha = 0.87))   # Adjust alpha for fill
  


