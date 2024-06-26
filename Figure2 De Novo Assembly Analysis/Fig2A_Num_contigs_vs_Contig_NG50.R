library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(openxlsx)
# Set the directory path where your folders are located
directory_path <- "./flye_1000G_assembly/1000g_quast_flye_samples"

# Get a list of all subdirectories in the specified directory
subdirectories <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

# Create an empty list to store the data frames
combined_data = read.delim('./flye_1000G_assembly/1000g_quast_flye_samples/HG00105.LSK110.R9/report.tsv', header = TRUE, sep = "\t") %>% select(Assembly)


# Loop through each subdirectory
for (subdirectory in subdirectories) {
  # Extract the folder name from the subdirectory path
  folder_name <- basename(subdirectory)
  
  # Construct the file path for the reports.tsv file in the current folder
  file_path <- file.path(subdirectory, "report.tsv")
  
 
  # Read the reports.tsv file into a data frame
  #rownames(data) <- NULL
  data <- read.delim(file_path, header = TRUE, sep = "\t")
  colnames(data)[2] <- folder_name
  
  # Rename the rows with the folder name
  #rownames(data) <- data[1]
  
  # Combine the data with the existing data using cbind
  combined_data <- cbind(combined_data, data[2])
  
}


df = data.table::transpose(combined_data)
colnames(df) <- df[1,]
df$sample = colnames(combined_data)

df = df %>%
  filter(!row_number() %in% c(1))

#write.table(df, file = './flye_1000g_quast_stats.tsv', row.names = FALSE, col.names = FALSE, sep = '\t')

#label samples as "ULK" if ultra long
df = df %>% mutate(sample_type = ifelse(grepl('*UL*', sample), 'UL', 'L')) %>%
  select(NG50, `# contigs (>= 0 bp)`, sample_type)

colnames(df) <- c("assembly_N50", "contigs", "sample_type")

df$assembly_N50 = as.numeric(df$assembly_N50)
df$contigs = as.numeric(df$contigs)


#CARD samples H1

directory_path <- "./1000g_quast_card_samples/dual_1/"

# Get a list of all subdirectories in the specified directory
subdirectories <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

# Create an empty list to store the data frames
combined_card_data = read.delim('./1000g_quast_card_samples/dual_1/HG00105.LSK110.R9_hapdup_dual_1.fasta/report.tsv', header = TRUE, sep = "\t") %>% select(Assembly)


# Loop through each subdirectory
for (subdirectory in subdirectories) {
  # Extract the folder name from the subdirectory path
  folder_name <- basename(subdirectory)
  
  # Construct the file path for the reports.tsv file in the current folder
  file_path <- file.path(subdirectory, "report.tsv")
  
  
  # Read the reports.tsv file into a data frame
  #rownames(data) <- NULL
  data <- read.delim(file_path, header = TRUE, sep = "\t")
  colnames(data)[2] <- folder_name
  
  # Rename the rows with the folder name
  #rownames(data) <- data[1]
  
  # Combine the data with the existing data using cbind
  combined_card_data <- cbind(combined_card_data, data[2])
  
}


df_card = data.table::transpose(combined_card_data)
colnames(df_card) <- df_card[1,]
df_card$sample = colnames(combined_card_data)
df_card = df_card %>%
  filter(!row_number() %in% c(1))



#label samples as "ULK" if ultra long
df_card = df_card %>% mutate(sample_type = ifelse(grepl('*UL*', sample), 'UL', 'L')) %>%
  select(NG50, `# contigs (>= 0 bp)`, sample_type)

colnames(df_card) <- c("assembly_N50", "contigs", "sample_type")

df_card$assembly_N50 = as.numeric(df_card$assembly_N50)
df_card$contigs = as.numeric(df_card$contigs)



#plot contigs vs assembly N50
df$assembler = rep("Flye", nrow(df))
df_card$assembler = rep("Shasta-Hapdup", nrow(df_card))
df = rbind(df, df_card)

df = df %>% mutate(technology = paste0(assembler, "_",sample_type))

# Create a scatter plot using ggplot2 only for long sequencing
df_long = df %>% filter(sample_type == "L")
shape_mapping <- c("Flye - ONT" = 16, "Shasta-Hapdup - ONT" = 17)
  

###################### BENCHMARK
#with benchmark datasets
# Set the directory path where your folders are located
directory_path <- "./flye_1000G_assembly/1000g_quast_benchmark_datasets/"

# Get a list of all subdirectories in the specified directory
subdirectories <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

# Create an empty list to store the data frames
combined_data = read.delim('./flye_1000G_assembly/1000g_quast_flye_samples/HG00105.LSK110.R9/report.tsv', header = TRUE, sep = "\t") %>% select(Assembly)


# Loop through each subdirectory
for (subdirectory in subdirectories) {
  # Extract the folder name from the subdirectory path
  folder_name <- basename(subdirectory)
  
  # Construct the file path for the reports.tsv file in the current folder
  file_path <- file.path(subdirectory, "report.tsv")
  
  
  # Read the reports.tsv file into a data frame
  #rownames(data) <- NULL
  data <- read.delim(file_path, header = TRUE, sep = "\t")
  colnames(data)[2] <- folder_name
  
  # Rename the rows with the folder name
  #rownames(data) <- data[1]
  
  # Combine the data with the existing data using cbind
  combined_data <- cbind(combined_data, data[2])
  
}


df = data.table::transpose(combined_data)
colnames(df) <- df[1,]
df$sample = colnames(combined_data)

df = df %>%
  filter(!row_number() %in% c(1))

#write.table(df, file = './flye_1000g_quast_stats.tsv', row.names = FALSE, col.names = FALSE, sep = '\t')

#label samples as "ULK" if ultra long
df = df %>% mutate(sample_type = ifelse(grepl('*UL*', sample), 'UL', 'L')) %>%
  select(NG50, `# contigs (>= 0 bp)`, sample_type)

colnames(df) <- c("assembly_N50", "contigs", "sample_type")

df$assembly_N50 = as.numeric(df$assembly_N50)
df$contigs = as.numeric(df$contigs)


#CARD samples H1

directory_path <- "./1000g_quast_card_samples/benchmark/dual_1/"

# Get a list of all subdirectories in the specified directory
subdirectories <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

# Create an empty list to store the data frames
combined_card_data = read.delim('./1000g_quast_card_samples/dual_1/HG00105.LSK110.R9_hapdup_dual_1.fasta/report.tsv', header = TRUE, sep = "\t") %>% select(Assembly)


# Loop through each subdirectory
for (subdirectory in subdirectories) {
  # Extract the folder name from the subdirectory path
  folder_name <- basename(subdirectory)
  
  # Construct the file path for the reports.tsv file in the current folder
  file_path <- file.path(subdirectory, "report.tsv")
  
  
  # Read the reports.tsv file into a data frame
  #rownames(data) <- NULL
  data <- read.delim(file_path, header = TRUE, sep = "\t")
  colnames(data)[2] <- folder_name
  
  # Rename the rows with the folder name
  #rownames(data) <- data[1]
  
  # Combine the data with the existing data using cbind
  combined_card_data <- cbind(combined_card_data, data[2])
  
}


df_card = data.table::transpose(combined_card_data)
colnames(df_card) <- df_card[1,]
df_card$sample = colnames(combined_card_data)
df_card = df_card %>%
  filter(!row_number() %in% c(1))

#write.table(df_card, file = './shasta_hapdup_1000g_quast_stats.tsv', row.names = FALSE, col.names = FALSE, sep = '\t')


#label samples as "ULK" if ultra long
df_card = df_card %>% mutate(sample_type = ifelse(grepl('*UL*', sample), 'UL', 'L')) %>%
  select(NG50, `# contigs (>= 0 bp)`, sample_type)

colnames(df_card) <- c("assembly_N50", "contigs", "sample_type")

df_card$assembly_N50 = as.numeric(df_card$assembly_N50)
df_card$contigs = as.numeric(df_card$contigs)



#plot contigs vs assembly N50
df$assembler = rep("Flye_Benchmark", nrow(df))
df_card$assembler = rep("Shasta-Hapdup_Benchmark", nrow(df_card))
df_benchmark = rbind(df, df_card)

df_benchmark = df_benchmark %>% mutate(technology = paste0(assembler, "_",sample_type))


#### merge
df_long_benchmark = rbind(df_long, df_benchmark)


########## PLOT with BENCHMARK
ggplot(data = df_long_benchmark, aes(x = assembly_N50, y = contigs)) +
  geom_point(aes(shape = factor(assembler), fill = factor(assembler)), size = 2.5, color = "black") +
  scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6)) +
  ylab("Number of contigs") +
  xlab("Contig NG50") +
  #labs(fill = "Assembler") +
  theme_bw() +
  scale_shape_manual(values = c("Flye" = 21, "Shasta-Hapdup" = 24, "Flye_Benchmark" = 21, "Shasta-Hapdup_Benchmark" = 24), guide = FALSE) +
  scale_fill_manual(values = c("Flye" = "#D55E00", "Shasta-Hapdup" = "#FFD700", "Flye_Benchmark" = "black", "Shasta-Hapdup_Benchmark" = "#0D0301"),
                    labels = c("Flye - ONT", "Shasta-Hapdup - ONT", "Flye Benchmark", "Shasta-Hapdup Benchmark")) +
  #scale_fill_manual(values = c("Flye - ONT" = "#D55E00", "Shasta-Hapdup - ONT" = "#FFD700"),
  #                  labels = c("Flye - ONT", "Shasta-Hapdup - ONT")) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "right", 
        legend.box = "horizontal",
        legend.background = element_rect(color = "black", size = 0.2)
  ) +
  guides(fill = guide_legend(title = "Assembler", override.aes = list(shape = c(21, 24, 21, 24), fill = c("#D55E00", "#FFD700", "black", "#0D0301"))))



