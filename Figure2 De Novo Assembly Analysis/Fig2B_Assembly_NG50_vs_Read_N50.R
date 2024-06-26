#script to plot assembly stats for flye 
library(dplyr)
library(tidyverse)
library(ggplot2)

# Analysis 1: assembly N50 vs read N50

# extract assembly stats for flye and clean
flye_stats <- read.table('./cleaned_flye_assembly_stats.txt', header = FALSE) %>%
  mutate(V1 = str_sub(V1, start = 18,end = -5),
         V2 = str_sub(V2, end = -2)) %>%
  dplyr::rename(sample = V1,
                metric = V2,
                value = V3) %>%
  pivot_wider(names_from = metric) %>%
  dplyr::rename(Assembly_N50 = Fragments_N50) %>%
  separate(`Reads_N50/N90`, into = c("Read_N50", "Read_N90"), sep = "/") %>%
  mutate_at(vars(2:10), as.numeric)

flye_stats$sample = str_sub(flye_stats$sample, start = 2)

#Quast Flye Stats
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
  select(sample,NG50, `# contigs (>= 0 bp)`, sample_type) %>%
  filter(sample_type != 'UL')

colnames(df) <- c("sample","assembly_N50", "contigs", "sample_type")

df$assembly_N50 = as.numeric(df$assembly_N50)
df$contigs = as.numeric(df$contigs)


# Analysis 2: flye vs shasta assembly N50
#quast_stats <- read_tsv('./1000g_quast_stats.tsv') %>% select(sample, N50) %>% filter(sample %in% flye_stats$sample)
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
  
  # dplyr::rename the rows with the folder name
  #rownames(data) <- data[1]
  
  # Combine the data with the existing data using cbind
  combined_card_data <- cbind(combined_card_data, data[2])
  
}


df_card = t(combined_card_data)
colnames(df_card) <- df_card[1,]
df_card = as.data.frame(df_card)
df_card$sample = colnames(combined_card_data)
df_card = df_card %>%
  filter(!row_number() %in% c(1))

#label samples as "ULK" if ultra long
df_card = df_card %>% mutate(sample_type = ifelse(grepl('*UL*', sample), 'UL', 'L')) %>%
  select(sample, NG50, `# contigs (>= 0 bp)`, sample_type) %>%
  filter(sample_type != "UL")





colnames(df_card) <- c("sample","Shasta_Hapdup_N50", "contigs", "sample_type")

df_card$Shasta_Hapdup_N50 = as.numeric(df_card$Shasta_Hapdup_N50)
df_card$contigs = as.numeric(df_card$contigs)

rownames(df_card) <- NULL

df_card = df_card %>% 
  mutate(sample = str_sub(sample, end = -21),
         sample = toupper(sample))


flye_stats$sample = toupper(flye_stats$sample)

#merging flye stats and df
flye_stats = left_join(flye_stats, df, by = "sample") %>% select(sample, Read_N50, assembly_N50)


combined_stats = merge(flye_stats, df_card, by = c("sample")) %>% 
  select(sample, Read_N50, assembly_N50, Shasta_Hapdup_N50, sample_type) %>%
  dplyr::rename(Flye_N50 = assembly_N50)

combined_stats = as.data.frame(combined_stats)

combined_stats$Flye_N50 = as.numeric(combined_stats$Flye_N50)
combined_stats$Read_N50 = as.numeric(combined_stats$Read_N50)


# plotting differences between along read N50s on X axis combined_stats =
ggplot(combined_stats, aes(x = Read_N50)) +
  geom_linerange(aes(ymin = Flye_N50, ymax = Shasta_Hapdup_N50), color = "black", alpha = 0.15) +
  geom_point(aes(y = Flye_N50), color = "orange") +
  geom_point(aes(y = Shasta_Hapdup_N50), color = "blue") +
  #geom_text(data = subset(combined_stats, Shasta_Hapdup_N50 < Flye_N50), 
  #          aes(y = Shasta_Hapdup_N50, label = sample), nudge_y = -0.5, vjust = 1.5, size = 2.5, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Read N50") +
  ylab("Assembly N50") +
  #scale_x_continuous(labels = function(x) paste0(x/1e6, "Mb")) +
  scale_y_continuous(labels = function(x) paste0(x/1e6, "Mb")) 



###############
#1000G figure:
###############
long_combined_stats = combined_stats 

#reading benchmarking datasets - quast Shasta 
directory_path <- "./1000g_quast_card_samples/benchmark/dual_1/"

# Get a list of all subdirectories in the specified directory
subdirectories <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

# Create an empty list to store the data frames
combined_card_benchmark_data = read.delim('./1000g_quast_card_samples/dual_1/HG00105.LSK110.R9_hapdup_dual_1.fasta/report.tsv', header = TRUE, sep = "\t") %>% select(Assembly)


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
  
  # dplyr::rename the rows with the folder name
  #rownames(data) <- data[1]
  
  # Combine the data with the existing data using cbind
  combined_card_benchmark_data <- cbind(combined_card_benchmark_data, data[2])
  
}


df_benchmark_card = t(combined_card_benchmark_data)
colnames(df_benchmark_card) <- df_benchmark_card[1,]
df_benchmark_card = as.data.frame(df_benchmark_card)
df_benchmark_card$sample = colnames(combined_card_benchmark_data)
df_benchmark_card = df_benchmark_card %>%
  filter(!row_number() %in% c(1))

#label samples as "ULK" if ultra long
df_benchmark_card = df_benchmark_card %>% mutate(sample_type = ifelse(grepl('*UL*', sample), 'UL', 'L')) %>%
  select(sample, NG50, `# contigs (>= 0 bp)`, sample_type)



colnames(df_benchmark_card) <- c("sample","Shasta_Hapdup_N50", "contigs", "sample_type")

df_benchmark_card$Shasta_Hapdup_N50 = as.numeric(df_benchmark_card$Shasta_Hapdup_N50)
df_benchmark_card$contigs = as.numeric(df_benchmark_card$contigs)

rownames(df_benchmark_card) <- NULL

df_benchmark_card = df_benchmark_card %>% 
  mutate(sample = str_sub(sample, end = -21),
         sample = toupper(sample))

df_benchmark_card$sample_type = rep("benchmark_shasta", 5)

#adding benchmark flye
flye_benchmark_stats <- read.table('./cleaned_flye_benchmark_assembly_stats.txt', header = FALSE) %>%
  mutate(V1 = str_sub(V1, start = 28,end = -5),
         V2 = str_sub(V2, end = -2)) %>%
  dplyr::rename(sample = V1,
                metric = V2,
                value = V3) %>%
  pivot_wider(names_from = metric) %>%
  dplyr::rename(Assembly_N50 = Fragments_N50) %>%
  separate(`Reads_N50/N90`, into = c("Read_N50", "Read_N90"), sep = "/") %>%
  mutate_at(vars(2:10), as.numeric)

#adding NG50 from quast
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
  select(sample,NG50, `# contigs (>= 0 bp)`, sample_type)

colnames(df) <- c("sample","assembly_N50", "contigs", "sample_type")

df$assembly_N50 = as.numeric(df$assembly_N50)
df$contigs = as.numeric(df$contigs)

#adding correct assembly N50
flye_benchmark_stats = left_join(flye_benchmark_stats, df, by = "sample") %>% select(sample, Read_N50, assembly_N50)

#merging
combined_benchmark_stats = merge(flye_benchmark_stats, df_benchmark_card, by = c("sample")) %>% 
  select(sample, Read_N50, assembly_N50, Shasta_Hapdup_N50, sample_type) %>%
  dplyr::rename(Flye_N50 = assembly_N50)

combined_benchmark_stats = as.data.frame(combined_benchmark_stats)

combined_benchmark_stats$Flye_N50 = as.numeric(combined_benchmark_stats$Flye_N50)
combined_benchmark_stats$Read_N50 = as.numeric(combined_benchmark_stats$Read_N50)
combined_benchmark_stats$sample_type = "benchmark"

benchmark_long_combined_stats = rbind(long_combined_stats, combined_benchmark_stats)

ggplot(benchmark_long_combined_stats, aes(x = Read_N50)) +
  geom_linerange(aes(ymin = Flye_N50, ymax = Shasta_Hapdup_N50), color = "black", alpha = 0.2) +
  geom_point(aes(y = Flye_N50, fill = "Flye"), size = 2, shape = "circle" ,color = "orange") +
  geom_point(aes(y = Shasta_Hapdup_N50, fill = "Shasta-Hapdup"), size = 2, shape = "triangle", color = "#00B4CC") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Read N50") +
  ylab("Assembly N50") +
  scale_y_continuous(labels = function(x) paste0(x/1e6, "Mb")) +
  guides(fill = guide_legend(title = "Assembler", override.aes = list(shape = c(16, 17), color = c("orange", "#00B4CC"))))
#scale_color_manual(values = c("orange", "blue"), name = "Assembler", labels = c("Flye", "Shasta-Hapdup"))
#scale_fill_manual(values = c("Flye" = "#D55E00", "Shasta-Hapdup" = "#00B4CC"),
#                  name = "Assembler",
#                labels = c("Flye - ONT", "Shasta-Hapdup - ONT")) 

benchmark_long_combined_stats$sample_type = as.factor(benchmark_long_combined_stats$sample_type)


ggplot(benchmark_long_combined_stats, aes(x = Read_N50)) +
  geom_linerange(aes(ymin = Flye_N50, ymax = Shasta_Hapdup_N50), color = "black", alpha = 0.2) +
  geom_point(aes(y = Flye_N50, fill = ifelse(sample_type == "L", 'orange', 'black'), color = ifelse(sample_type == "L", 'orange', 'black')), size = 2, shape = 16) +
  geom_point(aes(y = Shasta_Hapdup_N50, fill = ifelse(sample_type == "L", '#00B4CC', '#0D0301'), color = ifelse(sample_type == "L", '#00B4CC', '#0D0301')), size = 2, shape = 17) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_shape_identity() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Read N50") +
  ylab("Assembly NG50") +
  scale_y_continuous(labels = function(x) paste0(x/1e6, "Mb")) +
  scale_color_manual("Assembler",values = c( "#00B4CC", "black","#0D0301","orange"), 
                     labels = c("Flye - 1000G", "Shasta-Hapdup - benchmark", "Flye - benchmark", "Shasta-Hapdup - 1000G")) 

