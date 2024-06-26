# # script to look at breaks in regions
# 
# library(karyoploteR)
library(dplyr)
library(ggplot2)
library(patchwork)

#common files and scripts to all three plots
#ideogram 
ideogram_df <- read.table('../IBD_analyses/BAM_common_var_IBD/hgTables.csv', sep = ",", header = TRUE) %>% filter(chrom == "chr7") 

grayscale_palette <- c(
  "gneg" = "gray80",   # Color for gneg
  "gpos50" = "gray40", # Color for gpos50
  "gpos75" = "gray30", # Color for gpos75
  "gpos25" = "gray80", # Color for gpos25
  "gpos100" = "black", # Color for gpos100
  "acen" = "darkred",  # Color for acen
  "gvar" = "gray70",  # Color for gvar
  "stalk" = "gray20"  # Color for stalk
)

ideogram <- ggplot(ideogram_df) +
  geom_rect(data = ideogram_df,
            aes(xmin = chromStart, xmax = chromEnd, ymin = 0, ymax = 0.25, fill = gieStain), alpha = 1) +
  theme_void() +
  scale_fill_manual(values = grayscale_palette) +
  guides(fill = "none") +
  theme(plot.margin = margin(10, 10, 30, 10)) +
  scale_x_continuous(name = "Chromosome Position",
                     labels = scales::unit_format(unit = "Mb", scale = 1e-6)) +
  theme(axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.y = element_blank(),   # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        #axis.title.x = element_text(color = "black", size = 8),  # Show x-axis title
        axis.text.x = element_text(color = "black")) 

#segdup
segDup <- read.table('../../../../Desktop/LRS/IGVtracks/segDups.hg38.bed') %>% filter(V1 == "chr7") %>% select(V1, V2, V3)
segDup = unique(segDup)
#editing annotation based on 10 kb margin
new_segdup_annotation = segDup %>% mutate(V2 = V2 - 10000,
                                          V3 = V3 + 10000) 

new_segdup_annotation %>% arrange(V2)

gr = GenomicRanges::GRanges(seqnames = new_segdup_annotation$V1,
                            ranges = IRanges::IRanges(start = new_segdup_annotation$V2, end = new_segdup_annotation$V3))

compressed_gr = GenomicRanges::reduce(gr)

new_segdup_annotation = as.data.frame(compressed_gr) %>% select(-strand, -width) 

repeats <- read.table('../../../../Desktop/LRS/IGVtracks/repeats.hg38.sorted.bed', skip = 1) %>% filter(V1 == "chr7")


satellite_repeats = c("Satellite",
                      "Satellite/acro",
                      "Satellite/centr",
                      "Satellite/telo")

satellite_centromere_repeats = repeats %>% filter(V4 %in% satellite_repeats)
rm(repeats)

# Set the path to your folders
main_folder <- "./flye_1000G_assembly/align-flye-2.24-hg38/"

# Get a list of all subfolders within the main folder
subfolders <- list.dirs(main_folder, recursive = FALSE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each subfolder
for (subfolder in subfolders) {
  # Get the path to the bed file in the current subfolder
  bed_file <- file.path(subfolder, list.files(subfolder, pattern = "asm20_assembly.bed$", full.names = FALSE))
  
  # Read the bed file using read.table
  bed_data <- read.table(bed_file, header = FALSE, sep = "\t")
  
  
  # Concatenate columns 2 and 3 and store the result in a new column
  #bed_data$Concatenated_Columns <- paste(bed_data$V2, bed_data$V3, sep = "")
  
  bed_data = bed_data %>% filter(V5==60, V1 == "chr7") %>% group_by(V4) %>% summarize(V2 = min(V2), V3 = max(V3))
  
  # Append the modified data frame to the list
  data_list[[subfolder]] <- bed_data
}

# Combine all data frames into a single data frame
final_data <- do.call(rbind, data_list)
final_data$SampleID <- rownames(final_data)
# Extract all characters after the last /
final_data$SampleID <- gsub(".*/", "", final_data$SampleID)
final_data$contig_num <- sapply(strsplit(final_data$SampleID, "\\."), function(x) tail(x, 1))

# Remove the number after the last "."
final_data$SampleID <- gsub("\\.\\d+$", "", final_data$SampleID)

#filtering for long contigs only
final_filtered_data = final_data %>% filter((V3-V2)>=1000000)

# Filter the final data based on whether sample_name has "ULK" in it
final_L_data <- final_filtered_data %>%
  filter(!grepl("ULK", SampleID, ignore.case = TRUE))


breaks = c(final_L_data$V2, final_L_data$V3)


h <- hist(breaks, breaks = 12840, main = "Histogram of break locations for chr7", xlab = "position on x chromosome")
plot(h)

grouped_breaks_starts = final_L_data %>% group_by(V2) %>% summarize(counts = length(unique(SampleID)))
grouped_breaks_ends = final_L_data %>% group_by(V3) %>% summarize(counts = length(unique(SampleID)))
colnames(grouped_breaks_starts) = c("position", "counts")
colnames(grouped_breaks_ends) = c("position", "counts")
grouped_breaks = rbind(grouped_breaks_starts, grouped_breaks_ends) %>% arrange(position)

grouped_breaks = grouped_breaks %>% group_by(position) %>% summarize(counts = sum(counts))



chr_breaks <- ggplot(grouped_breaks, aes(x = position, y = counts, fill = position)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Contig Breaks for chr7",
       x = NULL,
       y = "Number of samples") +
  theme_bw() + theme(panel.grid = element_blank()) +
  guides(fill = FALSE) +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10))  # Adjust bottom margin to give space for labels





combined_plot <- (chr_breaks/ideogram) +
  plot_layout(ncol = 1, heights = c(1, 0.03))
# Print the combined plot
print(combined_plot)

#####################################################################################
#Highlighting based on a segDup
#####################################################################################


#satellite_centromere_repeats = repeats %>% filter(V4 == "Satellite/centr")

# Add a new column 'segdup' to grouped_breaks and initialize with 'no'
grouped_breaks$segdup <- 'no'

# Check if each position in grouped_breaks is within any segmental duplication interval
for (i in 1:nrow(grouped_breaks)) {
  position <- grouped_breaks$position[i]
  is_in_segdup <- any(new_segdup_annotation$start <= position & position <= new_segdup_annotation$end)
  grouped_breaks$segdup[i] <- ifelse(is_in_segdup, 'yes', 'no')
}

grouped_breaks$satellite_centr <- 'no'

for (i in 1:nrow(grouped_breaks)) {
  position <- grouped_breaks$position[i]
  is_in_satellite_centr <- any(satellite_centromere_repeats$V2 <= position & position <= satellite_centromere_repeats$V3)
  grouped_breaks$satellite_centr[i] <- ifelse(is_in_satellite_centr, 'yes', 'no')
}

grouped_breaks$region = "None"
grouped_breaks$region = ifelse(grouped_breaks$segdup == "yes" & grouped_breaks$satellite_centr == "no", "segdup", grouped_breaks$region)
grouped_breaks$region = ifelse(grouped_breaks$satellite_centr == "yes" & grouped_breaks$segdup == "no", "satellite_repeats", grouped_breaks$region)
grouped_breaks$region = ifelse(grouped_breaks$satellite_centr == "yes" & grouped_breaks$segdup == "yes", "segdup_and_satellite_repeats", grouped_breaks$region)


grouped_breaks$segdup = as.factor(grouped_breaks$segdup)
grouped_breaks$region = as.factor(grouped_breaks$region)

chr_breaks <- ggplot(grouped_breaks, aes(x = position, y = counts, fill = region, color = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("None" = "white", "segdup" = "white", "satellite_repeats" = "white", "segdup_and_satellite_repeats" = "white")) +  # Set colors for "no" and "yes"
  scale_color_manual(values = c("None" = "black", "segdup" = "orange", "satellite_repeats" = "red", "segdup_and_satellite_repeats" = "darkorange")) +  # Set colors for "no" and "yes"
  labs(title = "Contig Breaks for chr7 - Flye",
       x = NULL,
       y = "Number of samples") +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        # Add other theme settings as needed
  ) +
  #guides(fill = FALSE) +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10))  # Adjust bottom margin to give space for labels


combined_plot <- (chr_breaks/ideogram) +
  plot_layout(ncol = 1, heights = c(1, 0.03))
# Print the combined plot
print(combined_plot)



##########################################################################################################################################################################
##########################################################################################################################################################################
# Figure E: running breakpoint analysis for shasta h1 and h2 
##########################################################################################################################################################################
##########################################################################################################################################################################

# Set the path to your folders
main_folder <- "./shasta_1000g_assembly/"

# Get a list of all subfolders within the main folder
subfolders <- list.files(main_folder, pattern = "hapdiff_mat.bed", full.names = FALSE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each subfolder
for (subfolder in subfolders) {
  # Get the path to the bed file in the current subfolder
  bed_file <- paste0(main_folder,subfolder)
  
  # Read the bed file using read.table
  bed_data <- read.table(bed_file, header = FALSE, sep = "\t")
  
  
  # Concatenate columns 2 and 3 and store the result in a new column
  #bed_data$Concatenated_Columns <- paste(bed_data$V2, bed_data$V3, sep = "")
  
  bed_data = bed_data %>% filter(V5==60, V1 == "chr7") %>% group_by(V4) %>% summarize(V2 = min(V2), V3 = max(V3))
  
  # Append the modified data frame to the list
  data_list[[subfolder]] <- bed_data
}

# Combine all data frames into a single data frame
final_data <- do.call(rbind, data_list)
final_data$SampleID <- rownames(final_data)
# Extract all characters after the last /
final_data$SampleID <- gsub(".*/", "", final_data$SampleID)
final_data$contig_num <- sapply(strsplit(final_data$SampleID, "\\."), function(x) tail(x, 1))

# Remove the number after the last "."
final_data$SampleID <- gsub("\\.\\d+$", "", final_data$SampleID)

#filtering for long contigs only
final_filtered_data = final_data %>% filter((V3-V2)>=1000000)

# Filter the final data based on whether sample_name has "ULK" in it
final_L_data <- final_filtered_data %>%
  filter(!grepl("ULK", SampleID, ignore.case = TRUE))


breaks = c(final_L_data$V2, final_L_data$V3)


h <- hist(breaks, breaks = 12840, main = "Histogram of break locations for chr7", xlab = "position on x chromosome")
plot(h)

grouped_breaks_starts = final_L_data %>% group_by(V2) %>% summarize(counts = length(unique(SampleID)))
grouped_breaks_ends = final_L_data %>% group_by(V3) %>% summarize(counts = length(unique(SampleID)))
colnames(grouped_breaks_starts) = c("position", "counts")
colnames(grouped_breaks_ends) = c("position", "counts")
grouped_breaks = rbind(grouped_breaks_starts, grouped_breaks_ends) %>% arrange(position)

grouped_breaks = grouped_breaks %>% group_by(position) %>% summarize(counts = sum(counts))



chr_breaks2 <- ggplot(grouped_breaks, aes(x = position, y = counts, fill = position)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Contig Breaks for chr7",
       x = NULL,
       y = "Number of samples") +
  theme_bw() + theme(panel.grid = element_blank()) +
  guides(fill = FALSE) +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10))  # Adjust bottom margin to give space for labels



combined_plot2 <- (chr_breaks2/ideogram) +
  plot_layout(ncol = 1, heights = c(1, 0.03))
# Print the combined plot
print(combined_plot2)
#####################################################################################
#Highlighting based on a segDup
#####################################################################################

#segDup <- read.table('../../../../Desktop/LRS/IGVtracks/segDups.hg38.bed') %>% filter(V1 == "chr7") %>% select(V1, V2, V3)
#repeats <- read.table('../../../../Desktop/LRS/IGVtracks/repeats.hg38.sorted.bed', skip = 1) %>% filter(V1 == "chr7")

#satellite_centromere_repeats = repeats %>% filter(V4 == "Satellite/centr")

# Add a new column 'segdup' to grouped_breaks and initialize with 'no'
grouped_breaks$segdup <- 'no'

# Check if each position in grouped_breaks is within any segmental duplication interval
for (i in 1:nrow(grouped_breaks)) {
  position <- grouped_breaks$position[i]
  is_in_segdup <- any(new_segdup_annotation$start <= position & position <= new_segdup_annotation$end)
  grouped_breaks$segdup[i] <- ifelse(is_in_segdup, 'yes', 'no')
}

grouped_breaks$satellite_centr <- 'no'

for (i in 1:nrow(grouped_breaks)) {
  position <- grouped_breaks$position[i]
  is_in_satellite_centr <- any(satellite_centromere_repeats$V2 <= position & position <= satellite_centromere_repeats$V3)
  grouped_breaks$satellite_centr[i] <- ifelse(is_in_satellite_centr, 'yes', 'no')
}

grouped_breaks$region = "None"
grouped_breaks$region = ifelse(grouped_breaks$segdup == "yes" & grouped_breaks$satellite_centr == "no", "segdup", grouped_breaks$region)
grouped_breaks$region = ifelse(grouped_breaks$satellite_centr == "yes" & grouped_breaks$segdup == "no", "satellite_repeats", grouped_breaks$region)
grouped_breaks$region = ifelse(grouped_breaks$satellite_centr == "yes" & grouped_breaks$segdup == "yes", "segdup_and_satellite_repeats", grouped_breaks$region)


grouped_breaks$segdup = as.factor(grouped_breaks$segdup)
grouped_breaks$region = as.factor(grouped_breaks$region)

chr_breaks2 <- ggplot(grouped_breaks, aes(x = position, y = counts, fill = region, color = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("None" = "white", "segdup" = "white", "satellite_repeats" = "white", "segdup_and_satellite_repeats" = "white")) +  # Set colors for "no" and "yes"
  scale_color_manual(values = c("None" = "black", "segdup" = "orange", "satellite_repeats" = "red", "segdup_and_satellite_repeats" = "darkorange")) +  # Set colors for "no" and "yes"
  labs(title = "Contig Breaks for chr7 - Shasta hap1",
       x = NULL,
       y = "Number of samples") +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        # Add other theme settings as needed
  ) +
  #guides(fill = FALSE) +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10))  # Adjust bottom margin to give space for labels
# Adjust bottom margin to give space for labels


combined_plot2 <- (chr_breaks2/ideogram) +
  plot_layout(ncol = 1, heights = c(1, 0.03))
# Print the combined plot
print(combined_plot2)

grouped_breaks_hap1 <- grouped_breaks


# Set the path to your folders
main_folder <- "./shasta_1000g_assembly/"

# Get a list of all subfolders within the main folder
subfolders <- list.files(main_folder, pattern = "hapdiff_pat.bed", full.names = FALSE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each subfolder
for (subfolder in subfolders) {
  # Get the path to the bed file in the current subfolder
  bed_file <- paste0(main_folder,subfolder)
  
  # Read the bed file using read.table
  bed_data <- read.table(bed_file, header = FALSE, sep = "\t")
  
  
  # Concatenate columns 2 and 3 and store the result in a new column
  #bed_data$Concatenated_Columns <- paste(bed_data$V2, bed_data$V3, sep = "")
  
  bed_data = bed_data %>% filter(V5==60, V1 == "chr7") %>% group_by(V4) %>% summarize(V2 = min(V2), V3 = max(V3))
  
  # Append the modified data frame to the list
  data_list[[subfolder]] <- bed_data
}

# Combine all data frames into a single data frame
final_data <- do.call(rbind, data_list)
final_data$SampleID <- rownames(final_data)
# Extract all characters after the last /
final_data$SampleID <- gsub(".*/", "", final_data$SampleID)
final_data$contig_num <- sapply(strsplit(final_data$SampleID, "\\."), function(x) tail(x, 1))

# Remove the number after the last "."
final_data$SampleID <- gsub("\\.\\d+$", "", final_data$SampleID)

#filtering for long contigs only
final_filtered_data = final_data %>% filter((V3-V2)>=1000000)

# Filter the final data based on whether sample_name has "ULK" in it
final_L_data <- final_filtered_data %>%
  filter(!grepl("ULK", SampleID, ignore.case = TRUE))


breaks = c(final_L_data$V2, final_L_data$V3)


h <- hist(breaks, breaks = 12840, main = "Histogram of break locations for chr7", xlab = "position on x chromosome")
plot(h)

grouped_breaks_starts = final_L_data %>% group_by(V2) %>% summarize(counts = length(unique(SampleID)))
grouped_breaks_ends = final_L_data %>% group_by(V3) %>% summarize(counts = length(unique(SampleID)))
colnames(grouped_breaks_starts) = c("position", "counts")
colnames(grouped_breaks_ends) = c("position", "counts")
grouped_breaks = rbind(grouped_breaks_starts, grouped_breaks_ends) %>% arrange(position)

grouped_breaks = grouped_breaks %>% group_by(position) %>% summarize(counts = sum(counts))



chr_breaks3 <- ggplot(grouped_breaks, aes(x = position, y = counts, fill = position)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Contig Breaks for chr7",
       x = NULL,
       y = "Number of samples") +
  theme_bw() + theme(panel.grid = element_blank()) +
  guides(fill = FALSE) +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10))  # Adjust bottom margin to give space for labels



combined_plot3 <- (chr_breaks/ideogram) +
  plot_layout(ncol = 1, heights = c(1, 0.03))
# Print the combined plot
print(combined_plot3)

#####################################################################################
#Highlighting based on a segDup
#####################################################################################

#segDup <- read.table('../../../../Desktop/LRS/IGVtracks/segDups.hg38.bed') %>% filter(V1 == "chr7") %>% select(V1, V2, V3)
#repeats <- read.table('../../../../Desktop/LRS/IGVtracks/repeats.hg38.sorted.bed', skip = 1) %>% filter(V1 == "chr7")

#satellite_centromere_repeats = repeats %>% filter(V4 == "Satellite/centr")

# Add a new column 'segdup' to grouped_breaks and initialize with 'no'
grouped_breaks$segdup <- 'no'

# Check if each position in grouped_breaks is within any segmental duplication interval
for (i in 1:nrow(grouped_breaks)) {
  position <- grouped_breaks$position[i]
  is_in_segdup <- any(new_segdup_annotation$start <= position & position <= new_segdup_annotation$end)
  grouped_breaks$segdup[i] <- ifelse(is_in_segdup, 'yes', 'no')
}

grouped_breaks$satellite_centr <- 'no'

for (i in 1:nrow(grouped_breaks)) {
  position <- grouped_breaks$position[i]
  is_in_satellite_centr <- any(satellite_centromere_repeats$V2 <= position & position <= satellite_centromere_repeats$V3)
  grouped_breaks$satellite_centr[i] <- ifelse(is_in_satellite_centr, 'yes', 'no')
}

grouped_breaks$region = "None"
grouped_breaks$region = ifelse(grouped_breaks$segdup == "yes" & grouped_breaks$satellite_centr == "no", "segdup", grouped_breaks$region)
grouped_breaks$region = ifelse(grouped_breaks$satellite_centr == "yes" & grouped_breaks$segdup == "no", "satellite_repeats", grouped_breaks$region)
grouped_breaks$region = ifelse(grouped_breaks$satellite_centr == "yes" & grouped_breaks$segdup == "yes", "segdup_and_satellite_repeats", grouped_breaks$region)


grouped_breaks$segdup = as.factor(grouped_breaks$segdup)
grouped_breaks$region = as.factor(grouped_breaks$region)

chr_breaks3 <- ggplot(grouped_breaks, aes(x = position, y = counts, fill = region, color = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("None" = "white", "segdup" = "white", "satellite_repeats" = "white", "segdup_and_satellite_repeats" = "white")) +  # Set colors for "no" and "yes"
  scale_color_manual(values = c("None" = "black", "segdup" = "orange", "satellite_repeats" = "red", "segdup_and_satellite_repeats" = "darkorange")) +  # Set colors for "no" and "yes"
  labs(title = "Contig Breaks for chr7 - Shasta hap2",
       x = NULL,
       y = "Number of samples") +
  ylim(0,100)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        # Add other theme settings as needed
  ) +
  #guides(fill = FALSE) +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10))  # Adjust bottom margin to give space for labels


combined_plot3 <- (chr_breaks3/ideogram) +
  plot_layout(ncol = 1, heights = c(1, 0.03))
# Print the combined plot
print(combined_plot3)

grouped_breaks_hap2 <- grouped_breaks
##########################################################################################################################################################################
#########################################################  SHASTA H1+ H2 PLOT  ###############################################################################################################
##########################################################################################################################################################################
grouped_breaks_shasta = rbind(grouped_breaks_hap1, grouped_breaks_hap2)
chr_breaks_shasta <- ggplot(grouped_breaks_shasta, aes(x = position, y = counts, fill = region, color = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("None" = "white", "segdup" = "white", "satellite_repeats" = "white", "segdup_and_satellite_repeats" = "white")) +  # Set colors for "no" and "yes"
  scale_color_manual(values = c("None" = "black", "segdup" = "orange", "satellite_repeats" = "red", "segdup_and_satellite_repeats" = "darkorange")) +  # Set colors for "no" and "yes"
  labs(title = "Contig Breaks for chr7 - Shasta hap1+2",
       x = NULL,
       y = "Number of samples") +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        # Add other theme settings as needed
  ) +
  ylim(0,100)+
  #guides(fill = FALSE) +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10))  # Adjust bottom margin to give space for labels

combined_plot_shasta <- (chr_breaks_shasta/ideogram) +
  plot_layout(ncol = 1, heights = c(1, 0.03))
# Print the combined plot
print(combined_plot_shasta)

##########################################################################################################################################################################
#########################################################  COMBINING PLOTS  ###############################################################################################################
##########################################################################################################################################################################
fig2_plot <- (chr_breaks/chr_breaks_shasta/ideogram) +
  plot_layout(ncol = 1, nrow = 3, heights = c(1,1,0.1))
# Print the combined plot
print(fig2_plot)


#Shasta plot
shasta_plot <- (chr_breaks2/chr_breaks3/ideogram) +
  plot_layout(ncol = 1, nrow = 3, heights = c(1,1,0.1))
# Print the combined plot
print(shasta_plot)




big_combined_plot <- (chr_breaks/chr_breaks2/chr_breaks3/ideogram) +
  plot_layout(ncol = 1, nrow = 5, heights = c(1,1,1,0.2))
# Print the combined plot
print(big_combined_plot)


##########################################################################################################################################################################
##########################################################################################################################################################################
#########################################################  OLD CODE  ###############################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
#
# old code with rectangular segments showing contig breaks
# Example data
data1 <- read.table('./flye_1000G_assembly/align-flye-2.24-hg38/GM18856.LSK110.R9/asm20_assembly.bed', header = FALSE) %>% filter(V5 == 60)
data2 <- read.table('./flye_1000G_assembly/align-flye-2.24-hg38/GM18924.LSK110.R9/asm20_assembly.bed', header = FALSE) %>% filter(V5 == 60)

#removing duplicated rows
data1 = unique(data1)
data2 = unique(data2)


data1 = data1 %>% filter(V1 == "chr7") %>% group_by(V4) %>% summarize(V2 = min(V2), V3 = max(V3))
data1$contig_change = c(rep(seq(0:1), nrow(data1)/2),1)
d1 <- ggplot(data1, aes(x = V2, xend = V3, y = as.factor(1), yend = as.factor(1), color = as.factor(contig_change))) +
  geom_segment(size = 3) +
  scale_color_manual(values = c("orange", "yellow"), guide = FALSE) +  # Change colors as needed
  theme_void() +
  labs(x = "Position on Chromosome", y = "Group") 

data2 = data2 %>% filter(V1 == "chr7") %>% group_by(V4) %>% summarize(V2 = min(V2), V3 = max(V3))
data2$contig_change = c(rep(seq(0:1), nrow(data2)/2))
d2 <- ggplot(data2, aes(x = V2, xend = V3, y = as.factor(1), yend = as.factor(1), color = as.factor(contig_change))) +
  geom_segment(size = 3) +
  scale_color_manual(values = c("orange", "yellow"), guide = FALSE) +  # Change colors as needed
  theme_void() +
  labs(x = "Position on Chromosome", y = "Group")

d1 /d2/ideogram
