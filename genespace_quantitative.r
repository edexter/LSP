library(tidyverse)
library(ggplot2)
library(zoo)  # for rolling means

#Load data
df <- read.delim("C:/Users/ericd/Desktop/combBed.txt")

# Remove scaffold genomes
exclude_list <- c("XINB_scaffold","T1_scaffold","CH434_scaffold","CH14_scaffold")
df <- df[!(df$genome %in% exclude_list),  ]

# Step 1: Create a binary indicator for presence in each genome
df <- df %>%
  mutate(presence = 1)  # add a column with all 1s indicating presence

# Step 2: Transform the data to wide format
wide_df <- df %>%
  select(id, genome, presence) %>%
  distinct() %>%
  pivot_wider(names_from = genome, values_from = presence, values_fill = 0)

wide_df$totals <- rowSums(wide_df[,-1])

###############################################################################
#Reformat and count occurances

# 1. Count number of genomes each gene was found in:
gene_presence_count <- df %>%
  group_by(og) %>%
  summarise(presence_count = sum(presence > 0)) 

# 2. Filter for reference genome (CH14_scaffold):
reference_genome_data <- df %>%
  filter(genome == "t2_17_3_4i_13",chr == 31 | chr == 19) %>%
  left_join(gene_presence_count, by = "og")

# 2. Filter for chromosome 5 and reference genome (CH14_scaffold):
cutoff <- 0.5*length(unique(df$genome))

table(reference_genome_data$presence_count > cutoff)

# Calculate rolling mean of presence count with a specified window size:
window_size <- 25  # Define the size of the window for the rolling mean
reference_genome_data <- reference_genome_data %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(rolling_mean = rollmean(presence_count, window_size, fill = NA, align = 'right'))

# 3. Plotting with sliding window mean:
ggplot(reference_genome_data, aes(x = start, y = presence_count)) +
  geom_point(alpha = 0.2) +  # Points for each gene's location
  geom_line(aes(y = rolling_mean), size = 1, linetype = "dashed", col = "red") +  # Line for the rolling mean
  labs(title = "Gene Presence and Rolling Mean Across All Chromosomes",
       x = "Start Position", 
       y = "Number of Genomes / Rolling Mean",
       color = "Chromosome") +
  theme_minimal() +
  geom_hline(yintercept = 25,col = "blue")


################################################################################
# 2. Filter for chromosome 5 and reference genome (CH14_scaffold):
cutoff <- 0.5*length(unique(df$genome))

# Stats
counts <- df %>%
  filter(genome == "t2_17_3_4i_13",chr == 191 | chr == 351) %>%
  left_join(gene_presence_count, by = "og")

table(counts$presence_count > cutoff)
table(counts$presence_count > cutoff) / nrow(counts)

# Right arm
counts <- df %>%
  filter(genome == "t2_17_3_4i_13",chr == 161 | chr == 31) %>%
  left_join(gene_presence_count, by = "og")

table(counts$presence_count > cutoff)
table(counts$presence_count > cutoff) / nrow(counts)
