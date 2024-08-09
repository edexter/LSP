# Load required packages
library(tidyverse)
library(ggplot2)
library(zoo)  # for rolling means

#Load data
df <- read.delim("C:/Users/ericd/Desktop/combBed.txt")

# Remove scaffolded genome copies from the set(retain contig-level originals)
exclude_list <- c("XINB_scaffold","T1_scaffold","CH434_scaffold","CH14_scaffold")
df <- df[!(df$genome %in% exclude_list),  ]

# Create a binary indicator for presence of gene in each genome
df <- df %>%
  mutate(presence = 1)  # add a column with all 1s indicating presence

################################################################################
# CONSERVED GENE SLIDING WINDOW ACROSS CHROMOSOME 5

# Count number of genomes each gene was found in:
gene_presence_count <- df %>%
  group_by(og) %>%
  summarise(presence_count = sum(presence > 0)) 

# Filter for just the contig containing the HDH:
reference_genome_data <- df %>%
  filter(genome == "t2_17_3_4i_13",chr == 31) %>%
  left_join(gene_presence_count, by = "og")

# Calculate rolling mean of presence count with a specified window size:
window_size <- 5  # Define the size of the window for the rolling mean
reference_genome_data <- reference_genome_data %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(rolling_mean = rollmean(presence_count, window_size, fill = NA, align = 'right'))

# Plotting with sliding window mean:
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
# CONSERVED GENE COUNTS

# Set the proportion  the proportion of genomes that a gene must be found in to
# consider it conserved
cutoff <- 0.5*length(unique(df$genome))

# Count the number of conserved genes on the left chromosome 5 arm.
counts <- df %>%
  filter(genome == "t2_17_3_4i_13",chr == 191 | chr == 351) %>%
  left_join(gene_presence_count, by = "og")

# Raw numbers
table(counts$presence_count > cutoff)

# As proportion
table(counts$presence_count > cutoff) / nrow(counts)

# Count the number of conserved genes on the right chromosome 5 arm.
counts <- df %>%
  filter(genome == "t2_17_3_4i_13",chr == 161 | chr == 31) %>%
  left_join(gene_presence_count, by = "og")

# Raw numbers
table(counts$presence_count > cutoff)

# As proportion
table(counts$presence_count > cutoff) / nrow(counts)