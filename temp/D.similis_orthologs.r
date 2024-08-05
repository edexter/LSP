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

# Transform the data to wide format
wide_df <- df %>%
  select(og, genome, presence) %>%
  distinct() %>%
  pivot_wider(names_from = genome, values_from = presence, values_fill = 0)

# Count the number of genomes each gene was found in
wide_df$totals <- rowSums(wide_df[,-1])

################################################################################
#Count the number of orthologs in HDH contigs in D. similis

counts <- df %>%
  filter(genome == "D_similis_IL_SIM_A20",chr == 4 | chr == 93) %>%
  left_join(wide_df, by = "og")

res <- wide_df[ wide_df$og %in% counts$og,]
totals <- colSums(res[,-1])

results <- data.frame(total = totals)

################################################################################
#Count the number of orthologs in non-HDH contigs in D. similis

counts <- df %>%
  filter(genome == "D_similis_IL_SIM_A20",chr == 134 | 
           chr == 98 |
           chr == 482 |
           chr == 286 |
           chr == 158 |
           chr == 11 |
           chr == 81 
           ) %>%
  left_join(wide_df, by = "og")

res <- wide_df[ wide_df$og %in% counts$og,]
totals <- colSums(res[,-1])

results2 <- data.frame(total = totals)