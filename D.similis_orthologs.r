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
  select(og, genome, presence) %>%
  distinct() %>%
  pivot_wider(names_from = genome, values_from = presence, values_fill = 0)

wide_df$totals <- rowSums(wide_df[,-1])



################################################################################
#Find LSP orthologs in D similis

counts <- df %>%
  filter(genome == "D_similis_IL_SIM_A20",chr == 4 | chr == 93) %>%
  left_join(gene_presence_count, by = "og")

res <- wide_df[ wide_df$og %in% counts$og,]
totals <- colSums(res[,-1])
plot(totals)

results <- data.frame(total = totals)

################################################################################


counts <- df %>%
  filter(genome == "D_similis_IL_SIM_A20",chr == 134 | 
           chr == 98 |
           chr == 482 |
           chr == 286 |
           chr == 158 |
           chr == 11 |
           chr == 81 
           ) %>%
  left_join(gene_presence_count, by = "og")


counts <- df %>%
  filter(genome == "D_similis_IL_SIM_A20",chr == 134 | 
           chr == 81 |
           chr == 158

  ) %>%
  left_join(gene_presence_count, by = "og")

res <- wide_df[ wide_df$og %in% counts$og,]
totals <- colSums(res[,-1])
plot(totals)

results <- data.frame(total = totals)
