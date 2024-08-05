# Load required packages
library(ggplot2)
library(slider)

# Define function to load and process depth data
process_depth_data <- function(depth_file, bam_file, region_limits, output_prefix, lsp_label, cutoff1 = 0.8, cutoff2 = 1.1) {
  # Import mapping depth
  df <- read.delim(depth_file, header = FALSE, comment.char = "#")
  
  # Import the sample names
  ID <- read.table(bam_file, quote = "\"", comment.char = "")
  ID <- gsub("_.*", "", ID$V1)
  
  # Populate the data frame with more useful sample names
  colnames(df) <- c("Contig", "Position", ID)
  
  # Calculate mean mapping depth per site across all samples and check distribution
  siteMeans <- rowMeans(df[,-c(1:2)])
  
  # Get rid of sites with excessive depth and recheck distribution
  df <- df[siteMeans < 100,]
  siteMeans <- rowMeans(df[,-c(1:2)])
  
  # Calculate mean mapping depth per sample
  sampleMeans <- colMeans(df[,-c(1:2)])
  
  # Calculate depth over a sliding window
  Mean_slide <- slide_index_mean(x = siteMeans, i = df$Position, before = 100000, after = 100000, na_rm = TRUE)
  
  # Label each position as "LSP", "L_flank", or "R_flank" region
  region <- ifelse(df$Position < region_limits[1], "L_flank",
                   ifelse(df$Position > region_limits[2], "R_flank", "LSP"))
  
  # Get mean sample depth per region
  sampleMeansLSP <- colMeans(df[region == "LSP", -c(1:2)])
  sampleMeansL_flank <- colMeans(df[region == "L_flank", -c(1:2)])
  sampleMeansR_flank <- colMeans(df[region == "R_flank", -c(1:2)])
  df2 <- data.frame(sampleMeans, sampleMeansL_flank, sampleMeansLSP, sampleMeansR_flank)
  df2$LSPratio <- df2$sampleMeansLSP / ((sampleMeansR_flank + sampleMeansL_flank) / 2)
  
  # Plot
  p <- ggplot(df, aes(y = Mean_slide, x = Position)) + 
    theme_classic() +
    theme(legend.position = "none") +
    annotate("rect", xmin = region_limits[1], xmax = region_limits[2], ymin = 0, ymax = Inf, alpha = 0.2, fill = "red") +
    xlab("Genomic position") + ylab("Mean depth") +
    ggtitle(paste("Contig coverage (10 Kb sliding window)", output_prefix)) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point() +
    geom_hline(yintercept = c(25, 35), lty = 1)
  
  ggsave(paste0("C:/Users/ericd/Downloads/", output_prefix, ".png"), plot = p, width = 10, height = 6)
  
  # Label all samples with specified LSPs
  df2$LSP <- ifelse(df2$LSPratio > cutoff1, lsp_label, "Unknown")
  
  # Label Homozygotes, Heterozygotes, and homozygotes for not the specified LSP
  df2$state <- ifelse(df2$LSPratio >= cutoff2, paste("Homozygote", lsp_label),
                      ifelse(df2$LSPratio < cutoff2 & df2$LSPratio > cutoff1, paste("Heterozygote", lsp_label), "Homozygote X"))
  
  # Add sample IDs
  df2 <- cbind(ID, df2)
  
  # Export for IGV
  write.table(df2, paste0("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/", output_prefix, "_igv.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  return(df2)
}

# Define function to merge haplotype results
merge_haplotype_results <- function(haplotype_files, output_file) {
  # Load first haplotype results
  merged_data <- read.delim(haplotype_files[1])
  
  # Merge subsequent haplotype results
  for (i in 2:length(haplotype_files)) {
    df <- read.delim(haplotype_files[i])
    merged_data <- merge(merged_data, df, by = "ID")
  }
  
  # Clean up the dataframe
  merged_data <- merged_data[merged_data$sampleMeans >= 5,]
  
  # Create a final state variable
  merged_data$state <- "Unknown"
  merged_data$state <- ifelse(merged_data$stateCH14 == "Homozygote A", "A", merged_data$state)
  merged_data$state <- ifelse(merged_data$stateCH434 == "Homozygote B", "B", merged_data$state)
  merged_data$state <- ifelse(merged_data$stateT1 == "Homozygote C", "C", merged_data$state)
  merged_data$state <- ifelse(merged_data$stateCH14 == "Heterozygote A" & merged_data$stateCH434 == "Heterozygote B", "AB", merged_data$state)
  merged_data$state <- ifelse(merged_data$stateCH14 == "Heterozygote A" & merged_data$stateT1 == "Heterozygote C", "AC", merged_data$state)
  merged_data$state <- ifelse(merged_data$stateCH434 == "Heterozygote B" & merged_data$stateT1 == "Heterozygote C", "BC", merged_data$state)
  
  # Export the merged data for IGV
  write.table(merged_data, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
  
  return(merged_data)
}

# Process each haplotype
haplotype1 <- process_depth_data("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/LSP_depth_CH14",
                                 "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/bam.list",
                                 c(2000000, 7000000), "CH14", "A")

haplotype2 <- process_depth_data("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/LSP_depth_CH434",
                                 "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/bam.list",
                                 c(1800000, 3700000), "CH434", "B", cutoff1 = 1.09, cutoff2 = 1.19)

haplotype3 <- process_depth_data("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/LSP_depth_T1",
                                 "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/bam.list",
                                 c(1200000, 3500000), "T1", "C", cutoff1 = 0.65, cutoff2 = 0.85)

# Merge haplotype results
merged_haplotypes <- merge_haplotype_results(c("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/CH14_igv.txt",
                                               "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/CH434_igv.txt",
                                               "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/T1_igv.txt"),
                                             "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/merged_haplotypes_igv.txt")
