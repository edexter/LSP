# This R script makes a matrix of genome alignment dotplots using the data
# that is output from the D-genies software. Script written by Eric Dexter.

################################################################################
# Load required packages
################################################################################

#NOTE: Must also load the custom functions from the file "dotplot_matrix_utils.R"

library(pafr) # For reading and manipulating the alignments
library(gridExtra) # For plotting as a grid

################################################################################
# Get list of files to be analyzed
################################################################################

# Navigate to working directory where D-genies output is stored
setwd("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/dotplots/dgenies_output")

# List all dotplots stored in the folder
# Get all files and folders in the current working directory
folders <- list.files(getwd())

# Filter to keep only folders
folders <- folders[file.info(folders)$isdir]

################################################################################
# Function to read the input files and produce a plot for each alignment
################################################################################

generateSubplots <- function(numPlots, x_labels, y_labels, plot_titles,
                             contigLineColor, contigLineSize, highlightBoxes) {
  
  plotList <- list()
  
  for (i in 1:length(numPlots)) {
    folderName <- paste(folders[i], "/map.paf", sep="")
    ali <- read_paf(folderName)
    ali$qname <- sub("^0+", "", gsub("[^0-9]+", "", ali$qname))
    ali$tname <- sub("^0+", "", gsub("[^0-9]+", "", ali$tname))
    ali <- filter_secondary_alignments(ali)
    ali <- subset(ali, alen > 2e4 & mapq > 40)
    
    targetOrderName <- paste(folders[i], "/target.idx", sep="")
    targetOrder <- read.table(targetOrderName)
    targetOrder <- row.names(targetOrder)
    targetOrder <- sub("^0+", "", gsub("[^0-9]+", "", targetOrder))
    
    queryOrderName <- paste(folders[i], "/query.idx", sep="")
    queryOrder <- read.table(queryOrderName)
    queryOrder <- row.names(queryOrder)
    queryOrder <- sub("^0+", "", gsub("[^0-9]+", "", queryOrder))
    
    contigOrder <- list(queryOrder, targetOrder)
    
    plotList[[i]] <- customPlotFun(ali, order_by="provided", ordering=contigOrder, label_seqs=FALSE, 
                                   alignment_colour = "black", contigLineColor=contigLineColor,
                                   contigLineSize=contigLineSize, line_size = 1) +
      theme_bw() +
      coord_cartesian() +
      xlab(x_labels[i]) +
      ylab(y_labels[i]) +
      ggtitle(plot_titles[i]) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      geom_rect(data = highlightBoxes[[i]], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                color = "red", linetype = "dashed", fill = NA)
  }
  
  return(plotList)
}

################################################################################
# Apply to the function to the desired folder and plot the output
################################################################################
x_labels <- c("LSP1","LSP1","LSP1",
              "LSP2","LSP2","LSP2",
              "LSP3","LSP3","LSP3")

y_labels <- c("LSP1","LSP2","LSP3",
              "LSP1","LSP2","LSP3",
              "LSP1","LSP2","LSP3")

plot_titles <- c("Haplotype 1 vs 1","Haplotype 1 vs 2"," Haplotype 1 vs 3",
                 "Haplotype 2 vs 1","Haplotype 2 vs 2","Haplotype 2 vs 3",
                 "Haplotype 3 vs 1","Haplotype 3 vs 2","Haplotype 3 vs 3")

x_labels <- ""

y_labels <- ""

highlightBoxes <- list(
  data.frame(ymin = 13394725, ymax = 18970951, xmin = 12381521, xmax = 17957743), # First plot
  data.frame(ymin = 13394725, ymax = 18970951, xmin = 11468235, xmax = 14009082),
  data.frame(ymin = 13394725, ymax = 18970951, xmin = 11183119, xmax = 14248766),
  
  data.frame(ymin = 11468235, ymax = 14009082, xmin = 13394725, xmax = 18970951),
  data.frame(ymin = 11468235, ymax = 14009082, xmin = 12448162, xmax = 15309673),
  data.frame(ymin = 11468235, ymax = 14009082, xmin = 11183119, xmax = 14248766),
  
  data.frame(ymin = 11183119, ymax = 14248766, xmin = 13394725, xmax = 18970951),
  data.frame(ymin = 11183119, ymax = 14248766, xmin = 11468235, xmax = 14009082),
  data.frame(ymin = 11183119, ymax = 14248766, xmin = 11520838, xmax = 14003273)
)



################################################################################
# Create a list of labels for all the plots
################################################################################
results <- generateSubplots(numPlots = folders, x_labels = x_labels, y_labels = y_labels,
                            plot_titles = plot_titles, contigLineColor = "gray",
                            contigLineSize = 0.5, highlightBoxes = highlightBoxes)

# Combine the plots into a 3x3 grid
grid_plot <- do.call(grid.arrange, c(results, ncol=3, nrow=3))

# Save the combined plot
ggsave("C:/Users/ericd/Downloads/LD_chr5_rankB.png", plot = grid_plot, width = 7, height = 7, dpi = 900)
