library(pafr)

# Navigate to working directory where D-genies output is stored
setwd("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/dotplots/dgenies_output/6_T1_CH434")

# List all dotplots stored in the folder
# Get all files and folders in the current working directory
folders <- list.files(getwd())

# Filter to keep only folders
folders <- folders[file.info(folders)$isdir]
folders

ali <- read_paf("map.paf")
ali$qname <- sub("^0+", "", gsub("[^0-9]+", "", ali$qname))
ali$tname <- sub("^0+", "", gsub("[^0-9]+", "", ali$tname))
ali <- filter_secondary_alignments(ali)
ali <- subset(ali, alen > 1000 & mapq > 50)

# Import contig order and simplify contig names
queryOrder <- read.table("query.idx")
queryOrder <- row.names(queryOrder)
queryOrder <- sub("^0+", "", gsub("[^0-9]+", "", queryOrder))

# Import contig order and simplify contig names
targetOrder <- read.table("target.idx")
targetOrder <- row.names(targetOrder)
targetOrder <- sub("^0+", "", gsub("[^0-9]+", "", targetOrder))

# Save list of contig orders to pass to plot function
contigOrder <- list(queryOrder,targetOrder)
contigOrder <- list(23,12)
dotplot(ali, label_seqs=TRUE, order_by="provided", ordering=contigOrder) + theme_bw()

plot_coverage(ali,target=FALSE)   

# Get portion of contig that is covered
ali