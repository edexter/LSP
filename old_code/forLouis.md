````
# Extract low GC content region from contig 3 (narrow)

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta ptg000003l_1:8602558-8776527 > scratch/CH14_contig3_lowGC_region_right.fasta

# Extract low GC content region from contig 3 (wide)

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta ptg000003l_1:7773500-8776527 > scratch/CH14_contig3_lowGC_region_right.fasta

# Extract low GC content region from contig 19 
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta ptg000019l_1:6478672-7624328 > scratch/CH14_contig19_lowGC_region_left.fasta


````



````
library(ggplot2)
library(gridExtra)

scaleFactor <- 10^6
windowSize <- 50000

# Entire chromosome CH14
c19 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000019l_1/rc.tsv")
cutoff <- max(c19$X..Location) - windowSize
c19 <- c19[c19$X..Location < cutoff,]
maxPos <- max(c19$X..Location)

c35 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000035l_1.tsv")
cutoff <- max(c35$X..Location) - windowSize
c35 <- c35[c35$X..Location < cutoff,]
c35$X..Location <- c35$X..Location + maxPos
maxPos <- max(c35$X..Location)

c16 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000016l_1/rc.tsv")
cutoff <- max(c16$X..Location) - windowSize
c16 <- c16[c16$X..Location < cutoff,]
c16$X..Location <- c16$X..Location + maxPos
maxPos <- max(c16$X..Location)

c3 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000003l_1.tsv")
cutoff <- max(c3$X..Location) - windowSize
c3 <- c3[c3$X..Location < cutoff,]
c3$X..Location <- c3$X..Location + maxPos

CH14 <- rbind(c19,c35,c16,c3)

repeatStart <- 7773500+maxPos
repeatStop <- 8776527+maxPos

p <- ggplot(CH14,aes(x=X..Location / scaleFactor, y=X..GC))+
  geom_rect(aes(xmin = repeatStart / scaleFactor, xmax = repeatStop / scaleFactor, ymin = Inf, ymax = -Inf, alpha = 0.01), fill = "lightblue")+
  geom_line(size=0.5)+
  geom_hline(yintercept = 40, col = "red", lty =2)+
  labs(x="Position in chromosome (Mb)", y="GC %")+
  theme_bw()
p

````



````
# Load necessary libraries
library(ggplot2)

################################################################################
# Define variables
################################################################################

# File path - replace with the path to your .stk file
file_path <- "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/repeat_models/CH14-families.stk"

# Contig ID - replace with the ID of your contig
contig_id <- "ptg000003l_1"

################################################################################
# Define functions
################################################################################
# Updated function to parse a single line and extract position information along with the repeat ID
extract_position_and_id <- function(line, repeat_id) {
  parts <- strsplit(line, split="\\s+")[[1]]
  contig_and_positions <- strsplit(parts[1], split="[:-]")[[1]]
  contig <- contig_and_positions[1]
  start <- as.numeric(contig_and_positions[2])
  end <- as.numeric(contig_and_positions[3])
  return(c(RepeatID = repeat_id, Contig = contig, Start = start, End = end))
}

# Function to read .stk file and extract data for a specific contig including repeat IDs
extract_contig_positions_with_id <- function(file_name, contig_id) {
  lines <- readLines(file_name)
  data_list <- list()
  current_repeat_id <- NA
  for (line in lines) {
    if (grepl("^#=GF ID", line)) {
      current_repeat_id <- sub("^#=GF ID\\s+", "", line)
    } else if (grepl(contig_id, line)) {
      data_list[[length(data_list) + 1]] <- extract_position_and_id(line, current_repeat_id)
    }
  }
  data <- do.call(rbind, data_list)
  colnames(data) <- c("RepeatID", "Contig", "Start", "End")
  data <- data.frame(data, stringsAsFactors = FALSE)
  data$Start <- as.numeric(as.character(data$Start))
  data$End <- as.numeric(as.character(data$End))
  return(data)
}

################################################################################
# Run functions and produce plots
################################################################################
# Extract positions for the contig, including repeat IDs
positions <- extract_contig_positions_with_id(file_path, contig_id)

# Transforming 'positions' for plotting with geom_rect
positions$ymin <- 0 # Adding a constant ymin value for all rectangles
positions$ymax <- 1 # Adding a constant ymax value for all rectangles

# Plotting with ggplot2
ggplot() +
  geom_rect(data = positions, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 0.3) +
  theme_minimal() +
  geom_line( aes(c3$X..Location, y=c3$X..GC), color = "black")+
  labs(title = paste("Repeat Elements in Contig", contig_id),
       x = "Position in Contig",
       y = "")
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_y_continuous(breaks = NULL) # Remove Y-axis breaks and labels

unique(positions$RepeatID)

ggplot() +
  geom_rect(data = positions_with_id, aes(xmin = Start, xmax = End, ymin = Inf, ymax = -Inf), fill = "blue", alpha = 0.1) +
  geom_line(data = c3, aes(x=X..Location, y=X..GC), color = "black") +
  theme_minimal() +
  labs(title = paste("Repeat Elements and Additional Data in Contig", contig_id),
       x = "Position in Contig",
       y = "Value")+
  xlim(8602558,8776527)

p <- ggplot(data = c3,aes(x=X..Location / 10^6, y=X..GC))+
  geom_line(size=0.5)+
  geom_rect(data = positions, aes(xmin = Start, xmax = End, ymin = Inf, ymax = -Inf), fill = "blue", alpha = 0.3)+
  geom_hline(yintercept = 40, col = "red", lty =2)+
  labs(x="Position in chromosome (Mb)", y="GC %")+
  theme_bw()+
  xlim(7.5,8.8)

p

````

