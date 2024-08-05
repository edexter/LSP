# Calculate linkage disequilibrium (LD)

Here LD is calculated across D. magna chromosome 5 relative to reads mapped against the CH14 genome. Because the chromosome is broken into 4 contigs in the assembly, LD has to be calculated within and between all of the contigs separately. The results are then merged together in the subsequent R script. 



### LD calculations

This is a small task that can be performed interactively in Linux.

````bash
#Request interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

#Load BCFtools to fix missing genotype encoding problem introduced in current release of GATK
module purge
module load BCFtools

#Recode anything with less than DP 5 as missing (required since GATK update)
bcftools +setGT /scicore/home/ebertd/dexter0000/interlink/vcfsDaphnia_CH14/merged_daphnia_filtered.vcf -- -t q -n . -i 'FORMAT/DP<5' > CH14_filtered_fixed.vcf

#Load VCFtools for LD calculation
module purge
module load VCFtools

#Calculate LD (first within then across contigs - same settings for all)
vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig3 --chr "ptg000003l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig19 --chr "ptg000019l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig35 --chr "ptg000035l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig16 --chr "ptg000016l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --interchrom-geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14 --chr "ptg000016l_1" --chr "ptg000003l_1" --chr "ptg000035l_1" --chr "ptg000019l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30
````



### LD plots

The following R script performs some post-processing of the LD calculation files and produces a plot showing the final results.

````R
# Load required packages
library(ggplot2)

# Provide contig lengths (calculated from the CH14 genome assembly FASTA file)
contig3length <-  8776528
contig16length <- 1929634
contig19length <- 7624329
contig35length <-  2306407

# Load LD data for contig 19 and flip orientation
df1 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LD/CH14_contig19.geno.ld")
maxPos <- contig19length
df1$POS1 <- maxPos - df1$POS1
df1$POS2 <- maxPos - df1$POS2

# Load LD data for contig 3 (do not flip)
df2 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LD/CH14_contig3.geno.ld")

# Load LD data for contig 16 and flip orientation
df3 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LD/CH14_contig16.geno.ld")
maxPos <- contig16length
df3$POS1 <- maxPos - df3$POS1
df3$POS2 <- maxPos - df3$POS2

# Load LD data for contig 35 (do not flip)
df4 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LD/CH14_contig35.geno.ld")

# Merge together LD data from all 4 contigs
df <- rbind(df1,df2,df3,df4)

# Recalculate variant positions to position within chromosome (not within contig)
df$POS1 <- ifelse(df$CHR=="ptg000019l_1", df$POS1,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000019l_1", df$POS2,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000016l_1", df$POS1+contig19length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000016l_1", df$POS2+contig19length,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000003l_1", df$POS1+contig19length+contig16length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000003l_1", df$POS2+contig19length+contig16length,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000035l_1", df$POS1+contig19length+contig16length+contig3length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000035l_1", df$POS2+contig19length+contig16length+contig3length,df$POS2)

# Load between contig LD data and perform the same flipping and recaluation of positions as above
df5 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LD/CH14.interchrom.geno.ld")

df5$POS1 <- ifelse(df5$CHR1=="ptg000019l_1", contig19length - df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000019l_1", contig19length - df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000016l_1", contig16length - df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000016l_1", contig16length - df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000019l_1", df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000019l_1", df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000035l_1", df5$POS1+contig19length+contig16length+contig3length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000035l_1", df5$POS2+contig19length+contig16length+contig3length,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000016l_1", df5$POS1+contig19length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000016l_1", df5$POS2+contig19length,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000003l_1", df5$POS1+contig19length+contig16length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000003l_1", df5$POS2+contig19length+contig16length,df5$POS2)

# Refomat to match the within-contig LD data
df5<-df5[,-3]
colnames(df5)[1]<-"CHR"

# Concatenate all of the within and across contig LD data
df<-rbind(df,df5)

# Reformat to make the triangular LD matrix a square matrix
dfFlip<-df
dfFlip$POS2 <-df$POS1
dfFlip$POS1 <-df$POS2
df<-rbind(df,dfFlip)

# Define plotting colors
colors <- c("gray85", "yellow", "red")

# Downsample matrix to make plotting easier (Not possible with full set)
df7 <- df[sample(1:nrow(df), 5000000), ] #subsample for testing plot

# Convert positions to Mb
df7$pos1MB <- df7$POS1/10^6
df7$pos2MB <- df7$POS2/10^6


LSPleft <- (contig16length+contig19length+1534355) / 10^6
LSPright <- (contig16length+contig19length+7111581) / 10^6

# Make the plot
p <- ggplot(df7[df7$N_INDV>1 ,], aes(x = pos1MB, y = pos2MB, fill = R.2)) +
  geom_point(alpha = 0.9, shape = 22, size = 0.3, stroke = 0)+
  coord_fixed()+
  theme_bw()+
  ggtitle("Chromosome 5 linkage disequilibrium")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Position in chromosome 5 (Mb)") + ylab("Position in chromosome 5 (Mb)")+
  theme(legend.position = c(.12, .80))+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
  scale_fill_gradientn(
    colors = colors,
    values = c(0, 0.5, 1),  # Specify where the colors change in the gradient
    limits = c(0, 1),
    guide = guide_colorbar(title = expression("LD R"^"2"), title.position = "top", title.hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(.12, .80),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())+  # Remove minor grid lines
  geom_rect(aes(xmin = LSPleft , xmax = LSPright, ymin = LSPleft, ymax = LSPright), 
            inherit.aes = FALSE, fill = NA, color = "black", linetype = "dashed", size =0.3)

# Save the plot
ggsave("C:/Users/ericd/Downloads/LD_chr5.png", plot = p, width = 6, height = 6, dpi = 900)
````



