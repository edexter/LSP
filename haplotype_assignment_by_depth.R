# This script genotypes the hyper-variable region of chromosome 5, based on read
# mapping depth relative to three different reference genomes. All code written
# by Eric Dexter. Last revised on 17.01.2024

################################################################################
# SET GENOME SPECIFIC VARIABLES
################################################################################

# Assign all list of all the variables for each genome to be analyzed
sampleVariables <- list(
  CH14 = list(
    INPUT_DEPTH = "data/LSP_depth_stats/LSP_depth_CH14",
    LSP_START_POS = 2000000,
    LSP_STOP_POS = 6000000,
    HOMOZYGOTE_RATIO = 1.0,
    HETEROZYGOTE_RATIO = 0.8,
    OUTPUT_PLOT = "figures/Depth_LSP1.png",
    HAPLOTYPE_NAME = "Haplotype 1"
  ),
  CH434 = list(
    INPUT_DEPTH = "data/LSP_depth_stats/LSP_depth_CH434",
    LSP_START_POS = 180000,
    LSP_STOP_POS = 3700000,
    HOMOZYGOTE_RATIO = 1.3,
    HETEROZYGOTE_RATIO = 1.0,
    OUTPUT_PLOT = "figures/Depth_LSP2.png",
    HAPLOTYPE_NAME = "Haplotype 2"
  ),
  T1 = list(
    INPUT_DEPTH = "data/LSP_depth_stats/LSP_depth_T1",
    LSP_START_POS = 1200000,
    LSP_STOP_POS = 3500000,
    HOMOZYGOTE_RATIO = 0.85,
    HETEROZYGOTE_RATIO = 0.65,
    OUTPUT_PLOT = "figures/Depth_LSP3.png",
    HAPLOTYPE_NAME = "Haplotype 3"
))

################################################################################
# DEFINE THE ANALYSIS FUNCTIONS
################################################################################

# Define the function to set the genome
setVariablesForSample <- function(sampleName) {
  if (sampleName %in% names(sampleVariables)) {
    return(sampleVariables[[sampleName]])
  } else {
    stop("Sample not found")
  }
}

# Define the analysis function
performAnalysis <- function(vars) {
  
  # Set working directory
  setwd("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/")
  
  # Load mapping depth  
  df <- read.delim(vars$INPUT_DEPTH, header=FALSE, comment.char="#")
  
  # Load and format the sample names
  ID <- read.table("data/LSP_depth_stats/bam.list", quote="\"", comment.char="")
  ID <- gsub("_.*", "", ID$V1)
  colnames(df) <- c("Contig", "Position", ID)
  
  # Calculate mean mapping depth per site across all samples and check distribution
  siteMeans <- rowMeans(df[,-c(1:2)])
  hist(siteMeans, breaks=100)
  
  # Filter out sites with excessive depth and recheck distribution
  df <- df[siteMeans < 100,]
  siteMeans <- rowMeans(df[,-c(1:2)])
  hist(siteMeans, breaks=100)
  
  # Calculate mean mapping depth per sample
  sampleMeans <- colMeans(df[,-c(1:2)])
  hist(sampleMeans,breaks = 100)
  
  # Calculate depth over a sliding window
  Mean_slide <- slide_index_mean(x = siteMeans, i = df$Position, 
                                 before = 100000, after = 100000, na_rm = TRUE)
  
  # Label each position as "LSP", "L_flank", or "R_flank" region
  region <-ifelse(df$Position < vars$LSP_START_POS, "L_flank",
                  ifelse(df$Position > vars$LSP_STOP_POS, "R_flank","LSP"))
  
  # Calculate mean sample depth per region
  sampleMeansLSP <- colMeans(df[region == "LSP",-c(1:2)])
  sampleMeansL_flank <- colMeans(df[region == "L_flank",-c(1:2)])
  sampleMeansR_flank <- colMeans(df[region == "R_flank",-c(1:2)])
  
  # Calculate ratio of mapping depths between LSP and flanking regions
  df2 <- data.frame(cbind(sampleMeans,sampleMeansL_flank,sampleMeansLSP,sampleMeansR_flank))
  df2$LSPratio <- df2$sampleMeansLSP / ((sampleMeansR_flank+sampleMeansL_flank)/2)
  
  # Determine cutoffs to use for LSP ratios (interactive)
  plot(df2$LSPratio)+ abline(h=c(vars$HETEROZYGOTE_RATIO,vars$HOMOZYGOTE_RATIO),lty=2)
  
  # Assign genotypes
  df2$geno <- 0
  df2$geno <- ifelse(df2$LSPratio >= vars$HOMOZYGOTE_RATIO, paste("Homozygote",vars$HAPLOTYPE_NAME), df2$geno)
  df2$geno <- ifelse(df2$LSPratio < vars$HOMOZYGOTE_RATIO & df2$LSPratio >= vars$HETEROZYGOTE_RATIO , paste("Heterozgyote",vars$HAPLOTYPE_NAME), df2$geno)
  df2$geno <- ifelse(df2$LSPratio < vars$HETEROZYGOTE_RATIO, "Homozygote other LSP", df2$geno)
 
  # Produce plots 
  p <- ggplot(df2[df2$LSPratio<5,],aes(y=LSPratio, x=order(sampleMeans), color = geno ))+
    geom_point(alpha = 1)+
    theme_classic()+
    xlab("Sample Index") +ylab("Mapping ratio")+
    ggtitle(paste("Genotype by mapping depth (",vars$HAPLOTYPE_NAME,")", sep=""))+ theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = c("red", "blue", "black"))+
    theme(legend.position = "none")
  ggsave(vars$OUTPUT_PLOT,plot = p, width = 4, height = 3)
  
  # Label all samples containing LSP 1
  df2$LSP <- ifelse(df2$LSPratio >= vars$HETEROZYGOTE_RATIO, vars$HAPLOTYPE_NAME, "Unknown")
  
  # Label Homozygotes for LSP 1, Heterozygotes, and homozygotes for not LSP1
  df2$state <- ifelse(df2$LSPratio >= vars$HOMOZYGOTE_RATIO, 
                      paste("Homozygote",vars$HAPLOTYPE_NAME),
                      ifelse(df2$LSPratio < vars$HOMOZYGOTE_RATIO & df2$LSPratio >= vars$HETEROZYGOTE_RATIO, 
                             paste("Heterozgyote",vars$HAPLOTYPE_NAME),"Homozygote X"
                      ))
  
  # Add the sample IDs
  df2 <- cbind(ID,df2)
  
  # Return the results to the global environment
  return(df2)

}

################################################################################
# Run analysis for CH14 genome
################################################################################

vars <- setVariablesForSample("CH14")
CH14 <- performAnalysis(vars)

################################################################################
# Run analysis for CH434 genome
################################################################################

vars <- setVariablesForSample("CH434")
CH434 <- performAnalysis(vars)

################################################################################
# Run analysis for T1 genome
################################################################################

vars <- setVariablesForSample("T1")
T1 <- performAnalysis(vars)

################################################################################
# Merge the results together export the final genotype calls
################################################################################

# Extract the relevant fields from each data frame and merge together
ID <- CH14$ID
meanCov <- CH14$sampleMeans

stateCH14 <- CH14$state
dosCH14 <- CH14$LSP

stateCH434 <- CH434$state
dosCH434 <- CH434$LSP

stateT1 <- T1$state
dosT1 <- T1$LSP

res <- data.frame(cbind(ID, meanCov,stateCH14, stateCH434, stateT1, dosCH14, dosCH434, dosT1))

#Create a final state variable
res$state <- "Unknown"

#Assign homozygote states and verify dataframe (passes)
res$state <- ifelse(res$stateCH14 == "Homozygote Haplotype 1" & res$stateCH434 == "Homozygote X" & res$stateT1 == "Homozygote X" ,"Haplotype 1",res$state)
res$state <- ifelse(res$stateCH14 == "Homozygote X" & res$stateCH434 == "Homozygote Haplotype 2" & res$stateT1 == "Homozygote X","Haplotype 2",res$state)
res$state <- ifelse(res$stateCH14 == "Homozygote X" & res$stateCH434 == "Homozygote X" & res$stateT1 == "Homozygote Haplotype 3","Haplotype 3",res$state)

#Assign heterozygote states and verify dataframe (passes)
res$state <- ifelse(res$stateCH14 == "Heterozgyote Haplotype 1" & res$stateCH434 == "Heterozgyote Haplotype 2" & res$stateT1 == "Homozygote X","1/2",res$state)
res$state <- ifelse(res$stateCH14 == "Heterozgyote Haplotype 1" & res$stateCH434 == "Homozygote X" & res$stateT1 == "Heterozgyote Haplotype 3","1/3",res$state)
res$state <- ifelse(res$stateCH14 == "Homozygote X" & res$stateCH434 == "Heterozgyote Haplotype 2" & res$stateT1 == "Heterozgyote Haplotype 3","2/3",res$state)
res$state <- as.factor(res$state)
res$meanCov <- as.numeric(res$meanCov)
table(res$state)

#Export results
write.table(res[,c(1,2,6,7,8,9)], "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/LSP_depth_genotypes.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)


resOLD <- res
