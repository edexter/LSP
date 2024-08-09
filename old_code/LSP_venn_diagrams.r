################################################################################
#Load required packages
library(readr)
library(stringr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)

################################################################################
# Load required data

# Load orthogroups
Orthogroups <- read.delim2("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/Orthogroups.tsv")

################################################################################
# Prepare Venn Diagrams for left arm control region
################################################################################


#LSP 1 (CH14 assembly)
################################################################################
#Parameters for this genome
CONTIG <- 191
START <- 3000000
STOP <- 7001000

#Load and format input data
BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t2_17_3_4i_13.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)
BED$contig <- BED$V1
BED_CH14 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

#Append the orthogroup to each predicted gene
TARGET <- data.frame(gsub(" ","",Orthogroups$t2_17_3_4i_13))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t2_17_3_4i_13.,",",sep = ""))

#Initialize a vector to store the results
BED_CH14$ORTH <- 0

#Run the function
for(i in 1:nrow(BED_CH14)){
  GENE <- paste(",",BED_CH14$ID[i],",",sep = "")
  BED_CH14$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

#LSP 2 (CH434 assembly)
################################################################################
#Parameters for this genome
CONTIG <- 5
START <- 613227
STOP <- 4565149

#Load and format input data
BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/CH_434-inb3-a-1.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)
BED$contig <- BED$V1
BED_CH434 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

#Append the orthogroup to each predicted gene
TARGET <- data.frame(gsub(" ","",Orthogroups$CH_434_inb3_a_1))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.CH_434_inb3_a_1.,",",sep = ""))

#Remove a few singleton genes that cause the function to return an error
BED_CH434 <- BED_CH434[BED_CH434$ID != "g5074" & BED_CH434$ID != "g5075",]

#Initialize a vector to store the results
BED_CH434$ORTH <- 0

#Run the function
for(i in 1:nrow(BED_CH434)){
  GENE <- paste(",",BED_CH434$ID[i],",",sep = "")
  BED_CH434$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

#LSP 3 (T1 assembly)
################################################################################
#Parameters for this genome
CONTIG <- 5
START <- 613227
STOP <- 4565149

#Load and format input data
BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t1_10_3_2.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)
BED$contig <- BED$V1
BED_T1 <- BED[BED$contig == 35 & BED$V2 >= 400736 | BED$contig == 82,]

#Append the orthogroup to each predicted gen
TARGET <- data.frame(gsub(" ","",Orthogroups$t1_10_3_2))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t1_10_3_2.,",",sep = ""))

#Initialize a vector to store the results
BED_T1$ORTH <- 0

#Run the function
for(i in 1:nrow(BED_T1)){
  GENE <- paste(",",BED_T1$ID[i],",",sep = "")
  BED_T1$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
#Plot Venn diagram comparison across genomes
x <- list(BED_CH14$ORTH, BED_CH434$ORTH, BED_T1$ORTH)

p1 <- ggVennDiagram(x, c("Hap. 1","Hap. 2","Hap. 3"),label_alpha=0.7)+
  scale_fill_gradient(low="gray95",high = "blue")+
  theme(legend.position = c(1, 0.8))+
          ggtitle("Chr. 5 left arm - control region")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill='Gene count')+scale_color_manual(values = c("black", "black", "black"))
ggsave("C:/Users/ericd/Downloads/Venn_control.png",plot = p1, width = 6, height = 5)


################################################################################
# Prepare diagrams for right-arm LSP region
################################################################################

#CH14
################################################################################
CONTIG <- 31
START <- 1534355
STOP <- 7111581

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t2_17_3_4i_13.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_CH14 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

#This function appends the orthogroup to each predicted gene
TARGET <- data.frame(gsub(" ","",Orthogroups$t2_17_3_4i_13))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t2_17_3_4i_13.,",",sep = ""))

BED_CH14$ORTH <- 0

for(i in 1:nrow(BED_CH14)){
  GENE <- paste(",",BED_CH14$ID[i],",",sep = "")
  BED_CH14$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

#CH434
################################################################################
CONTIG <- 12
START <- 537051
STOP <- 3078906

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/CH_434-inb3-a-1.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_CH434 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

#This function appends the orthogroup to each predicted gene
TARGET <- data.frame(gsub(" ","",Orthogroups$CH_434_inb3_a_1))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.CH_434_inb3_a_1.,",",sep = ""))

BED_CH434 <- BED_CH434[BED_CH434$ID != "g5074" & BED_CH434$ID != "g5075",]
BED_CH434$ORTH <- 0

for(i in 1:nrow(BED_CH434)){
  GENE <- paste(",",BED_CH434$ID[i],",",sep = "")
  BED_CH434$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

#T1
################################################################################
CONTIG <- 23
START <- 627218
STOP <- 3693875

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t1_10_3_2.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_T1 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

#This function appends the orthogroup to each predicted gene
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$t1_10_3_2))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t1_10_3_2.,",",sep = ""))

BED_T1 <- BED_T1[BED_T1$ID != "g4210",]
BED_T1$ORTH <- 0

for(i in 1:nrow(BED_T1)){
  GENE <- paste(",",BED_T1$ID[i],",",sep = "")
  BED_T1$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
#Comparison across genomes
x <- list(BED_CH14$ORTH, BED_CH434$ORTH, BED_T1$ORTH)


p2 <- ggVennDiagram(x, c("Hap. 1","Hap. 2","Hap. 3"),label_alpha=0.7)+
  scale_fill_gradient(low="gray95",high = "blue",limits = c(0, 400))+
  theme(legend.position = c(1, 0.8))+
  ggtitle("Chr. 5 right arm - HDH region")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill='Gene count')+scale_color_manual(values = c("black", "black", "black"))

ggsave("C:/Users/ericd/Downloads/Venn_LSP.png",plot = p2, width = 6, height = 5)

################################################################################
# Both sub-plots in a single image
################################################################################
p <- grid.arrange(p1, p2, ncol = 2) # Equivalent to nrow = 1
ggsave("C:/Users/ericd/Downloads/Venn_LSP.png",plot = p, width = 9, height = 4.5)

################################################################################
# Find private orthologs in Haplotype 1

privateBED <- BED_CH14[!(BED_CH14$ORTH %in% BED_CH434$ORTH) & !(BED_CH14$ORTH %in% BED_T1$ORTH),]

geneList <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/tsv_files/t2_17_3_4i_13_functional.tsv", header = FALSE)

geneList <- geneList[geneList$V1 %in% privateBED$ID,]

geneList <-merge(geneList, privateBED, by.x = "V1", by.y = "ID")

contig <- geneList$V1.y
startPos <- geneList$V2.y
stopPos <- geneList$V3.y
model <- geneList$V4.x
hitName <- geneList$V6.x
geneID <- geneList$V1

# All private orthologs
result <- data.frame(geneID,contig,startPos,stopPos,model,hitName)
write.csv(result, file = "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/CH14_private.csv")

# Private fucosyltransferase
result2 <- result[grepl("fuc",result$hitName, ignore.case = TRUE),]
write.csv(result2, file = "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/CH14_private_fucosyl.csv")

################################################################################
# Find private orthologs in Haplotype 2

privateBED <- BED_CH434[!(BED_CH434$ORTH %in% BED_CH14$ORTH) & !(BED_CH434$ORTH %in% BED_T1$ORTH),]

geneList <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/tsv_files/CH_434-inb3-a-1_functional.tsv", header = FALSE)

geneList <- geneList[geneList$V1 %in% privateBED$ID,]

geneList <-merge(geneList, privateBED, by.x = "V1", by.y = "ID")

contig <- geneList$V1.y
startPos <- geneList$V2.y
stopPos <- geneList$V3.y
model <- geneList$V4.x
hitName <- geneList$V6.x
geneID <- geneList$V1

# All private orthologs
result <- data.frame(geneID,contig,startPos,stopPos,model,hitName)
write.csv(result, file = "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/CH434_private.csv")

# Private fucosyltransferase
result2 <- result[grepl("fuc",result$hitName, ignore.case = TRUE),]

################################################################################
# Find private orthologs in Haplotype 3

privateBED <- BED_T1[!(BED_T1$ORTH %in% BED_CH14$ORTH) & !(BED_T1$ORTH %in% BED_CH434$ORTH),]

geneList <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/tsv_files/t3_12_3_1i_12_functional.tsv", header = FALSE)

geneList <- geneList[geneList$V1 %in% privateBED$ID,]

geneList <-merge(geneList, privateBED, by.x = "V1", by.y = "ID")

contig <- geneList$V1.y
startPos <- geneList$V2.y
stopPos <- geneList$V3.y
model <- geneList$V4.x
hitName <- geneList$V6.x
geneID <- geneList$V1

# All private orthologs
result <- data.frame(geneID,contig,startPos,stopPos,model,hitName)
write.csv(result, file = "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/T1_private.csv")

# Private fucosyltransferase
result2 <- result[grepl("fuc",result$hitName, ignore.case = TRUE),]