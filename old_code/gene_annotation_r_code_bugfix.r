#Load required packages
library(readr)
library(stringr)
library(ggplot2)
library(ggVennDiagram)

################################################################################
#CH14
################################################################################
CONTIG <- 191
START <- 3000000
STOP <- 7001000

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

################################################################################
#T3
################################################################################
CONTIG <- 33
START <- 0
STOP <- 10^9

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t3_12_3_1i_12.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_T3 <- BED[BED$contig == 33 & BED$V2 >= 2619543 | BED$contig == 87
              | BED$contig == 20 & BED$V2 <= 606860,]

#This function appends the orthogroup to each predicted gene
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$t3_12_3_1i_12))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t3_12_3_1i_12.,",",sep = ""))

for(i in 1:nrow(BED_T3)){
  GENE <- paste(",",BED_T3$ID[i],",",sep = "")
  BED_T3$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
#CH434
################################################################################
CONTIG <- 5
START <- 613227
STOP <- 4565149

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/CH_434-inb3-a-1.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_CH434 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

#This function appends the orthogroup to each predicted gene
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$CH_434_inb3_a_1))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.CH_434_inb3_a_1.,",",sep = ""))

BED_CH434 <- BED_CH434[BED_CH434$ID != "g5074" & BED_CH434$ID != "g5075",]
BED_CH434$ORTH <- 0

for(i in 1:nrow(BED_CH434)){
  GENE <- paste(",",BED_CH434$ID[i],",",sep = "")
  BED_CH434$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
#T1
################################################################################
CONTIG <- 5
START <- 613227
STOP <- 4565149

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t1_10_3_2.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_T1 <- BED[BED$contig == 35 & BED$V2 >= 400736 | BED$contig == 82,]

#This function appends the orthogroup to each predicted gene
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$t1_10_3_2))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t1_10_3_2.,",",sep = ""))

#BED_T1 <- BED_T1[BED_T1$ID != "g5074" & BED_CH434$ID != "g5075",]
BED_T1$ORTH <- 0

for(i in 1:nrow(BED_T1)){
  GENE <- paste(",",BED_T1$ID[i],",",sep = "")
  BED_T1$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################

#Comparison across genomes
x <- list(BED_CH14$ORTH, BED_CH434$ORTH, BED_T1$ORTH)

ggVennDiagram(x, c("CH14","CH434","T1"))

################################################################################
################################################################################
################################################################################

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

################################################################################
#T3
################################################################################
CONTIG <- 23
START <- 1533405
STOP <- 7110627

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t3_12_3_1i_12.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_T3 <- BED[BED$contig == CONTIG | BED$contig == 87 | BED$contig == 20 & BED$V2 >= START & BED$V3 <= STOP, ]

#This function appends the orthogroup to each predicted gene
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$t3_12_3_1i_12))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t3_12_3_1i_12.,",",sep = ""))

for(i in 1:nrow(BED_T3)){
  GENE <- paste(",",BED_T3$ID[i],",",sep = "")
  BED_T3$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
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
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$CH_434_inb3_a_1))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.CH_434_inb3_a_1.,",",sep = ""))

BED_CH434 <- BED_CH434[BED_CH434$ID != "g5074" & BED_CH434$ID != "g5075",]
BED_CH434$ORTH <- 0

for(i in 1:nrow(BED_CH434)){
  GENE <- paste(",",BED_CH434$ID[i],",",sep = "")
  BED_CH434$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
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


ggVennDiagram(x, c("CH14","CH434","T1"))

################################################################################
################################################################################