#Load required packages
library(readr)
library(stringr)
library(ggplot2)
library(ggVennDiagram)

################################################################################
#CH14
################################################################################
CONTIG <- 3
START <- 1535355
STOP <- 7111581

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t2_17_3_4i_13.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- str_sub(BED$V1, start = 1, end = -2)

BED_CH14 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

TSV <- read_delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/tsv_files/t2_17_3_4i_13_functional.tsv",
                   delim = "\t", escape_double = FALSE, 
                   col_names = FALSE, trim_ws = TRUE)

test <- TSV[TSV$X1 %in% BED_CH14$ID,]

#This function appends the orthogroup to each predicted gene
TARGET <- data.frame(gsub(" ","",Orthogroups$CH14_scaffold))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.CH14_scaffold.,",",sep = ""))

for(i in 1:nrow(BED_CH14)){
  GENE <- paste(",",BED_CH14$ID[i],",",sep = "")
  BED_CH14$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
#CHT3
################################################################################
CONTIG <- 23
START <- 0
STOP <- 10^9

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/t3_12_3_1i_12.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_T3 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

TSV <- read_delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/tsv_files/t3_12_3_1i_12_functional.tsv",
                  delim = "\t", escape_double = FALSE, 
                  col_names = FALSE, trim_ws = TRUE)

test <- TSV[TSV$X1 %in% BED_CH14$ID,]

#This function appends the orthogroup to each predicted gene
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$t3_12_3_1i_12))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.t3_12_3_1i_12.,",",sep = ""))

for(i in 1:nrow(BED_T3)){
  GENE <- paste(",",BED_T3$ID[i],",",sep = "")
  BED_T3$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

#CH434
CONTIG <- 12
START <- 1508248
STOP <- 4050114

BED <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/CH_434-inb3-a-1.bed", header=FALSE)
BED <- BED[BED$V8 == "gene",]
BED$ID <- gsub("ID=","",BED$V10)
BED$ID <- gsub(";","",BED$ID)

BED$contig <- BED$V1

BED_CH434 <- BED[BED$contig == CONTIG & BED$V2 >= START & BED$V3 <= STOP, ]

#This function appends the orthogroup to each predicted gene
#################################################
TARGET <- data.frame(gsub(" ","",Orthogroups$CH434_scaffold))
TARGET <- data.frame(paste(",",TARGET$gsub..........Orthogroups.CH434_scaffold.,",",sep = ""))

for(i in 1:nrow(BED_CH434)){
  GENE <- paste(",",BED_CH434$ID[i],",",sep = "")
  BED_CH434$ORTH[i] <- Orthogroups$Orthogroup[grep(GENE,unlist(TARGET))]
}

################################################################################
#Comparison across genomes
x <- list(BED_CH14$ORTH, BED_CH434$ORTH, BED_T3$ORTH)
ggVennDiagram(x, c("CH14", "CH434","T3"),
              fill_color = c("#0073C2FF", "#EFC000FF"))
str(x)
?ggVennDiagram
##########################################################################
LSP1 <- read_delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/tsv_files/t2_17_3_4i_13_functional.tsv", 
                   delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

BED1 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/bed_files/CH14_scaffold.bed", header=FALSE)

GENES1 <- BED1[BED1$V1 == 3,]

hist(BED1$V1,breaks = 3000,xlim = c(0,30))
LSP2 <- read_delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/CH434_functional.tsv",
                   delim = "\t", escape_double = FALSE, 
                   col_names = FALSE, trim_ws = TRUE)

#Note that the TSV has to be converted to a CSV for some cases
LSP3 <- read_csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/gene_prediction/T1_functional.csv", 
                 col_names = FALSE)

library(stringr)

# Rename the new column to "GO"

library(ggVennDiagram)

#Create vectors to plot in venn diagram
LSP1vec <- LSP1$X6
LSP1vec <- LSP1vec[LSP1vec != "-"]

LSP2vec <- LSP2$X6
LSP2vec <- LSP2vec[LSP2vec != "-"]

LSP3vec <- LSP3$X6
LSP3vec <- LSP3vec[LSP3vec != "-"]

# List of items

#x <- list(LSP2vec, LSP3vec)

# 3D Venn diagram
ggVennDiagram(x, c("LSP 1", "LSP 2", "LSP 3"),
              fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))

LSP1vec <- unlist(strsplit(LSP1$X14,split = "\\|"))
LSP1vec <- LSP1vec[LSP1vec != "-"]
LSP1vec <- LSP1vec[complete.cases(LSP1vec)]

LSP2vec <- unlist(strsplit(LSP2$X14,split = "\\|"))
LSP2vec <- LSP2vec[LSP2vec != "-"]
LSP2vec <- LSP2vec[complete.cases(LSP2vec)]

LSP3vec <- unlist(strsplit(LSP3$X14,split = "\\|"))
LSP3vec <- LSP3vec[LSP3vec != "-"]
LSP3vec <- LSP3vec[complete.cases(LSP3vec)]

x <- list(LSP1vec, LSP2vec, LSP3vec)
ggVennDiagram(x, c("LSP 1", "LSP 2", "LSP 3"),
              fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))

library(ggplot2)
barplot(sort(table(LSP1vec),decreasing=TRUE))
sort(table(LSP1vec))

barplot(sort(table(LSP2vec),decreasing=TRUE))
sort(table(LSP2vec))

#Get number of unique GO terms per LSP

#Get unique GO terms per LSP
LSP1uni <- unique(LSP1vec)
LSP2uni <- unique(LSP2vec)
LSP3uni <- unique(LSP3vec)

#Get length of private GO terms
length(LSP1uni[! LSP1uni %in% LSP2uni & ! LSP1uni %in% LSP3uni]) 
length(LSP2uni[! LSP2uni %in% LSP1uni & ! LSP2uni %in% LSP3uni]) 
length(LSP3uni[! LSP3uni %in% LSP1uni & ! LSP3uni %in% LSP2uni]) 

#Get shared GO terms
LSP1uni[LSP1uni %in% LSP2uni & LSP1uni %in% LSP3uni] 

length(unique(LSP3$X1))

x <- list(LSP1uni, LSP2uni,LSP3uni)

# 3D Venn diagram
ggVennDiagram(x, c("LSP 1", "LSP 2", "LSP 3"),
              fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))

##Same with gene descriptions
################################################################################
LSP1vec <- unlist(strsplit(LSP1$X13,split = "\\|"))
LSP1vec <- LSP1vec[LSP1vec != "-"]
LSP1vec <- LSP1vec[complete.cases(LSP1vec)]

LSP2vec <- unlist(strsplit(LSP2$X13,split = "\\|"))
LSP2vec <- LSP2vec[LSP2vec != "-"]
LSP2vec <- LSP2vec[complete.cases(LSP2vec)]

LSP3vec <- unlist(strsplit(LSP3$X13,split = "\\|"))
LSP3vec <- LSP3vec[LSP3vec != "-"]
LSP3vec <- LSP3vec[complete.cases(LSP3vec)]

x <- list(LSP1vec, LSP2vec, LSP3vec)
ggVennDiagram(x, c("LSP 1", "LSP 2", "LSP 3"),
              fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))

library(ggplot2)
barplot(sort(table(LSP1vec),decreasing=TRUE))
sort(table(LSP1vec))

barplot(sort(table(LSP2vec),decreasing=TRUE))
sort(table(LSP2vec))

#Get number of unique GO terms per LSP

#Get unique GO terms per LSP
LSP1uni <- unique(LSP1vec)
LSP2uni <- unique(LSP2vec)
LSP3uni <- unique(LSP3vec)

length(LSP1uni)
length(LSP2uni)
length(LSP3uni)

#Get length of private GO terms
length(LSP1uni[! LSP1uni %in% LSP2uni & ! LSP1uni %in% LSP3uni]) 
length(LSP2uni[! LSP2uni %in% LSP1uni & ! LSP2uni %in% LSP3uni]) 
length(LSP3uni[! LSP3uni %in% LSP1uni & ! LSP3uni %in% LSP2uni]) 

#Get shared GO terms
LSP1uni[LSP1uni %in% LSP2uni & LSP1uni %in% LSP3uni] 

length(unique(LSP3$X1))

x <- list(LSP1uni, LSP2uni,LSP3uni)

# 3D Venn diagram
ggVennDiagram(x, c("LSP 1", "LSP 2", "LSP 3"),
              fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))

LSP3uni

table(LSP1vec[! LSP1vec %in% LSP2vec & ! LSP1vec %in% LSP3vec]) 
table(LSP2vec[! LSP2vec %in% LSP1vec & ! LSP2vec %in% LSP3vec]) 
table(LSP3vec[! LSP3vec %in% LSP1vec & ! LSP3vec %in% LSP2vec]) 
