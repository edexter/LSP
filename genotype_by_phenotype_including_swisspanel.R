################################################################################
#Load data
################################################################################
#Corrected sample names list
newNames <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/samples_to_rename.txt", quote="\"")
colnames(newNames) <- c("old", "new")

################################################################################
##MARKER 1 (Interlink)
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_1_genotypes_interlink.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M1"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M1Q"

#Merge formatted genotype and quality data for marker
M1A <- merge(GT, GQ , by = "ID")

##MARKER 1 (Swisspanel)
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_1_genotypes_swisspanel.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M1"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M1Q"

#Merge formatted genotype and quality data for marker
M1B <- merge(GT, GQ , by = "ID")

M1 <- rbind(M1A,M1B)

#Remove duplicates
M1 <- M1[!duplicated(M1$ID),]

#Format names correctly
M1$ID <- ifelse(nchar(M1$ID) == 4,paste("CH-H-2019-",M1$ID, sep =""),M1$ID)

for(i in 1:nrow(newNames)){
  M1$ID[M1$ID == newNames[i,1]] <- newNames[i,2]
}

##MARKER 2
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_2_genotypes_interlink.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M2"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M2Q"

#Merge formatted genotype and quality data for marker
M2A <- merge(GT, GQ , by = "ID")

##MARKER 2 (Swisspanel)
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_2_genotypes_swisspanel.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M2"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M2Q"

#Merge formatted genotype and quality data for marker
M2B <- merge(GT, GQ , by = "ID")

M2 <- rbind(M2A,M2B)

#Remove duplicates
M2 <- M2[!duplicated(M2$ID),]

#Format names correctly
M2$ID <- ifelse(nchar(M2$ID) == 4,paste("CH-H-2019-",M2$ID, sep =""),M2$ID)

for(i in 1:nrow(newNames)){
  M2$ID[M2$ID == newNames[i,1]] <- newNames[i,2]
}

##MARKER 3
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_3_genotypes_interlink.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M3"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M3Q"

#Merge formatted genotype and quality data for marker
M3A <- merge(GT, GQ , by = "ID")

##MARKER 3 (Swisspanel)
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_3_genotypes_swisspanel.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M3"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M3Q"

#Merge formatted genotype and quality data for marker
M3B <- merge(GT, GQ , by = "ID")

M3 <- rbind(M3A,M3B)

#Remove duplicates
M3 <- M3[!duplicated(M3$ID),]

#Format names correctly
M3$ID <- ifelse(nchar(M3$ID) == 4,paste("CH-H-2019-",M3$ID, sep =""),M3$ID)

for(i in 1:nrow(newNames)){
  M3$ID[M3$ID == newNames[i,1]] <- newNames[i,2]
}

markers <- merge(M1, M2, by = "ID")
markers <- merge(markers, M3, by = "ID")

##LSP DEPTH
#Load LSP depth data
LSP <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/LSP_depth_genotypes.txt")
LSP$ID <- ifelse(nchar(LSP$ID) == 4,paste("CH-H-2019-",LSP$ID, sep =""),LSP$ID)
markers <- merge(markers, LSP, by = "ID")

temp <- markers[markers$M1Q > 1 & markers$state != "Unknown",]
table(temp$M1, temp$state) #6 mismatches

temp <- markers[markers$M2Q > 1 & markers$state != "Unknown",]
table(temp$M2, temp$state) #10 mismatches

temp <- markers[markers$M3Q > 1 & markers$state != "Unknown",]
table(temp$M3, temp$state) #6 mismatches

#P20 RESISTANCE
#Load and format phenotype data
#P20 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/IGVgroupsElocusAll.txt", header=TRUE)
P20 <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/phenotypes/resistotypes_binary.csv")
P20 <- P20[,c(1,2,3,12)]
P20$phenotype <- paste(P20$C1,P20$C19,P20$P20,sep="")
colnames(P20)[1] <- "ID"

markers <- merge(markers, P20, by = "ID")

temp <- markers[markers$phenotype == "RRS" | markers$phenotype == "RRR",]

table(temp$M1,temp$phenotype) # 6 mismatches
table(temp$M2,temp$phenotype) # 6 mismatches
table(temp$M3,temp$phenotype) # 5 mismatches
table(temp$state,temp$phenotype) #5 mismatches

#Get number of samples in table
length(temp$phenotype)
#Pasteuria lineages
past <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/PastDAPCgroups.txt")
markers <- merge(markers, past, by = "ID")

temp <- markers[markers$phenotype == "RRS" | markers$phenotype == "RRR",]
table(temp$group, temp$state)
table(temp$group, temp$M3)
