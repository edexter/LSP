````
#Request interactive node
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

#Navigate to project folder
cd LSP

#Load required module
module purge
module load GATK


#Interlink dataset

gatk IndexFeatureFile \
-I merged_daphnia_filtered.vcf.gz

#Candidate marker 1
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/interlink/vcfsDaphnia/merged_daphnia_filtered.vcf \
  --intervals "000027F|quiver:425845" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_1_genotypes.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta

#Candidate marker 2
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/interlink/vcfsDaphnia/merged_daphnia_filtered.vcf \
  --intervals "000027F|quiver:1477569" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_2_genotypes.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta
  
#Candidate marker 3
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/interlink/vcfsDaphnia/merged_daphnia_filtered.vcf \
  --intervals "000025F|quiver:1606630" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_3_genotypes.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta
 
 #Swisspond dataset
 
 gatk IndexFeatureFile \
-I /scicore/home/ebertd/dexter0000/EAWAG/VCF/daphnia_unfiltered.vcf.gz

 #Candidate marker 1
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/EAWAG/VCF/daphnia_unfiltered.vcf.gz \
  --intervals "000027F|quiver:425845" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_1_genotypes_swisspanel.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta

#Candidate marker 2
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/EAWAG/VCF/daphnia_unfiltered.vcf.gz \
  --intervals "000027F|quiver:1477569" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_2_genotypes_swisspanel.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta
  
#Candidate marker 3
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/EAWAG/VCF/daphnia_unfiltered.vcf.gz \
  --intervals "000025F|quiver:1606630" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_3_genotypes_swisspanel.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta
 
 #Diversity panel
 
 #Candidate marker 1
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/diversity_panel/dmdp_fabienne.snps.vcf.gz \
  --intervals "000027F|quiver:405845-425845" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_1_genotypes_divpanel.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta

#Candidate marker 2
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/diversity_panel/dmdp_fabienne.snps.vcf.gz \
  --intervals "000027F|quiver:1477569" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_2_genotypes_divpanel.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta
  
#Candidate marker 3
  gatk VariantsToTable \
  -V /scicore/home/ebertd/dexter0000/diversity_panel/dmdp_fabienne.snps.vcf.gz \
  --intervals "000025F|quiver:1606630" -F ID -F REF -F ALT -GF GT -GF GQ \
  -O scratch/LSP_marker_3_genotypes_divpanel.txt \
  -R /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta
 
````



# R code for genotype-phenotype analysis

````
#### Interlink ----

##MARKER 1
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
M1 <- merge(GT, GQ , by = "ID")

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
M2 <- merge(GT, GQ , by = "ID")

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
M3 <- merge(GT, GQ , by = "ID")

markers <- merge(M1, M2, by = "ID")
markers_interlink <- merge(markers, M3, by = "ID")

#### Swisspanel ----

##MARKER 1
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
M1 <- merge(GT, GQ , by = "ID")

##MARKER 2
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
M2 <- merge(GT, GQ , by = "ID")

##MARKER 3
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
M3 <- merge(GT, GQ , by = "ID")

markers <- merge(M1, M2, by = "ID")
markers_swisspanel <- merge(markers, M3, by = "ID")
#### Divpanel ----

##MARKER 1
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_1_genotypes_divpanel.txt")

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
M1 <- merge(GT, GQ , by = "ID")

##MARKER 2
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_2_genotypes_divpanel.txt")

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
M2 <- merge(GT, GQ , by = "ID")

##MARKER 3
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_3_genotypes_divpanel.txt")

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
M3 <- merge(GT, GQ , by = "ID")

markers <- merge(M1, M2, by = "ID")
markers_divpanel <- M3

markers <- rbind(markers_interlink[,c(1,6,7)],markers_swisspanel[c(1,6,7)],
                 markers_divpanel)

length(unique(markers$ID[markers$M3Q!="./."]))

#### Phenotypes ----
#Import stick test data
pheno <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/resistotypes_binary.csv")

#Fix ID formating
markers$ID <- gsub("\\.","-",markers$ID)
markers$ID <- ifelse(nchar(markers$ID) == 4, paste("CH-H-2019-",markers$ID,sep=""),markers$ID)
markers$ID[!(markers$ID %in% pheno$host)]

df <- merge(markers, pheno, by.x = "ID", by.y = "host")
df <- df[df$C1 == "R" & df$C19 == "R",]
df <- df[unique(df$ID)]
table(df$P20, df$M3)

temp <- df[df$P20 == "R", ]
````



