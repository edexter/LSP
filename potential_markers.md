# Potential LSP markers

The goal here is to find a potential size-polymorphic marker that is diagnostic of the LSP genotype using bioinformatic methods. The following criteria should apply:

* Variants length differ by at least 5 bp
* Variants are present in at least 95% of samples



### Prepare VCF files for indel visualization using IGV

````bash
#Request interactive node
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

#Navigate to project folder
cd LSP

#Load required module
module purge
module load BCFtools

#Reduce raw VCF to only indels >4 bp
#Note: Need to specify that indel can be 4bp longer OR shorter
bcftools view --types "indels" /scicore/home/ebertd/dexter0000/interlink/vcfsDaphnia/merged_daphnia_filtered.vcf \
-i "ILEN>4 | ILEN<-4" \
-o scratch/merged_daphnia_filtered_indels.vcf

#Apply a depth filter to remove obviously spurious variants
bcftools view scratch/merged_daphnia_filtered_indels.vcf \
-i "INFO/DP<10000" \
-o scratch/merged_daphnia_filtered_indels_depth.vcf

#Convert to PLINK format for more convenient processing
#Note: Need to use new PLINK2 PGEN file format to handle multi-allelic variants
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --vcf scratch/merged_daphnia_filtered_indels_depth.vcf \
--make-pgen --out scratch/merged_daphnia_filtered_indels --allow-extra-chr

#Step 2: Calculate and filter missingness for retained variants
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --pfile scratch/merged_daphnia_filtered_indels --allow-extra-chr --geno 0.05 --make-pgen --out scratch/merged_daphnia_filtered_indels_missing

#Step 2: Calculate and filter maf (sum of non-major alleles)
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --pfile scratch/merged_daphnia_filtered_indels_missing --allow-extra-chr --maf 0.10 "minor" --make-pgen --out scratch/merged_daphnia_filtered_indels_missing_maf

#Export vcf
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --pfile scratch/merged_daphnia_filtered_indels_missing_maf --export vcf --allow-extra-chr --out scratch/merged_daphnia_filtered_indels_missing_maf --min-alleles 3
````



## Candidate markers based on visual search in IGV

````
#Possible tri-allelic single variant markers
Chr: 000027F|quiver [1,4,8]
Position: 425845
Reference: A*
Alternate: AAAAAACG,AACG

Chr: 000027F|quiver [1,3,8]
Position: 1477569
Reference: T*
Alternate: TTGAAAAA,TTG

Chr: 000025F|quiver [1,6,11] (K2 locus)
Position: 1606630
Reference: C*
Alternate: CTTTTT,CTTCTCTTTTT
````



#### Get genotypes at candidate positions

````
#Request interactive node
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

#Navigate to project folder
cd LSP

#Load required module
module purge
module load GATK
  
gatk IndexFeatureFile \
-I scratch/merged_daphnia_filtered_indels_missing_maf.vcf

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
 
 less scratch/LSP_marker_1_genotypes.table
  
````

### Make table of genotype LSP correspondance

````
################################################################################
#Load data
################################################################################

##MARKER 1
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_1_genotypes.txt")

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
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_2_genotypes.txt")

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
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_3_genotypes.txt")

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
markers <- merge(markers, M3, by = "ID")

##LSP DEPTH
#Load LSP depth data
LSP <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_Ch14_igv.txt")

markers <- merge(markers, LSP, by = "ID")

temp <- markers[markers$M1Q > 1 & markers$state != "Unknown",]
table(temp$M1, temp$state) #6 mismatches

temp <- markers[markers$M2Q > 1 & markers$state != "Unknown",]
table(temp$M2, temp$state) #10 mismatches

temp <- markers[markers$M3Q > 1 & markers$state != "Unknown",]
table(temp$M3, temp$state) #6 mismatches

#P20 RESISTANCE
#Load P20 resistance data
P20 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/IGVgroupsElocusAll.txt", header=TRUE)
P20 <- P20[,c(1,13)]
colnames(P20)[1] <- "ID"

markers <- merge(markers, P20, by = "ID")

temp <- markers[markers$phenotype == "RRS" | markers$phenotype == "RRR",]
table(temp$M1,temp$phenotype) # 6 mismatches
table(temp$M2,temp$phenotype) # 6 mismatches
table(temp$M3,temp$phenotype) # 5 mismatches
table(temp$state,temp$phenotype) #5 mismatches

#Pasteuria lineages
past <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/PastDAPCgroups.txt")
markers <- merge(markers, past, by = "ID")

temp <- markers[markers$phenotype == "RRS" | markers$phenotype == "RRR",]
table(temp$group, temp$state)
table(temp$group, temp$M3)
````





## Get unfiltered VCF for regions around candidate markers

We need to look for potential PCR primer binding sites

````
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash
module load BCFtools
cd LSP

#The VCF needs to be compressed and indexedfor further BCFtools actions
bgzip -c /scicore/home/ebertd/dexter0000/interlink/vcfsDaphnia/merged_daphnia_filtered.vcf > scratch/filtered.vcf.gz

bcftools index scratch/filtered.vcf.gz

bcftools view scratch/filtered.vcf.gz --regions "000027F|quiver:423345-428345" --output-type v > scratch/LSP_marker_candidate_1.vcf

bcftools view scratch/filtered.vcf.gz --regions "000027F|quiver:1475069-1480069" --output-type v > scratch/LSP_marker_candidate_2.vcf

bcftools view scratch/filtered.vcf.gz --regions "000025F|quiver:1604130-1609130" --output-type v > scratch/LSP_marker_candidate_3.vcf
````

​      

````
#Marker candidate 1
Wide-window 425820 < marker > 425855 (35 total). Can go up to 425730 < window > 425900 (170 total)

#Marker candidate 2
Wide-window 1477555 < marker > 1477580 (25 total). Can go up to 1477457 < marker > 1477662 (205 total).

#Marker candidate 3
Wide-window 1606596 < marker > 1606649 (53 total). Can go up to 1606477 < marker > 1606800 (323 total).
````

# Get predicted fragment lengths

````
#Iterate fragment length function across this index
module load SAMtools

cp /scicore/home/ebertd/dexter0000/interlink/scripts/index scratch/index

REF="/scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta"

#Marker 1
cat scratch/index | while read line 
do
   samtools faidx $REF "000027F|quiver:425730-425900"  | bcftools consensus --sample $line scratch/filtered.vcf.gz | tail -n +2 | tr -d '\n' | wc -c
done > scratch/fragment_lengths.txt

#Append sample names to fragment list
paste -d' ' scratch/index scratch/fragment_lengths.txt > scratch/fragement_lengths_LSP_marker_1.txt

#Marker 2
cat scratch/index | while read line 
do
   samtools faidx $REF "000027F|quiver:1477457-1477662"  | bcftools consensus --sample $line scratch/filtered.vcf.gz | tail -n +2 | tr -d '\n' | wc -c
done > scratch/fragment_lengths.txt

#Append sample names to fragment list
paste -d' ' scratch/index scratch/fragment_lengths.txt > scratch/fragement_lengths_LSP_marker_2.txt

#Marker 3
cat scratch/index | while read line 
do
   samtools faidx $REF "000025F|quiver:1606477-1606800"  | bcftools consensus --sample $line scratch/filtered.vcf.gz | tail -n +2 | tr -d '\n' | wc -c
done > scratch/fragment_lengths.txt

#Append sample names to fragment list
paste -d' ' scratch/index scratch/fragment_lengths.txt > scratch/fragement_lengths_LSP_marker_3.txt
````

#R code for plotting

````
#Load required packages
library(ggplot2)

#Load predicted LSP data
LSP <- LSP_depth_Ch14_igv <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_Ch14_igv.txt")

#Load and format predicted marker length data
lengths <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/fragement_lengths_LSP_marker_3.txt", quote="\"", comment.char="")
colnames(lengths) <- c("ID","length")
table(lengths$length)

df <- merge(LSP,lengths, by="ID")
df <- df[df$state == "A" | df$state == "B" | df$state == "C",]
df <- df[df$length != 0,]
df$state <- as.factor(df$state)
table(df$length,df$state)

#Add more data

#Nice plot showing marker 3
set.seed(3) #For nice jittering
ggplot(data=df, ggplot(data=df, ggplot(data=df, aes(x=state,y=length))+
  geom_jitter(height = 0,width = 0.30, aes(color=state))+
  theme_bw()+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("Predicted LSP genotype")+ylab("Predicted fragment length")+ggtitle("Expected fragment lengths for LSP marker (homozygotes)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(name = "Haplotype", labels = c("LSP1", "LSP2", "LSP3"))+
  scale_x_discrete(labels = c("LSP1","LSP2","LSP3"))+
  annotate("text", x=c(1,2,3), y=c(325,335,340), label = c("n = 10","n = 42","n = 41"))+
  annotate("text", x=1.2, y=319, label = c("*likely incorrect*"))

#Load resistotype data
Eloc <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/IGVgroupsElocusAll.txt", header=TRUE)
Eloc <- Eloc[,c(1,13)]
colnames(Eloc) <- c("ID","Eloc")

df2 <- merge(LSP,Eloc)
df2 <- df2[df2$Eloc=="RRR" | df2$Eloc == "RRS",]
table(df2$Eloc,df2$state)

#Load DAPC data
clusters <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups3.txt")
colnames(clusters)[1] <- "ID"
df3 <- merge(LSP,clusters, by = "ID")
df3 <- merge(df3,df2, by = "ID")
table(df3$E.Eclus,df3$state)
table(df3$Eloc,df3$E.Eclus)

#Load month data
month <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/CH14/vcf/IGV_resistotypes.txt")
colnames(month)[1] <-"ID"
df4 <- merge(LSP,month,by = "ID")
plot(df4$month~as.factor(df4$stateT1))
table(df4$month,df4$stateT1)

#Load lineage data
past <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/PastDAPCgroups.txt")
df4 <- merge(LSP,past,by = "ID")

table(df4$group,df4$state)

df4 <- merge(df4,resistance_loci_clusters_IGV, by = "ID")
table(df4$group,df4$stateCH14,df4$ABCbinary)
````



#### Get Nucleotide sequence for candidate marker

````
module load SAMtools
samtools faidx Xinb3_ref.fasta "000025F|quiver:1606130-1607130"
#Target variant is 1606630 (500 bp from start of sequence)
1606630

>000025F|quiver:1606130-1607130
TTTTACGAAGACCGAGAAAGAGACCTATTACAATCAATGGTGTCTTTGTAGTTTAAATAT
AAAACCTTTGGGTTTCTTTAAGAAAACAATTACTCACGAGGGTATAAGATCAACAGTGTC
AGACTGTCAGTGAACCATTCAGTCATTCACTTGTCAACCAACAGCAAAGCGTGTTAGATT
TGAGATAGAACTAAGAAGTTAGCACTAGCTTCTGAGTGTATCTGCACTTCTGTAGTATGC
ATTGTACGTCAAACTGACGCTGTGAAGCTTTGTAAGATTTAACGTTAACGTCACGTTGAC
GTTAACGTTATAATAACGCTATTAGACAACAGAGGGAACCACGCGCGAAAAAATAGTCAT
GTGACAAAATTATCTTTTTAAATATCTTTAAGAAAACCGATAAGATAAGATAAGATATCA
CCAGAATTTTGAAAGATAAAAGATATCAGAAGATACCAAAATTATCTTATCTTATCTTAT
CTTATCTTATCTTATCTTCTCGAAAGATGCAAAACACTGTTACCCCACCCACATCGTGCC
CCTCTATTCGGTTTAGGGTTGTCGAGAGCATGGATCAAAAGGTATTAATAGATCTGTAGG
CCTTACCCTCCTTGTGCAGATCCGAGAGAAAATCAAGGACAACGGTTAAAGAAGCAGACA
GGGGATCCTGATCCCGTCTAGCCTACCAAGTGAGCCATCGTTCCATGCGTCTGGTAGCAA
GCTGATGTGGATTATCGGTTTCCACCCATGAGAAGCTGGACCACTCGTTTCGAAAAGCCT
CTGGCTTTGAAGCGATCCCGGATAACGCCCAAGCGGCAACCAGGAGGGATCCTCGCGTCA
GAAGAGGGTGAGGGTTGCCCAGTGGGTACAAAAGCAGTTTCATGCTCGGTCGAATTACCT
GAGAAGGTTGGTGGGCTAGCTCCATAATTGTCAGCCACCGTGGTTGCGCTTGCCAAACCG
GTGCCACTAGAATCGGAGCGGCTTTTTCCTTCCTTATTTTT
````



