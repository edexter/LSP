````
#Request interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

#Load BCFtools to fix missing genotype encoding problem
module purge
module load BCFtools

#Recode anything with less than DP 5 as missing (required since GATK update)
bcftools +setGT /scicore/home/ebertd/dexter0000/interlink/vcfsDaphnia_CH14/merged_daphnia_filtered.vcf -- -t q -n . -i 'FORMAT/DP<5' > CH14_filtered_fixed.vcf

#Load VCFtools for file conversion
module purge
module load VCFtools

#Create plink input files
vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig3 --chr "ptg000003l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig19 --chr "ptg000019l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig35 --chr "ptg000035l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig16 --chr "ptg000016l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30

vcftools --interchrom-geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14 --chr "ptg000016l_1" --chr "ptg000003l_1" --chr "ptg000035l_1" --chr "ptg000019l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30
````

````
##############
library(ggplot2)

contig3length <-  8776528
contig16length <- 1929634
contig19length <- 7624329
contig35length <-  2306407

df1 <- read.delim("C:/Users/ericd/Downloads/CH14_contig19.geno.ld")
maxPos <- contig19length
df1$POS1 <- maxPos - df1$POS1
df1$POS2 <- maxPos - df1$POS2

df2 <- read.delim("C:/Users/ericd/Downloads/CH14_contig3.geno.ld")
df3 <- read.delim("C:/Users/ericd/Downloads/CH14_contig16.geno.ld")

maxPos <- contig16length
df3$POS1 <- maxPos - df3$POS1
df3$POS2 <- maxPos - df3$POS2

df4 <- read.delim("C:/Users/ericd/Downloads/CH14_contig35.geno.ld")

#maxPos <- contig35length
#df4$POS1 <- maxPos - df4$POS1
#df4$POS2 <- maxPos - df4$POS2

df <- rbind(df1,df2,df3,df4)

df$POS1 <- ifelse(df$CHR=="ptg000019l_1", df$POS1,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000019l_1", df$POS2,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000016l_1", df$POS1+contig19length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000016l_1", df$POS2+contig19length,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000003l_1", df$POS1+contig19length+contig16length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000003l_1", df$POS2+contig19length+contig16length,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000035l_1", df$POS1+contig19length+contig16length+contig3length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000035l_1", df$POS2+contig19length+contig16length+contig3length,df$POS2)

##########################
df5 <- read.delim("C:/Users/ericd/Downloads/CH14.interchrom.geno.ld")

df5$POS1 <- ifelse(df5$CHR1=="ptg000019l_1", contig19length - df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000019l_1", contig19length - df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000016l_1", contig16length - df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000016l_1", contig16length - df5$POS2,df5$POS2)

#df5$POS1 <- ifelse(df5$CHR1=="ptg000035l_1", contig35length - df5$POS1,df5$POS1)
#df5$POS2 <- ifelse(df5$CHR2=="ptg000035l_1", contig35length - df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000019l_1", df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000019l_1", df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000035l_1", df5$POS1+contig19length+contig16length+contig3length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000035l_1", df5$POS2+contig19length+contig16length+contig3length,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000016l_1", df5$POS1+contig19length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000016l_1", df5$POS2+contig19length,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000003l_1", df5$POS1+contig19length+contig16length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000003l_1", df5$POS2+contig19length+contig16length,df5$POS2)

df5<-df5[,-3]
colnames(df5)[1]<-"CHR"

df<-rbind(df,df5)

#############
dfFlip<-df
dfFlip$POS2 <-df$POS1
dfFlip$POS1 <-df$POS2
df<-rbind(df,dfFlip)

colour_breaks <- c(0, 0.5, 0.9)
colours <- c("gray85", "yellow", "red")

df7 <- df[sample(1:nrow(df), 10000000), ]
p<-ggplot(df7[df7$N_INDV>1 ,], aes(x = rank(POS1), y = rank(POS2), color = R.2)) +
  #geom_point(size=0.001,alpha =0.5, shape = 0) +
  geom_tile(width=10000,height=10000,alpha = 0.5)+
  coord_fixed()+
  theme_bw()+
  ggtitle("CH14 Chromosome 5 LD")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Position in chromosome 5 (MB)") + ylab("Position in chromosome 5 (MB)")+
  scale_colour_gradientn(
    limits  = c(0,1),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = c(0,1)), 1),)
ggsave("C:/Users/ericd/Downloads/LD_chr5_rank.png",plot = p, width = 6, height = 6, dpi = 3000)

````



# Subset to just some samples

````bash
#Request interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

#Load VCFtools for file conversion
module purge
module load VCFtools

#Create plink input files
vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig3_A --chr "ptg000003l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30 --keep A_list.txt

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig19_A --chr "ptg000019l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30 --keep A_list.txt

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig35_A --chr "ptg000035l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30 --keep A_list.txt

vcftools --geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_contig16_A --chr "ptg000016l_1" --thin 2000 --min-r2 0.0 --maf 0.05 --max-missing 0.5 --minGQ 30 --keep A_list.txt

vcftools --interchrom-geno-r2 --vcf CH14_filtered_fixed.vcf --out CH14_A --chr "ptg000016l_1" --chr "ptg000003l_1" --chr "ptg000035l_1" --chr "ptg000019l_1" --thin 2000 --min-r2 0.2 --maf 0.05 --max-missing 0.5 --minGQ 30 --keep A_list.txt
````

````R
##############
library(ggplot2)

contig3length <-  8776528
contig16length <- 1929634
contig19length <- 7624329
contig35length <-  2306407

df1 <- read.delim("C:/Users/ericd/Downloads/CH14_contig19_A.geno.ld")
maxPos <- contig19length
df1$POS1 <- maxPos - df1$POS1
df1$POS2 <- maxPos - df1$POS2

df2 <- read.delim("C:/Users/ericd/Downloads/CH14_contig3_A.geno.ld")
df3 <- read.delim("C:/Users/ericd/Downloads/CH14_contig16_A.geno.ld")

maxPos <- contig16length
df3$POS1 <- maxPos - df3$POS1
df3$POS2 <- maxPos - df3$POS2

df4 <- read.delim("C:/Users/ericd/Downloads/CH14_contig35_A.geno.ld")

#maxPos <- contig35length
#df4$POS1 <- maxPos - df4$POS1
#df4$POS2 <- maxPos - df4$POS2

df <- rbind(df1,df2,df3,df4)

df$POS1 <- ifelse(df$CHR=="ptg000019l_1", df$POS1,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000019l_1", df$POS2,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000016l_1", df$POS1+contig19length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000016l_1", df$POS2+contig19length,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000003l_1", df$POS1+contig19length+contig16length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000003l_1", df$POS2+contig19length+contig16length,df$POS2)

df$POS1 <- ifelse(df$CHR=="ptg000035l_1", df$POS1+contig19length+contig16length+contig3length,df$POS1)
df$POS2 <- ifelse(df$CHR=="ptg000035l_1", df$POS2+contig19length+contig16length+contig3length,df$POS2)

##########################
df5 <- read.delim("C:/Users/ericd/Downloads/CH14.interchrom_A.geno.ld")

df5$POS1 <- ifelse(df5$CHR1=="ptg000019l_1", contig19length - df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000019l_1", contig19length - df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000016l_1", contig16length - df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000016l_1", contig16length - df5$POS2,df5$POS2)

#df5$POS1 <- ifelse(df5$CHR1=="ptg000035l_1", contig35length - df5$POS1,df5$POS1)
#df5$POS2 <- ifelse(df5$CHR2=="ptg000035l_1", contig35length - df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000019l_1", df5$POS1,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000019l_1", df5$POS2,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000035l_1", df5$POS1+contig19length+contig16length+contig3length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000035l_1", df5$POS2+contig19length+contig16length+contig3length,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000016l_1", df5$POS1+contig19length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000016l_1", df5$POS2+contig19length,df5$POS2)

df5$POS1 <- ifelse(df5$CHR1=="ptg000003l_1", df5$POS1+contig19length+contig16length,df5$POS1)
df5$POS2 <- ifelse(df5$CHR2=="ptg000003l_1", df5$POS2+contig19length+contig16length,df5$POS2)

df5<-df5[,-3]
colnames(df5)[1]<-"CHR"

df<-rbind(df,df5)

#############
dfFlip<-df
dfFlip$POS2 <-df$POS1
dfFlip$POS1 <-df$POS2
df<-rbind(df,dfFlip)

colour_breaks <- c(0, 0.5, 0.9)
colours <- c("gray85", "yellow", "red")

df7 <- df[sample(1:nrow(df), 10000000), ]
p<-ggplot(df7[df7$N_INDV>1 ,], aes(x = POS1, y = POS2, color = R.2)) +
  #geom_point(size=0.001,alpha =0.5, shape = 0) +
  geom_tile(width=10000,height=10000,alpha = 0.5)+
  coord_fixed()+
  theme_bw()+
  ggtitle("CH14 Chromosome 5 LD")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Position in chromosome 5 (MB)") + ylab("Position in chromosome 5 (MB)")+
  scale_colour_gradientn(
    limits  = c(0,1),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = c(0,1)), 1),)
ggsave("C:/Users/ericd/Downloads/LD_chr5_A.png",plot = p, width = 6, height = 6, dpi = 3000)
````

