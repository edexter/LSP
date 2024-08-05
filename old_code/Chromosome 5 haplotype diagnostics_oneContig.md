# Chromosome 5 diagnostics (Haplotype 1)

This project maps reads to three different references genomes, each containing a different variant forms of the chromosome 5 HDH (highly divergent haplotype region) - sometimes previously known as the LSP (large structural polymorphism). Missingness and read depth are then used to determine which versions of the HDH are present for many different D. magna clones as we have sequence data for. There are some idiosyncrasies in the formatting of each of the genome assemblies and so each one is processed with a separate script.



### Prepare reference genome and project folder

````bash
# Request an interactive job
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

# Create the project folder
mkdir bams_CH14_C3 bams_CH14_C3/stats_flag bams_CH14_C3/stats_cov bams_CH14_C3/stats_dup bams_CH14_C3/ref bams_CH14_C3/stats_LSP_depth

# Navigate to the main project folder
cd bams_CH14_C3

# Create a smaller reference genome with only contig of interest
module purge
module load SAMtools
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta "ptg000003l_1" > ref/CH14_C3.fasta

# Assign the new reference path as a variable
REF="ref/CH14_C3.fasta"

# Index the reference genome for BWA tool
module purge
module load bwa-mem2
bwa-mem2 index -p "$REF" "$REF"

# Index the reference genome for SAMtools
module purge 
module load SAMtools
samtools faidx "$REF"
````



### Short read mapping to the Haplotype 1 reference

````bash
#!/bin/bash

#SBATCH --job-name=bam_prep_daphnia_CH14		#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=4G              			#Memory reserved per core.
												#Total memory reserved: 32GB
#SBATCH --time=24:00:00	        				#Maximum time the job will run
#SBATCH --qos=1day	           					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/CH14_C3_map_out_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/CH14_C3_map_err_%A_%a

#Specifies an array of jobs from 1-8 with 100 max simultaneous
#SBATCH --array=1-258%100

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.
#Make sure reference is indexed before running (only needs to be done once)
#bwa index "$REF"

###############################################################################
#Define variables

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

REF="ref/CH14_C3.fasta"

################################################################################

#load required modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no
module load bwa-mem2
module load SAMtools
module load picard
################################################################################
#Part 1

cd /scicore/home/ebertd/dexter0000/interlink/bams_CH14_C3

#Map reads to reference genome, convert SAM output to BAM, and sort BAM.
#It's best to combine these steps to avoid writing large intermediate files
#to disk

echo "mapping reads of "$SAMP""
bwa-mem2 mem -t 8 -M "$REF" \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R1_trimmed.fq.gz \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R2_trimmed.fq.gz \
| samtools sort -@ 8 -o "$SAMP"_daphnia_sorted.bam

# print mapping statistics to screen
echo "Mapping stats for "$SAMP"_daphnia are:"
samtools flagstat "$SAMP"_daphnia_sorted.bam | tee stats_flag/"$SAMP"_daphnia_flagstat.txt

#Append read group metadata
echo "Adding readgroup to "$SAMP"_daphnia_sorted.bam"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I="$SAMP"_daphnia_sorted.bam \
O="$SAMP"_daphnia_R.bam \
RGID="$SAMP"  \
RGLB="$SAMP"  \
RGPL=illumina  \
RGPU="$SAMP"  \
RGSM="$SAMP"

#Index BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R.bam
	
#remove duplicates
echo "duplicates are being removed"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I="$SAMP"_daphnia_R.bam \
O="$SAMP"_daphnia_R_rdup.bam \
M=stats_dup/"$SAMP"_daphnia_duplicate_metrics.txt
	
#Index new BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_R_rdup.bam ]; then
	echo "removing intermediate files. ALMOST DONE!"
	rm "$SAMP"_daphnia_sorted.bam
	rm "$SAMP"_daphnia_R.bam
	rm "$SAMP"_daphnia_R.bam.bai
else 
	echo "Something went wrong! the sorted bam file has size 0" 
fi

################################################################################
#Part 2: Prepre BAM files for variant calling
################################################################################

#Remove unmapped reads from bam
echo "removing unmapped reads of "$SAMP" bam file"
samtools view -b -F 4 -@ 8 "$SAMP"_daphnia_R_rdup.bam \
> "$SAMP"_daphnia_Rm_rdup.bam

#Index new bam
echo "indexing new "$SAMP" bam file"
samtools index -@ 8 "$SAMP"_daphnia_Rm_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_Rm_rdup.bam ]; then
	echo "removing intermediate bam files."
	rm "$SAMP"_daphnia_R_rdup.bam
	rm "$SAMP"_daphnia_R_rdup.bam.bai
else
	echo "ERROR: the sorted bam file has size 0"
fi

#Get coverage stats per contig
samtools coverage "$SAMP"_daphnia_Rm_rdup.bam -o stats_cov/"$SAMP"_daphnia_coverage.txt
````



### Collect mapping depth stats from the HDH and flanking regions

This can be performed interactively.

````bash
#Request an interactive job
srun --nodes=1 --cpus-per-task=16 --mem=4G --pty bash

#Navigate to project folder
cd bams_CH14_C3

#Make an index of the bams
ls *.bam > bam.list

#Calculate depth per site per sample over a specific interval (MQ parameter is important). Multi-threading helps a lot here. This file is too large so we subsample every 100th base before writing to disk. This task takes 15-20 minutes to run. 
module purge
module load SAMtools
samtools depth -@ 16 -a -f bam.list -r "ptg000003l_1" --min-MQ 40 | awk 'NR % 100 == 0' > stats_LSP_depth/LSP_depth_CH14
````



### R code for plotting the results and further analysis

````R
#Load required packages
library(ggplot2)
library(slider)

#Import mapping depth
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_CH14", header=FALSE, comment.char="#")

#Import the sample names
ID <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/bam.list", quote="\"", comment.char="")
ID <- gsub("_.*", "", ID$V1)

#Populate the data frame with more useful sample names
colnames(df) <- c("Contig", "Position", ID)

#Calculate mean mapping depth per site across all samples and check distribution
siteMeans <- rowMeans(df[,-c(1:2)])
hist(siteMeans, breaks=100)

#Get rid of sites with excessive depth and recheck distribution
df <- df[siteMeans < 100,]
siteMeans <- rowMeans(df[,-c(1:2)])
hist(siteMeans, breaks=100)

#Calculate mean mapping depth per sample
sampleMeans <- colMeans(df[,-c(1:2)])
hist(sampleMeans,breaks = 100)

#Calculate depth over a sliding window
#df <- df[df$Depth < 50,]
Mean_slide <- slide_index_mean(x = siteMeans, i = df$Position, before = 100000, after = 100000, na_rm = TRUE)

#Label each position as "LSP", "L_flank", or "R_flank" region
region <-ifelse(df$Position < 2000000, "L_flank",
                   ifelse(df$Position > 7000000, "R_flank","LSP"))

#Get mean sample depth per region
sampleMeansLSP <- colMeans(df[region == "LSP",-c(1:2)])
sampleMeansL_flank <- colMeans(df[region == "L_flank",-c(1:2)])
sampleMeansR_flank <- colMeans(df[region == "R_flank",-c(1:2)])
df2 <- data.frame(cbind(sampleMeans,sampleMeansL_flank,sampleMeansLSP,sampleMeansR_flank))
df2$LSPratio <- df2$sampleMeansLSP / ((sampleMeansR_flank+sampleMeansL_flank)/2)
plot(df2$sampleMeansLSP ~ df2$sampleMeansL_flank)
plot(df2$sampleMeansLSP ~ df2$sampleMeansR_flank)
plot(df2$sampleMeansL_flank ~ df2$sampleMeansR_flank)
plot(df2$LSPratio)
hist(df2$LSPratio, breaks = 200)

#Plots
p <- ggplot(df[], aes(y=Mean_slide, x=Position))+ 
  theme_classic() +
  theme(legend.position = "none")+
  #annotate("rect", xmin=41759, xmax=43666, ymin=20, ymax=Inf, alpha=0.2, fill="blue")+ #repeats 
  #annotate("rect", xmin=1338469, xmax=1869000, ymin=20, ymax=Inf, alpha=0.2, fill="blue")+ #repeats
  #annotate("rect", xmin=64297, xmax=77384, ymin=0, ymax=25, alpha=0.2, fill="red")+ #repeats
  annotate("rect", xmin=7100000, xmax=7500000, ymin=0, ymax=Inf, alpha=0.2, fill="red")+ #repeats
  xlab("Genomic position") +ylab("Mean depth")+
  ggtitle("Contig 25 coverage (10 Kb sliding window)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point()+
  geom_hline(yintercept = c(25,35), lty = 1)
ggsave("C:/Users/ericd/Downloads/contig25.png",plot = p, width = 10, height = 6)

#Produce line plot of mapping depths showing position of LSP
df2 <- df[,3:length(ID)] / sampMeans[1:length(ID)]


#Calculate mean mapping depths per region per sample

#Determine if mean mapping depths are different between regions

#Determine mapping ratio per sample relative to the flanking regions

#Label all samples with CH14 like LSPs as having version "A"
df2$LSP <- ifelse(df2$LSPratio> 0.8, "A", "Unknown")

#Label Homozygotes for "A", Heterozygotes, and homozygotes for not "A"
df2$state <- ifelse(df2$LSPratio >= 1.1, "Homozygote A",
                    ifelse(df2$LSPratio < 1.1 & df2$LSPratio > 0.7, "Heterozgyote A", "Homozygote X"))

#Add sample IDs
df2 <- cbind(ID,df2)

#Export for IGV
write.table(df2, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_Ch14_igv.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
````



# Mapping reads to just one contig for CH434

Reference and folder preparation

````
#Request an interactive job
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

#Create the project folder
mkdir bams_CH434_C12 bams_CH434_C12/stats_flag bams_CH434_C12/stats_cov bams_CH434_C12/stats_dup bams_CH434_C12/ref bams_CH434_C12/stats_LSP_depth

#Navigate to the main project folder
cd bams_CH434_C12

#Create a smaller reference genome with only contig of interest
module purge
module load SAMtools
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/alt2/CH_434-inb3-a-1.falcon.polish.mask.fasta "000012F" > ref/CH434_C12.fasta

#Assign the new reference path as a variable
REF="ref/CH434_C12.fasta"

#Index the reference genome for BWA tool
module purge
module load bwa-mem2
bwa-mem2 index -p "$REF" "$REF"

#Index the reference genome for SAMtools
module purge 
module load SAMtools
samtools faidx "$REF"
````

Mapping to just LSP contig

````
#!/bin/bash

#SBATCH --job-name=bam_prep_daphnia_CH434		#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=4G              			#Memory reserved per core.
												#Total memory reserved: 32GB
#SBATCH --time=24:00:00	        			#Maximum time the job will run
#SBATCH --qos=1day	           					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/CH434_C3_map_out_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/CH434_C3_map_err_%A_%a

#Specifies an array of jobs from 1-8 with 100 max simultaneous
#SBATCH --array=1-258%100

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.
#Make sure reference is indexed before running (only needs to be done once)
#bwa index "$REF"

###############################################################################
#Define variables

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

REF="ref/CH434_C12.fasta"

################################################################################

#load required modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no
module load bwa-mem2
module load SAMtools
module load picard
################################################################################
#Part 1

cd /scicore/home/ebertd/dexter0000/interlink/bams_CH434_C12

#Map reads to reference genome, convert SAM output to BAM, and sort BAM.
#It's best to combine these steps to avoid writing large intermediate files
#to disk

echo "mapping reads of "$SAMP""
bwa-mem2 mem -t 8 -M "$REF" \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R1_trimmed.fq.gz \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R2_trimmed.fq.gz \
| samtools sort -@ 8 -o "$SAMP"_daphnia_sorted.bam

# print mapping statistics to screen
echo "Mapping stats for "$SAMP"_daphnia are:"
samtools flagstat "$SAMP"_daphnia_sorted.bam | tee stats_flag/"$SAMP"_daphnia_flagstat.txt

#Append read group metadata
echo "Adding readgroup to "$SAMP"_daphnia_sorted.bam"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I="$SAMP"_daphnia_sorted.bam \
O="$SAMP"_daphnia_R.bam \
RGID="$SAMP"  \
RGLB="$SAMP"  \
RGPL=illumina  \
RGPU="$SAMP"  \
RGSM="$SAMP"

#Index BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R.bam
	
#remove duplicates
echo "duplicates are being removed"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I="$SAMP"_daphnia_R.bam \
O="$SAMP"_daphnia_R_rdup.bam \
M=stats_dup/"$SAMP"_daphnia_duplicate_metrics.txt
	
#Index new BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_R_rdup.bam ]; then
	echo "removing intermediate files. ALMOST DONE!"
	rm "$SAMP"_daphnia_sorted.bam
	rm "$SAMP"_daphnia_R.bam
	rm "$SAMP"_daphnia_R.bam.bai
else 
	echo "Something went wrong! the sorted bam file has size 0" 
fi

################################################################################
#Part 2: Prepre BAM files for variant calling
################################################################################

#Remove unmapped reads from bam
echo "removing unmapped reads of "$SAMP" bam file"
samtools view -b -F 4 -@ 8 "$SAMP"_daphnia_R_rdup.bam \
> "$SAMP"_daphnia_Rm_rdup.bam

#Index new bam
echo "indexing new "$SAMP" bam file"
samtools index -@ 8 "$SAMP"_daphnia_Rm_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_Rm_rdup.bam ]; then
	echo "removing intermediate bam files."
	rm "$SAMP"_daphnia_R_rdup.bam
	rm "$SAMP"_daphnia_R_rdup.bam.bai
else
	echo "ERROR: the sorted bam file has size 0"
fi

#Get coverage stats per contig
samtools coverage "$SAMP"_daphnia_Rm_rdup.bam -o stats_cov/"$SAMP"_daphnia_coverage.txt
````

````
#Request an interactive job
srun --nodes=1 --cpus-per-task=16 --mem=4G --pty bash

#Navigate to project folder
cd bams_CH434_C12

#Make an index of the bams
ls *.bam > bam.list

#Calculate depth per site per sample over a specific interval (MQ parameter is important). Multi-threading helps a lot here. This file is too large so we subsample every 100th base before writing to disk. This task takes 15-20 minutes to run. 
module purge
module load SAMtools
samtools depth -@ 16 -a -f bam.list -r "000012F" --min-MQ 40 | awk 'NR % 100 == 0' > stats_LSP_depth/LSP_depth_CH434
````

R-code for determining LSP state relative to CH434

````
#Load required packages
library(ggplot2)
library(slider)

#Import mapping depth
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_CH434", header=FALSE, comment.char="#")

#Import the sample names
ID <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/bam.list", quote="\"", comment.char="")
ID <- gsub("_.*", "", ID$V1)

#Populate the data frame with more useful sample names
colnames(df) <- c("Contig", "Position", ID)

#Calculate mean mapping depth per site across all samples and check distribution
siteMeans <- rowMeans(df[,-c(1:2)])
hist(siteMeans, breaks=100)

#Get rid of sites with excessive depth and recheck distribution
df <- df[siteMeans < 100,]
siteMeans <- rowMeans(df[,-c(1:2)])
hist(siteMeans, breaks=100)

#Calculate mean mapping depth per sample
sampleMeans <- colMeans(df[,-c(1:2)])
hist(sampleMeans,breaks = 100)

#Calculate depth over a sliding window
#df <- df[df$Depth < 50,]
Mean_slide <- slide_index_mean(x = siteMeans, i = df$Position, before = 100000, after = 100000, na_rm = TRUE)

#Label each position as "LSP", "L_flank", or "R_flank" region
region <-ifelse(df$Position < 1800000, "L_flank",
                ifelse(df$Position > 3700000, "R_flank","LSP"))

#Get mean sample depth per region
sampleMeansLSP <- colMeans(df[region == "LSP",-c(1:2)])
sampleMeansL_flank <- colMeans(df[region == "L_flank",-c(1:2)])
sampleMeansR_flank <- colMeans(df[region == "R_flank",-c(1:2)])
df2 <- data.frame(cbind(sampleMeans,sampleMeansL_flank,sampleMeansLSP,sampleMeansR_flank))
df2$LSPratio <- df2$sampleMeansLSP / ((sampleMeansR_flank+sampleMeansL_flank)/2)
plot(df2$sampleMeansLSP ~ df2$sampleMeansL_flank)
plot(df2$sampleMeansLSP ~ df2$sampleMeansR_flank)
plot(df2$sampleMeansL_flank ~ df2$sampleMeansR_flank)
plot(df2$LSPratio)
hist(df2$LSPratio, breaks = 200)+abline(v=c(1.09,1.19))

#Plots
p <- ggplot(df[], aes(y=Mean_slide, x=Position))+ 
  theme_classic() +
  theme(legend.position = "none")+
  #annotate("rect", xmin=41759, xmax=43666, ymin=20, ymax=Inf, alpha=0.2, fill="blue")+ #repeats 
  #annotate("rect", xmin=1338469, xmax=1869000, ymin=20, ymax=Inf, alpha=0.2, fill="blue")+ #repeats
  #annotate("rect", xmin=64297, xmax=77384, ymin=0, ymax=25, alpha=0.2, fill="red")+ #repeats
  annotate("rect", xmin=1800000, xmax=3700000, ymin=0, ymax=Inf, alpha=0.2, fill="red")+ #repeats
  xlab("Genomic position") +ylab("Mean depth")+
  ggtitle("Contig 25 coverage (10 Kb sliding window)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point()+
  geom_hline(yintercept = c(25,35), lty = 1)
ggsave("C:/Users/ericd/Downloads/CH434_c12.png",plot = p, width = 10, height = 6)

#Produce line plot of mapping depths showing position of LSP
df2 <- df[,3:length(ID)] / sampMeans[1:length(ID)]


#Calculate mean mapping depths per region per sample

#Determine if mean mapping depths are different between regions

#Determine mapping ratio per sample relative to the flanking regions

#Label all samples with CH14 like LSPs as having version "A"
df2$LSP <- ifelse(df2$LSPratio> 1.09, "B", "Unknown")

#Label Homozygotes for "A", Heterozygotes, and homozygotes for not "A"
df2$state <- ifelse(df2$LSPratio >= 1.19, "Homozygote B",
                    ifelse(df2$LSPratio < 1.19 & df2$LSPratio > 1.09, "Heterozgyote B", "Homozygote X"))

#Add sample IDs
df2 <- cbind(ID,df2)


#Import previous haplotype results
CH14 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_Ch14_igv.txt")

CH14$LSP434 <- df2$LSP
CH14$state434 <- df2$state

#Remove sample with too low coverage
CH14 <- CH14[CH14$sampleMeans >= 5,]
plot(CH14$sampleMeans)
table(CH14$state,CH14$state434)

#Export for IGV
write.table(CH14, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_Ch14_igv.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
````



# Determine haplotypes relative to T1 (t1_10.3_2) clone

Reference and folder preparation

````
#Request an interactive job
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

#Create the project folder
mkdir bams_T1_C23 bams_T1_C23/stats_flag bams_T1_C23/stats_cov bams_T1_C23/stats_dup bams_T1_C23/ref bams_T1_C23/stats_LSP_depth

#Navigate to the main project folder
cd bams_T1_C23

#Create a smaller reference genome with only contig of interest
module purge
module load SAMtools
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000023l" > ref/T1_C23.fasta

#Assign the new reference path as a variable
REF="ref/T1_C23.fasta"

#Index the reference genome for BWA tool
module purge
module load bwa-mem2
bwa-mem2 index -p "$REF" "$REF"

#Index the reference genome for SAMtools
module purge 
module load SAMtools
samtools faidx "$REF"
````

Mapping to one LSP contig

````
#!/bin/bash

#SBATCH --job-name=bam_prep_daphnia_T1			#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=4G              			#Memory reserved per core.
												#Total memory reserved: 32GB
#SBATCH --time=24:00:00	        			#Maximum time the job will run
#SBATCH --qos=1day	           					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/T1_C23_map_out_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/T1_C23_map_err_%A_%a

#Specifies an array of jobs from 1-8 with 100 max simultaneous
#SBATCH --array=1-258%100

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.
#Make sure reference is indexed before running (only needs to be done once)
#bwa index "$REF"

###############################################################################
#Define variables

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

REF="ref/T1_C23.fasta"

################################################################################

#load required modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no
module load bwa-mem2
module load SAMtools
module load picard
################################################################################
#Part 1

cd /scicore/home/ebertd/dexter0000/interlink/bams_T1_C23

#Map reads to reference genome, convert SAM output to BAM, and sort BAM.
#It's best to combine these steps to avoid writing large intermediate files
#to disk

echo "mapping reads of "$SAMP""
bwa-mem2 mem -t 8 -M "$REF" \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R1_trimmed.fq.gz \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R2_trimmed.fq.gz \
| samtools sort -@ 8 -o "$SAMP"_daphnia_sorted.bam

# print mapping statistics to screen
echo "Mapping stats for "$SAMP"_daphnia are:"
samtools flagstat "$SAMP"_daphnia_sorted.bam | tee stats_flag/"$SAMP"_daphnia_flagstat.txt

#Append read group metadata
echo "Adding readgroup to "$SAMP"_daphnia_sorted.bam"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I="$SAMP"_daphnia_sorted.bam \
O="$SAMP"_daphnia_R.bam \
RGID="$SAMP"  \
RGLB="$SAMP"  \
RGPL=illumina  \
RGPU="$SAMP"  \
RGSM="$SAMP"

#Index BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R.bam
	
#remove duplicates
echo "duplicates are being removed"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I="$SAMP"_daphnia_R.bam \
O="$SAMP"_daphnia_R_rdup.bam \
M=stats_dup/"$SAMP"_daphnia_duplicate_metrics.txt
	
#Index new BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_R_rdup.bam ]; then
	echo "removing intermediate files. ALMOST DONE!"
	rm "$SAMP"_daphnia_sorted.bam
	rm "$SAMP"_daphnia_R.bam
	rm "$SAMP"_daphnia_R.bam.bai
else 
	echo "Something went wrong! the sorted bam file has size 0" 
fi

################################################################################
#Part 2: Prepre BAM files for variant calling
################################################################################

#Remove unmapped reads from bam
echo "removing unmapped reads of "$SAMP" bam file"
samtools view -b -F 4 -@ 8 "$SAMP"_daphnia_R_rdup.bam \
> "$SAMP"_daphnia_Rm_rdup.bam

#Index new bam
echo "indexing new "$SAMP" bam file"
samtools index -@ 8 "$SAMP"_daphnia_Rm_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_Rm_rdup.bam ]; then
	echo "removing intermediate bam files."
	rm "$SAMP"_daphnia_R_rdup.bam
	rm "$SAMP"_daphnia_R_rdup.bam.bai
else
	echo "ERROR: the sorted bam file has size 0"
fi

#Get coverage stats per contig
samtools coverage "$SAMP"_daphnia_Rm_rdup.bam -o stats_cov/"$SAMP"_daphnia_coverage.txt
````

Calculate coverage

````
#Request an interactive job
srun --nodes=1 --cpus-per-task=16 --mem=4G --pty bash

#Navigate to project folder
cd bams_T1_C23

#Make an index of the bams
ls *.bam > bam.list

#Calculate depth per site per sample over a specific interval (MQ parameter is important). Multi-threading helps a lot here. This file is too large so we subsample every 100th base before writing to disk. This task takes 15-20 minutes to run. 
module purge
module load SAMtools
samtools depth -@ 16 -a -f bam.list -r "ptg000023l" --min-MQ 40 | awk 'NR % 100 == 0' > stats_LSP_depth/LSP_depth_T1
````

R-code for further analysis and plots

````
#Load required packages
library(ggplot2)
library(slider)

#Import mapping depth
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_T1", header=FALSE, comment.char="#")

#Import the sample names
ID <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/bam.list", quote="\"", comment.char="")
ID <- gsub("_.*", "", ID$V1)

#Populate the data frame with more useful sample names
colnames(df) <- c("Contig", "Position", ID)

#Calculate mean mapping depth per site across all samples and check distribution
siteMeans <- rowMeans(df[,-c(1:2)])
hist(siteMeans, breaks=100)

#Get rid of sites with excessive depth and recheck distribution
df <- df[siteMeans < 100,]
siteMeans <- rowMeans(df[,-c(1:2)])
hist(siteMeans, breaks=100)

#Calculate mean mapping depth per sample
sampleMeans <- colMeans(df[,-c(1:2)])
hist(sampleMeans,breaks = 100)

#Calculate depth over a sliding window
#df <- df[df$Depth < 50,]
Mean_slide <- slide_index_mean(x = siteMeans, i = df$Position, before = 100000, after = 100000, na_rm = TRUE)

#Label each position as "LSP", "L_flank", or "R_flank" region
region <-ifelse(df$Position < 1200000, "L_flank",
                ifelse(df$Position > 3500000, "R_flank","LSP"))

#Get mean sample depth per region
sampleMeansLSP <- colMeans(df[region == "LSP",-c(1:2)])
sampleMeansL_flank <- colMeans(df[region == "L_flank",-c(1:2)])
sampleMeansR_flank <- colMeans(df[region == "R_flank",-c(1:2)])
df2 <- data.frame(cbind(sampleMeans,sampleMeansL_flank,sampleMeansLSP,sampleMeansR_flank))
df2$LSPratio <- df2$sampleMeansLSP / ((sampleMeansR_flank+sampleMeansL_flank)/2)

#Exclude low coverage samples
df2 <- df2[df2$sampleMeans >= 5,]
plot(df2$sampleMeansLSP ~ df2$sampleMeansL_flank)
plot(df2$sampleMeansLSP ~ df2$sampleMeansR_flank)
plot(df2$sampleMeansL_flank ~ df2$sampleMeansR_flank)
plot(df2$LSPratio)
cutoff1 <- 0.65
cutoff2 <- 0.85
hist(df2$LSPratio, breaks = 200)+abline(v=c(cutoff1,cutoff2))

#Plots
p <- ggplot(df[], aes(y=Mean_slide, x=Position))+ 
  theme_classic() +
  theme(legend.position = "none")+
  #annotate("rect", xmin=41759, xmax=43666, ymin=20, ymax=Inf, alpha=0.2, fill="blue")+ #repeats 
  #annotate("rect", xmin=1338469, xmax=1869000, ymin=20, ymax=Inf, alpha=0.2, fill="blue")+ #repeats
  #annotate("rect", xmin=64297, xmax=77384, ymin=0, ymax=25, alpha=0.2, fill="red")+ #repeats
  annotate("rect", xmin=1800000, xmax=3700000, ymin=0, ymax=Inf, alpha=0.2, fill="red")+ #repeats
  xlab("Genomic position") +ylab("Mean depth")+
  ggtitle("Contig 25 coverage (10 Kb sliding window)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point()+
  geom_hline(yintercept = c(25,35), lty = 1)
ggsave("C:/Users/ericd/Downloads/CH434_c12.png",plot = p, width = 10, height = 6)

#Produce line plot of mapping depths showing position of LSP
df2 <- df[,3:length(ID)] / sampMeans[1:length(ID)]


#Calculate mean mapping depths per region per sample

#Determine if mean mapping depths are different between regions

#Determine mapping ratio per sample relative to the flanking regions

#Label all samples with CH14 like LSPs as having version "A"
df2$LSP <- ifelse(df2$LSPratio> cutoff1, "C", "Unknown")

#Label Homozygotes for "A", Heterozygotes, and homozygotes for not "A"
df2$state <- ifelse(df2$LSPratio >= cutoff2, "Homozygote C",
                    ifelse(df2$LSPratio < cutoff2 & df2$LSPratio > cutoff1, "Heterozgyote C", "Homozygote X"))

#Add sample IDs
df2$ID <- row.names(df2)


#Import previous haplotype results
CH14 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_Ch14_igv.txt")
CH14 <- CH14[CH14$sampleMeans >= 5,]
df <- merge(CH14, df2, by = "ID")

#Make life simpler
table(df$state.x,df$state434,df$state.y)

#Remove sample with too low coverage
CH14 <- CH14[CH14$sampleMeans >= 5,]
plot(CH14$sampleMeans)
table(CH14$state,CH14$state434)


#Clean up the dataframe
stateCH14 <- df$state.x
stateCH434 <- df$state434
stateT1 <- df$state.y
dosCH14 <- df$LSP.x
dosCH434 <- df$LSP434
dosT1 <- df$LSP.y
meanCov <- df$sampleMeans.x
ID <- df$ID
res <- data.frame(cbind(ID, meanCov,stateCH14, stateCH434, stateT1, dosCH14, dosCH434, dosT1))

#Create a final state variable
res$state <- "Unknown"

#Assign homozygote states and verify dataframe (passes)
res$state <- ifelse(res$stateCH14 == "Homozygote A","A",res$state)
res$state <- ifelse(res$stateCH434 == "Homozygote B","B",res$state)
res$state <- ifelse(res$stateT1 == "Homozygote C","C",res$state)

#Assign heterozygote states and verify dataframe (passes)
res$state <- ifelse(res$stateCH14 == "Heterozgyote A" & res$stateCH434 == "Heterozgyote B","AB",res$state)
res$state <- ifelse(res$stateCH14 == "Heterozgyote A" & res$stateT1 == "Heterozgyote C","AC",res$state)
res$state <- ifelse(res$stateCH434 == "Heterozgyote B" & res$stateT1 == "Heterozgyote C","BC",res$state)
res$state <- as.factor(res$state)
res$meanCov <- as.numeric(res$meanCov)
table(res$state)
boxplot(res$meanCov ~ res$state)

#Export for IGV
write.table(res, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/LSP_depth_stats/LSP_depth_Ch14_igv.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
````

