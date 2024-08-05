# Chromosome 5 haplotype diagnostics

This project maps reads to three different references genomes, each containing a different variant of the large structural polymorphism (LSP) on chromosome 5. Missingness and read depth are then used to determine which versions of the LSP are present in as many D. magna clones as we have sequence data for.

# Map reads to LSP 1 reference (CH14 genome) 

Reference and folder preparation. Can be run interactively.

````bash
#Request an interactive job
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

#Create the project folder
mkdir bams_CH14_C3 bams_CH14_C3/stats_flag bams_CH14_C3/stats_cov bams_CH14_C3/stats_dup bams_CH14_C3/ref bams_CH14_C3/stats_LSP_depth

#Navigate to the main project folder
cd bams_CH14_C3

#Create a smaller reference genome with only contig of interest
module purge
module load SAMtools
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta "ptg000003l_1" > ref/CH14_C3.fasta

#Assign the new reference path as a variable
REF="ref/CH14_C3.fasta"

#Index the reference genome for BWA tool
module purge
module load bwa-mem2
bwa-mem2 index -p "$REF" "$REF"

#Index the reference genome for SAMtools
module purge 
module load SAMtools
samtools faidx "$REF"
````



Mapping to LSP 1 reference. Should be submitted as a script.

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

#print mapping statistics to screen
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



Get coverage stats per sample for LSP and flanking regions. Performed interactively.

````
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



# Map reads to LSP 2 reference (CH434 genome) 

Reference and folder preparation. Can be run interactively.

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



Mapping to LSP 2 reference. Should be submitted as a script.

````
#!/bin/bash

#SBATCH --job-name=bam_prep_daphnia_CH434		#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=4G              			#Memory reserved per core.
												#Total memory reserved: 32GB
#SBATCH --time=24:00:00	        				#Maximum time the job will run
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



Get coverage stats per sample for LSP and flanking regions. Performed interactively.

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



# Map reads to LSP 3 reference (T1 genome) 

Reference and folder preparation. Can be run interactively.

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



Mapping to LSP 2 reference. Should be submitted as a script.

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



Get coverage stats per sample for LSP and flanking regions. Performed interactively.

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

