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

### Short read mapping to the contigs of interest

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





