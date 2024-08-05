# Gene annotation and synteny plots with whole D. magna genomes

### Format_genomes.sh

First upload all genome assemblies to [genomes_original] directory and manually clean up the file names. This script then performs various formatting operations required for  the GENESPACE program. This script completes in a few minutes.

````
#!/bin/bash

#SBATCH --job-name=format_genomes		#Job name
#SBATCH --cpus-per-task=1	        	#Number of cores reserved
#SBATCH --mem-per-cpu=4G            	#Memory reserved per core.
										#Total memory reserved: 64GB (required!)
#SBATCH --time=24:00:00	        		#Maximum time the job will run
#SBATCH --qos=1day           			#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/LSP/logs/genome_format_%A_%a.out

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/LSP/logs/genome_format_%A_%a.err

#Specifies an array of jobs from 1-n with n max simultaneous
#SBATCH --array=1-20%20

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

########################################################################
#Define variables

#Interlink IDs list
INDEXFILE=/scicore/home/ebertd/dexter0000/LSP/scripts/genome_index.txt

GENOME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)
################################################################################

#Navigte to project directory
cd /scicore/home/ebertd/dexter0000/LSP

#Remove any previous genome masking by converting all nucleotides to uppercase
zcat genomes_original/"$GENOME".gz | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' | gzip > genomes/"$GENOME".gz

#Remove any non-numeric characters from the fasta headers and leading zeros
zcat genomes/"$GENOME".gz | sed -E '/^>/ { s/[^0-9>]//g; s/>0+/>/; }' | gzip > genomes/"$GENOME".2.gz

#Make sure all contigs are numbered
zcat genomes/"$GENOME".2.gz  | sed -e '/^>$/s/^>$/>0/' | gzip > genomes/"$GENOME".3.gz

#Clean up intermediate files and rename final output
#Note that some of the downstream applications will require an unzipped reference
rm genomes/"$GENOME".gz genomes/"$GENOME".2.gz
mv genomes/"$GENOME".3.gz genomes/"$GENOME".fa.gz
gunzip genomes/"$GENOME".fa.gz

if [ -s genomes/"$GENOME".fa ]
then 
        printf ""Finished formatting "$GENOME""
else 
        printf ""there was a problem formatting "$GENOME""
fi
````



## Annotate_genomes.sh

Predict genes using AUGUSTUS and predict gene function using INTERPROSCAN. This completes in a little under 24 hours when submitted as an array job.

````bash
#!/bin/bash

#SBATCH --job-name=annotate_genomes		#Job name
#SBATCH --cpus-per-task=16	        	#Number of cores reserved
#SBATCH --mem-per-cpu=4G            	#Memory reserved per core.
										#Total memory reserved: 64GB (required!)
#SBATCH --time=24:00:00	        		#Maximum time the job will run
#SBATCH --qos=1day           			#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/LSP/logs/gene_predict_%A_%a.out

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/LSP/logs/gene_predict_%A_%a.err

#Specifies an array of jobs from 1-n with n max simultaneous
#SBATCH --array=1-20%20

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

########################################################################
#Define variables

#Interlink IDs list
INDEXFILE=/scicore/home/ebertd/dexter0000/LSP/scripts/genome_index.txt

GENOME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)
################################################################################
#Navigate to project folder
cd /scicore/home/ebertd/dexter0000/LSP

#Load required module
module load RepeatMasker

#Repeat masking
#[-dir] specifies the output directory
#[-nolow] is specific to gene prediction downstream usage
#[-small] soft-masks repetitive regions instead of replacement with N's
#Note that you cannot use RepeatMasker with a gzipped assembly fasta file

echo Repeat masking genome "$GENOME"

RepeatMasker genomes/"$GENOME".fa -pa 16 -gff -xsmall -lib RM_120310.ThuApr201445402023/consensi.fa.classified -dir FASTA 

if [ -s FASTA/"$GENOME".fa.masked  ]
then 
        echo "Repeat masking genome "$GENOME" completed"
else 
        printf "Repeat masking genome "$GENOME" did not complete successfully"
fi

#Load required module
module purge
module load AUGUSTUS

#Required setup steps before initial run
#cp -r /scicore/soft/apps/AUGUSTUS/3.4.0-foss-2021a/config .

#Needs to be added at every session if not permanently added to bash.rc
export AUGUSTUS_CONFIG_PATH=/scicore/home/ebertd/dexter0000/LSP/config/

#Run Augustus
echo Starting augustus gene prediction for genome "$GENOME"

augustus --species=Dmagna_iso FASTA/"$GENOME".fa.masked > annotations/"$GENOME".gtf

gtf2gff.pl <annotations/"$GENOME".gtf --out=annotations/"$GENOME".gff3 --gff3

if [ -s annotations/"$GENOME".gff3 ]
then 
        echo "Augustus gene prediction for "$GENOME" completed"
else 
        printf "Augustus gene prediction for "$GENOME" did not complete successfully"
fi

#Extract protein sequences from annotation (note that the order of the flags matters!)

echo "Extracting protein sequences from annotation"

scripts/protein_extract.sh -p annotations/"$GENOME".protein.faa -c annotations/"$GENOME".CDS.faa annotations/"$GENOME".gtf

if [ -s annotations/"$GENOME".protein.faa ]
then 
        echo "Protein sequences extracted"
else 
        printf "The augustus_extract.sh script was not able to extract protein sequences"
fi

#Functional annotation with interproscan

#Swap module
module purge
module load InterProScan

#Run interproscan. The [-dp] flag disables pre-calculated match lookup service.
#Note that his requires a lot of memory and will crash with less than 64 GB

echo "Starting InterProScan functional annotation"

interproscan.sh -i annotations/"$GENOME".protein.faa -cpu 16 -dp -goterms --output-file-base annotations/"$GENOME"_functional

if [ -s annotations/"$GENOME"_functional.gff3 ]
then 
        echo "Completed InterProScan functional annotation. Gene Annotation complete."
else 
        printf "InterProScan functional annotation encountered an error"
fi
````

# Convert GFF to BED

````
srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash
module load BEDOPS

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/t2_17_3_4i_13.gff3 > t2_17_3_4i_13.bed

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/CH_434-inb3-a-1.gff3 > CH_434-inb3-a-1.bed

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/t1_10_3_2.gff3 > t1_10_3_2.bed

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/t3_12_3_1i_12.gff3 > t3_12_3_1i_12.bed
````

### Run genespace

Perform synteny analysis using GENESPACE. This script performs some preliminary setup for the GENESPACE run and the calls and R script to perform the analysis. An entire run takes XXX.

````bash
#!/bin/bash

#SBATCH --job-name=genespace			#Job name
#SBATCH --cpus-per-task=32	        	#Number of cores reserved
#SBATCH --mem-per-cpu=4G            	#Memory reserved per core.
										#Total memory reserved: 128GB (required!)
#SBATCH --time=168:00:00	        	#Maximum time the job will run
#SBATCH --qos=1week           			#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/LSP/logs/genespace.out

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/LSP/logs/genespace.err

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

#Activate conda environment
conda init bash
source ~/.bashrc #This line is mandatory
conda activate orthofinder

#Temporarily allow loading of multiple modules simultaneously (SciCore specific)
LMOD_DISABLE_SAME_NAME_AUTOSWAP="no"

#Load DIAMOND module (required)
module load DIAMOND/2.0.15-GCC-10.3.0

#Navigate to project directory
cd /scicore/home/ebertd/dexter0000/LSP

#Make working directories if they don't exist
mkdir -p genespace/genomeRepo/*
mkdir -p genespace/working_dir

#Make sure working directories are empty (old runs can cause conflicts)
rm -r genespace/genomeRepo/*
rm -r genespace/working_dir/*

#Organize genome files (this structure is mandatory and requires an index file of genome names)
while read GENOME; do
	mkdir -p genespace/genomeRepo/"$GENOME"
	cp annotations/"$GENOME".gff3 genespace/genomeRepo/"$GENOME"/"$GENOME".gff3
	cp annotations/"$GENOME".protein.faa genespace/genomeRepo/"$GENOME"/"$GENOME".faa
done <scripts/genome_index.txt

#Replace all dashes (-) in genome names with an underscore (_). Required by genespace.
find genespace/genomeRepo/ -depth -name '*-*' -execdir bash -c 'mv -i "$1" "${1//-/_}"' bash {} \;

#Load required module
module load R

#Start R session
Rscript /scicore/home/ebertd/dexter0000/LSP/scripts/genespace.r
````

### Genespace.r

This R script is called by the previous bash script. The completed run files are all saved to disk and can be reloaded in an interactive session for further plotting.

````R
#Load genespace package
library(GENESPACE)

#Set directory paths
# -- change paths to those valid on your system
genomeRepo <- "/scicore/home/ebertd/dexter0000/LSP/genespace/genomeRepo"
genespaceWd <- "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir"
path2mcscanx <- "/scicore/home/ebertd/dexter0000/software/MCScanX-master"

#Create list of genomes to use
genomeList <- list.dirs(path = "/scicore/home/ebertd/dexter0000/LSP/genespace/genomeRepo", full.names = FALSE, recursive = FALSE)

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = genomeList,
  genomeIDs = genomeList,
  genespaceWd = genespaceWd,
  troubleShoot = TRUE,
  headerEntryIndex = 1,
  gffIdColumn = "ID")

##########################
# -- initalize the run and QC the inputs
#Disable dotplot generation because it is very very slow.
#Oneway blast also speeds things up

  gpar <- init_genespace(
  wd = genespaceWd, 
  path2mcscanx = path2mcscanx,
  blkSize = 2,
  blkRadius = 50,
  dotplots = "never",
  onewayBlast = TRUE
  )
  
#Run orthofinder
out <- run_genespace(gpar, overwrite=TRUE)

#Add a custom background color
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "black"))

#Specify a region of interest for plotting
roi <- data.frame( 
genome = c(rep("t2_17_3_4i_13",4)), 
chr = c( "19","35","16","3")
)
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/chr5.pdf",
    width = 5,
    height = 5)
    
ripDat <- plot_riparian(
  out, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c( "19","35","16","3"),
  highlightBed = roi,
  backgroundColor = NULL)
  
dev.off()

pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/all.pdf",
    width = 5,
    height = 5)
````

### Curated riparian plot

````
srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

conda activate orthofinder

module load R

R

library(GENESPACE)

#Load genespace results
load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
################################################################
#Manual curation of genomes
################################################################
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(12,37,5))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(8,39,31,4,42,25))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,114,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(75,33,157,20,17))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(26,76,29,3,43))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(77,85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,108,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,448,64,81,43,31,338,104,111))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(68,76,35,151))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(46))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,208,82,662,35,56,98))

CH_2016_H_34<-data.frame(
  genome = c("CH_2016_H_34"), 
  chr = c(55,164,21,74,105))

US_D_3<-data.frame(
  genome = c("US_D_3"), 
  chr = c(126,7,98,27,38,106,77,9))

NCBI_scaffold<-data.frame(
  genome = c("NCBI_scaffold"), 
  chr = c(361461))

ET_C_1<-data.frame(
  genome = c("ET_C_1"), 
  chr = c(28,64,2))

DZ_JV_2<-data.frame(
  genome = c("DZ_JV_2"), 
  chr = c(26,189,73,116,75,103,19))

IL_TY_10<-data.frame(
  genome = c("IL_TY_10"), 
  chr = c(25,56,15))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, t3_12_3_1i_12, CH_t4_12_3_3, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49,NCBI_scaffold,CH_2016_H_34,US_D_3,ET_C_1,DZ_JV_2,IL_TY_10)

genomeIds <- c("D_similis_IL_SIM_A20",             
               "US_SP_221_1",
               "RU_RM1_2",
               "CN_W1_1",               
               "CH_H_2299",               
               "CH_434_inb3_a_1",
               "FI_SK_58_2_18_4",               
               "CH_H_2015_59",              
               "t3_12_3_1i_12",
               "CH_H_2015_49",
               "t2_17_3_4i_13"
              )

##############################################################################
#plots
##############################################################################
roi <- data.frame(
  genome = c("t2_17_3_4i_13"), 
  chr = c("191","161","31","351"))

# export plot
png(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/riparian.png",
    width = 8,
    height = 9, units="in", res = 900)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  syntenyWeight = 1,
  xlabel = sprintf("D. magna chromosome 5"),
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000,
  forceRecalcBlocks = FALSE,
  scalePlotWidth = 1,
  gapProp = 0.002,
  braidAlpha = 0.8,
  customRefChrOrder=c("191","161","31","351")
)
  
  dev.off()

````

### Minimal riparian plot

````
genomeIds <- c("t1_10_3_2",             
               "CH_434_inb3_a_1",
               "t2_17_3_4i_13"
              )

##############################################################################
#plots
##############################################################################
roi <- data.frame(
  genome = c("t2_17_3_4i_13"), 
  chr = c("191","161","31","351"))

# export plot
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/riparian_swisspond.pdf",
    width = 7,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  syntenyWeight = 1,
  xlabel = sprintf("D. magna chromosome 5"),
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000,
  forceRecalcBlocks = FALSE,
  scalePlotWidth = 1,
  gapProp = 0.002,
  braidAlpha = 0.8,
  customRefChrOrder=c("191","161","31","351")
)
  
  dev.off()
````



