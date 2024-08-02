# Gene annotation and synteny plots with whole D. magna genomes

### Format_genomes.sh

First upload all genome assemblies to [genomes_original]  directory and manually clean up the file names. This script then performs various formatting operations required for  the GENESPACE program. This script completes in a few minutes.

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

#Improve the plotting by inverting chromosomes as needed
#invchr <- data.frame(
#  genome = c("CH434"), 
#  chr = c("000012F/rc:1509248-4049114"))

#Add a custom color palette
#customPal <- colorRampPalette(
#  c("lightblue"))

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
    
ripDat <- plot_riparian(
  out, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  addThemes = ggthemes,
  syntenyWeight = 1,
  backgroundColor = NULL)
  
dev.off()
````



### Curated  chr5 riparian plot

````R
srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

conda activate orthofinder

module load R

R

library(GENESPACE)

load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
############################################################
#SMALLER SUBSET OF GOOD ONES
###########################################################
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(5,12,37))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,46,641,261,30,248,47,56,5,107,180,482))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49)

#T1 looks too bad "t1_10_3_2"
genomeIds <- c("D_similis_IL_SIM_A20","US_SP_221_1","Xinb3","NO_V_7","RU_RM1_2","CN_W1_1","t1_10_3_2","CH_H_2299","CH_434_inb3_a_1","FI_SK_58_2_18_4","CH_H_2015_49","t3_12_3_1i_12","CH_H_2015_59","t2_17_3_4i_13")

roi <- data.frame(
  genome = c("t2_17_3_4i_13", "t2_17_3_4i_13"), 
  chr = c("191","351","161","31"))
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/chr5_curated.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c("191","351","161","31"),
  #highlightBed = roi,
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
````

### LSP 1+2+3 riparian plot

````bash
R

library(GENESPACE)

load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
############################################################
#SMALLER SUBSET OF GOOD ONES
###########################################################
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(5,12,37))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,46,641,261,30,248,47,56,5,107,180,482))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49)

#T1 looks too bad "t1_10_3_2"
genomeIds <- c("t1_10_3_2","CH_434_inb3_a_1","t2_17_3_4i_13")

roi <- data.frame(
  genome = c("t2_17_3_4i_13", "t2_17_3_4i_13"), 
  chr = c("191","351","161","31"))
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/LSPSwisspond.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c("191","351","161","31"),
  #highlightBed = roi,
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
  
````



### All genomes

````
R

library(GENESPACE)

load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
############################################################
#SMALLER SUBSET OF GOOD ONES
###########################################################
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(5,12,37))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,46,641,261,30,248,47,56,5,107,180,482))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49)

#T1 looks too bad "t1_10_3_2"
genomeIds <- c("t1_10_3_2","CH_434_inb3_a_1","t2_17_3_4i_13")

roi <- data.frame(
  genome = c("t2_17_3_4i_13", "t2_17_3_4i_13"), 
  chr = c("191","351","161","31"))
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/chr5ALL.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c("191","351","161","31"),
  #highlightBed = roi,
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  #inversionColor = "white",
  #genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
````



### ABC locus

````bash
srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

conda activate orthofinder

module load R

R

library(GENESPACE)

load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
############################################################
#SMALLER SUBSET OF GOOD ONES
###########################################################
t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(46,161,51,241,11))
    
CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31,29,68,33,186,165,98,89))

t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157,89,186,33,164,37,126,38))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35,33,43,218,120,44,24,52))

FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42,14,2))
  
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(5,12,37,6,27))

CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97,2,3,82))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,46,641,261,30,248,47,56,5,107,180,482,75,140,103))

CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71,4))
  
RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,73,25))
  
NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11,4,86))
  
US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37,8))
  
CH_2016_H_34<-data.frame(
  genome = c("CH_2016_H_34"), 
  chr = c(175,331,69,242,75,61,269,17,132,251,83,87))
  
D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93,165,29,42,440,78,38,6,206,78))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129,49,83,32))
  
invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49, CH_t4_12_3_3,CH_2016_H_34)

#T1 looks too bad "t1_10_3_2"
genomeIds <- c("D_similis_IL_SIM_A20","US_SP_221_1","CH_2016_H_34","NO_V_7","CH_t4_12_3_3","CN_W1_1","CH_H_2299","CH_434_inb3_a_1","FI_SK_58_2_18_4","t3_12_3_1i_12","CH_H_2015_59","t2_17_3_4i_13","Xinb3")

roi <- data.frame(
  genome = c("Xinb3"), 
  chr = c("28","225","26","64","51","87","69","92","65","60","23","11"))
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/chr4_curated_H.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "Xinb3",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c("28","225","26","64","51","87","69","92","65","60","23","11"),
  #highlightBed = roi,
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
  
````

# Scaffolded LSP 1+2+3+refs

````
srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

conda activate orthofinder

module load R

R

library(GENESPACE)

load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
############################################################
#SMALLER SUBSET OF GOOD ONES
###########################################################
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(12,37,11,17))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,43,64,93,69,17,158))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(191,46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,35,261,30,248,47,56,107,180,482,5,218,64,159,156,155,82))

NCBI_scaffold<-data.frame(
  genome = c("NCBI_scaffold"), 
  chr = c(361461,100000111))

CH14_scaffold<-data.frame(
  genome = c("CH14_scaffold"), 
  chr = c(9))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49,NCBI_scaffold,CH14_scaffold)

#T1 looks too bad "t1_10_3_2"
genomeIds <- c("D_similis_IL_SIM_A20","t1_10_3_2","CH_434_inb3_a_1","t2_17_3_4i_13","CH14_scaffold","NCBI_scaffold","XINB_scaffold")

roi <- data.frame(
  genome = c("XINB_scaffold"), 
  chr = c("5","5_1"))
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/CHR5_1MB.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "XINB_scaffold",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c("5","5_1"),
  #highlightBed = roi,
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 500000
)
  
  dev.off()
  
````





#CHR6

````
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(12,37,11,17))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93,69,17))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(191,46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,35,261,30,248,47,56,107,180,482,5,218,64,159,156,155))

NCBI_scaffold<-data.frame(
  genome = c("NCBI_scaffold"), 
  chr = c(361461,100000111))

CH14_scaffold<-data.frame(
  genome = c("CH14_scaffold"), 
  chr = c(9))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49,NCBI_scaffold,CH14_scaffold)

genomeIds <- c("D_similis_IL_SIM_A20","t1_10_3_2","CH_434_inb3_a_1","t2_17_3_4i_13","CH14_scaffold","NCBI_scaffold","XINB_scaffold")

roi <- data.frame(
  genome = c("XINB_scaffold"), 
  chr = c("6"))

pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/CHR6_1MB.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "XINB_scaffold",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  #customRefChrOrder = c("5","5_1"),
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 500000
)
  
  dev.off()
````



# CH5 10 best genomes

````

srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

conda activate orthofinder

module load R

R

library(GENESPACE)

load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
############################################################
#SMALLER SUBSET OF GOOD ONES
###########################################################
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(5,12,37))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,46,641,261,30,248,47,56,5,107,180,482))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49)

#T1 looks too bad "t1_10_3_2"
genomeIds <- c("D_similis_IL_SIM_A20","US_SP_221_1","RU_RM1_2","CN_W1_1","CH_H_2299","CH_434_inb3_a_1","FI_SK_58_2_18_4","CH_H_2015_49","t3_12_3_1i_12","CH_H_2015_59","t2_17_3_4i_13")

roi <- data.frame(
  genome = c("t2_17_3_4i_13", "t2_17_3_4i_13"), 
  chr = c("191","351","161","31"))
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/chr5_curated_10.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c("191","351","161","31"),
  #highlightBed = roi,
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
````

# Chr 6 top 10

````
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(12,37,11,17))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42,13))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71,13,9))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157,10,137,76))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97,29,24,38,25))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93,69,17))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(191,46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,35,261,30,248,47,56,107,180,482,5,218,64,159,156,155))

NCBI_scaffold<-data.frame(
  genome = c("NCBI_scaffold"), 
  chr = c(361461,100000111))

CH14_scaffold<-data.frame(
  genome = c("CH14_scaffold"), 
  chr = c(9))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49,NCBI_scaffold,CH14_scaffold)

genomeIds <- c("D_similis_IL_SIM_A20","US_SP_221_1","RU_RM1_2","CN_W1_1","CH_H_2299","CH_434_inb3_a_1","FI_SK_58_2_18_4","CH_H_2015_49","t3_12_3_1i_12","CH_H_2015_59","t2_17_3_4i_13")

roi <- data.frame(
  genome = c("t2_17_3_4i_13"), 
  chr = c("151","131"))

pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/CHR6b.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  #customRefChrOrder = c("15_1","13_1"),
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 6"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
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

# Extract LSP flanking markers

````
samtools faidx /scicore/home/ebertd/dexter0000/LSP/genomes/t2_17_3_4i_13.fa "31:1534355-1535355" > query1.fa
samtools faidx /scicore/home/ebertd/dexter0000/LSP/genomes/t2_17_3_4i_13.fa "31:7110581-7111581" > query2.fa
````

````
makeblastdb -in /scicore/home/ebertd/dexter0000/LSP/genomes/t2_17_3_4i_13.fa -out daphnia_blast_db -parse_seqids -dbtype nucl

blastn -query query.fa -db daphnia_blast_db -out blastres.txt
blastn -query query2.fa -db daphnia_blast_db -out blastres2.txt

makeblastdb -in /scicore/home/ebertd/dexter0000/LSP/genomes/t3_12_3_1i_12.fa -out daphnia_blast_db_t3 -parse_seqids -dbtype nucl

blastn -query query.fa -db daphnia_blast_db_t3 -out blastres.txt
blastn -query query2.fa -db daphnia_blast_db_t3 -out blastres2.txt

makeblastdb -in /scicore/home/ebertd/dexter0000/LSP/genomes/CH_434-inb3-a-1.fa -out daphnia_blast_db_CH434 -parse_seqids -dbtype nucl

blastn -query query.fa -db daphnia_blast_db_CH434 -out blastres.txt
blastn -query query2.fa -db daphnia_blast_db_CH434 -out blastres2.txt

makeblastdb -in /scicore/home/ebertd/dexter0000/LSP/genomes/t1_10_3_2.fa -out daphnia_blast_db_T1 -parse_seqids -dbtype nucl

blastn -query query.fa -db daphnia_blast_db_T1
blastn -query query2.fa -db daphnia_blast_db_T1
````

### Control range

````
samtools faidx /scicore/home/ebertd/dexter0000/LSP/genomes/t2_17_3_4i_13.fa "191:3000000-3001000" > query3.fa
samtools faidx /scicore/home/ebertd/dexter0000/LSP/genomes/t2_17_3_4i_13.fa "191:7000000-7001000" > query4.fa

blastn -query query3.fa -db daphnia_blast_db
blastn -query query4.fa -db daphnia_blast_db

blastn -query query3.fa -db daphnia_blast_db_t3
blastn -query query4.fa -db daphnia_blast_db_t3

blastn -query query3.fa -db daphnia_blast_db_CH434
blastn -query query4.fa -db daphnia_blast_db_CH434

blastn -query query3.fa -db daphnia_blast_db_T1
blastn -query query4.fa -db daphnia_blast_db_T1
````

