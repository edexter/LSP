#Gene annotation whole chromosome 5

````bash
#!/bin/bash

#SBATCH --job-name=gene_predict			#Job name
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
#SBATCH --array=1-3%3

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

echo Repeat masking genome "$GENOME"

RepeatMasker FASTA/"$GENOME".fa -pa 16 -gff -xsmall -lib RM_120310.ThuApr201445402023/consensi.fa.classified  -dir FASTA 

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

Now copy files to EvoMartes server

````
#Request compute node for interactive session
srun --nodes=1 --cpus-per-task=16 --mem=4G --pty bash

conda activate GENESPACE
R
rm(list = ls())

library(GENESPACE)
###############################################
# -- change paths to those valid on your system
genomeRepo <- "/home/eric/scratch/genomeRepo"
genespaceWd <- "/home/eric/scratch/genespaceChr5"
path2mcscanx <- "/home/eric/software/MCScanX-master/"
###############################################

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("T1","CH14","CH434"),
  genomeIDs = c("T1","CH14","CH434"),
  genespaceWd = genespaceWd,
  troubleShoot = TRUE,
  headerEntryIndex = 1,
  gffIdColumn = "ID")

##########################
# -- initalize the run and QC the inputs

  gpar <- init_genespace(
  wd = genespaceWd, 
  path2mcscanx = path2mcscanx,
  blkSize = 2,
  blkRadius = 50
  )
 #blksize 5 also works fine 
# -- accomplish the run
out <- run_genespace(gpar, overwrite=TRUE)

#Improve the plotting
invchr <- data.frame(
  genome = c("CH434"), 
  chr = c("000012F/rc:1509248-4049114"))

customPal <- colorRampPalette(
  c("lightblue"))

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "black"))
  
roi <- data.frame(
genome = c("CH14", "CH14"), 
chr = c("ptg000019l_1/rc","ptg000035l_1"), 
color = c("#FAAA1D", "#17B5C5"))
  
pdf(file = "/home/eric/scratch/genespaceLSP/My Plot.pdf",
    width = 5,
    height = 6)
    
ripDat <- plot_riparian(
  out, 
  refGenome = "CH14",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  addThemes = ggthemes,
  syntenyWeight = 1,
  customRefChrOrder = c("ptg000019l_1/rc","ptg000035l_1","ptg000016l_1/rc","ptg000003l_1"),
  #highlightBed = roi,
  backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  inversionColor = "green")
  
  dev.off()
  

ripDat <- plot_riparian(
  gsParam = out, 
  highlightBed = roi, 
  backgroundColor = NULL, 
  genomeIDs = c("sandLizard", "chicken", "human", "mouse", "platypus"),
  refGenome = "human", 
  customRefChrOrder = c("X", 1:22))
````

/scicore/home/ebertd/dexter0000/software/MCScanX-master

```r
install.packages("devtools", repos='http://cran.us.r-project.org')
```