# 3 geneome genespace run

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
Rscript /scicore/home/ebertd/dexter0000/LSP/scripts/genespace_small.r
````

### Genespace_small.r

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
genomeList <- c("t2_17_3_4i_13","CH_434_inb3_a_1","t1_10_3_2")

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



### #### custom plots

````
load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
###WORKS
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/test.pdf",
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
  #customRefChrOrder = c("19","35","16","3"),
  #highlightBed = roi,
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 5"),
  inversionColor = "green",
  genomeIDs = gsParam$genomeIDs[1:3] 
)
  
  dev.off()
  
###TEST2

roi <- data.frame(
  genome = c("t2_17_3_4i_13", "t2_17_3_4i_13"), 
  chr = c("191","351","161","31"))

invchr <- data.frame(
  genome = c("t1_10_3_2", "t1_10_3_2","t1_10_3_2","t1_10_3_2", "t1_10_3_2", "CH_434_inb3_a_1","CH_434_inb3_a_1","t1_10_3_2"), 
  chr = c(46, 23,641, 107,180, 5, 12,5))
  
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/test2.pdf",
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
  genomeIDs = c("t1_10_3_2","CH_434_inb3_a_1","t2_17_3_4i_13"),
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
````