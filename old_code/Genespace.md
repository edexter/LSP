Genespace

```r
install.packages(/scicore/home/ebertd/dexter0000/software/devtools_2.4.5.tar.gz, repos = NULL, type="source")
```

````
srun --nodes=1 --cpus-per-task=16 --mem=2G --pty bash

conda config --set ssl_verify no
conda create -n orthofinder
conda install -c bioconda orthofinder 

conda activate orthofinder
R #Need to open fro inside the conda environment

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("remotes")
devtools::install_github("jtlovell/GENESPACE")
````

````
bash install Anaconda
remotes::install_github("jtlovell/GENESPACE")

remotes::install_local(""/home/eric/software/GENESPACE-master.zip")
````

````
conda env create -n GENESPACE --file GENESPACE.yaml
````



### Run GENESPACE

````
conda activate GENESPACE
R
rm(list = ls())

library(GENESPACE)
###############################################
# -- change paths to those valid on your system
genomeRepo <- "/home/eric/scratch/genomeRepo"
genespaceWd <- "/home/eric/scratch/genespaceLSP"
path2mcscanx <- "/home/eric/software/MCScanX-master/"
###############################################

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("CH14","CH434"),
  genomeIDs = c("CH14","CH434"),
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
  blkRadius = 25)
  
# -- accomplish the run
out <- run_genespace(gpar, overwrite=TRUE)

#Improve the plotting
invchr <- data.frame(
  genome = c("CH434"), 
  chr = c("000012F/rc:1509248-4049114"))
  
ripDat <- plot_riparian(
  out, 
  refGenome = "CH14",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  invertTheseChrs = invchr)

````



#New genespace

````
conda activate GENESPACE
R
rm(list = ls())

library(GENESPACE)
###############################################
# -- change paths to those valid on your system
genomeRepo <- "/home/eric/scratch/genomeRepo"
genespaceWd <- "/home/eric/scratch/genespaceLSP"
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

#"T1","CH434","CH14" works
##########################
# -- initalize the run and QC the inputs

  gpar <- init_genespace(
  wd = genespaceWd, 
  path2mcscanx = path2mcscanx,
  blkSize = 2,
  blkRadius = 50
  )
  
# -- accomplish the run
out <- run_genespace(gpar, overwrite=TRUE)

#Improve the plotting
invchr <- data.frame(
  genome = c("CH434"), 
  chr = c("000012F/rc:1509248-4049114"))

customPal <- colorRampPalette(
  c("skyblue"))
  
pdf(file = "/home/eric/scratch/genespaceLSP/My Plot7.pdf",
    width = 4,
    height = 4)
ripDat <- plot_riparian(
  out, 
  refGenome = "CH434",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  palette = customPal)
  
  dev.off()

````

2 25 is pairwise

````
ripDat <- plot_riparian(
  out, 
  refGenome = "CH14",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE)
````



### Convert GFF to BED and perform additional formatting

````
srun --nodes=1 --cpus-per-task=16 --mem=2G --pty bash

#Load required module
module load BEDOPS

# Input file path
GENOME=CH434

#Convert GFF to BED
gff2bed < "$GENOME"_LSP.gff > temp_"$GENOME".bed

input_file=temp_"$GENOME".bed

# Set input and output file names
input_file=temp_"$GENOME".bed
output_file="$GENOME".bed

# Remove lines that don't contain "gene" in field 8 and replace field 1 with "Chr5"
awk '$8 ~ /gene/ {gsub($1, "Chr5", $1); print}' $input_file |

# Replace field 4 with last value in field 10 and print only fields 1-4
awk '{split($10, arr, "."); $4 = arr[length(arr)]; print $1"\t"$2"\t"$3"\t"$4}' > $output_file
````



````
rm(list = ls())

setwd("/home/eric/scratch/gggenomes")

library(gggenomes)

gggenomes(genes=emale_genes, seqs=emale_seqs) 

little_genes<-read_feats(list.files(("./"), "*.gff3", full.names=TRUE),fix_augustus_cds=TRUE)
little_genes <- little_genes %>%
  mutate(seq_id = feat_id)
  
  little_genes <- little_genes %>%
  filter_all(any_vars(str_detect(., "gene")))


# Read in fasta files & set order in the plot from top to bottom
little_seqs<-read_seqs(list.files(("./"), "*.fa", full.names=TRUE), parse_desc=FALSE)



gggenomes(genes=little_genes) 

#Works!
s0 <- read_seqs(ex("emales/emales.fna"))
g0 <- read_feats(ex("emales/emales.gff"))
gggenomes(g0, s0) +
  geom_seq() + geom_gene()
  
#My data (works!)
s0 <- read_seqs(list.files(("./"), "*.fna", full.names=TRUE), parse_desc=FALSE)
g0 <- little_genes<-read_feats(list.files(("./"), "*.gff3", full.names=TRUE),fix_augustus_cds=TRUE)

gggenomes(g0[1:100,], s0, wrap=2e5) +
  geom_seq() + geom_gene()
````

````
export AUGUSTUS_CONFIG_PATH=/scicore/home/ebertd/dexter0000/LSP/config/
module load AUGUSTUS

gtf2gff.pl <CH14_LSP.gff --out=CH14_LSP.gff3 --gff3
gtf2gff.pl <CH434_LSP.gff --out=CH434_LSP.gff3 --gff3
gtf2gff.pl <T1_LSP.gff --out=T1_LSP.gff3 --gff3
````

ripDat <- plot_riparian(
  out,
  refGenome = "CH434",
  useOrder = FALSE,
  useRegions = FALSE)

