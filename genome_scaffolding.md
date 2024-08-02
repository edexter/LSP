# Installation

The installation can be tricky because of the usual issues with dependencies.

````bash
conda install mamba python=3.9
mamba install -c bioconda -c conda-forge longstitch
export PATH="/home/eric/anaconda3:$PATH" #If not already in path
````



## Scaffolding CH2299 with CH2299 nanopore reads

````
#Request node to run in interactive mode
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

#Navigate to working directory
cd scratch

#Create a directory for run
mkdir CH2299
cd CH2299

#Activate conda environment
conda activate TEST

#Make symlinks for input files (required). Note that file input names are extremely strict including suffixes. The read file must be named end with .fq or .fq.gz and the reference must end with .fa or .fa.gz, but only the base name can be passed to longstitch

ln -s /scicore/home/ebertd/GROUP/oxford_nanopore_pascal/CH-H-2014-2299_Mitosporidium_daphniae/CH-H-2014-2299_Md-2.fastq.gz reads.fq.gz

ln -s /scicore/home/ebertd/dexter0000/LSP/FASTA/CH-H-2299.fa.masked reference.fa

longstitch run draft="reference" reads="reads" G=200000000

####Scaffolding only mode
#Create a directory for run
mkdir CH2299_B
cd CH2299_B

#Activate conda environment
conda activate TEST

#Make symlinks for input files (required). Note that file input names are extremely strict including suffixes. The read file must be named end with .fq or .fq.gz and the reference must end with .fa or .fa.gz, but only the base name can be passed to longstitch

ln -s /scicore/home/ebertd/GROUP/oxford_nanopore_pascal/CH-H-2014-2299_Mitosporidium_daphniae/CH-H-2014-2299_Md-2.fastq.gz reads.fq.gz

ln -s /scicore/home/ebertd/dexter0000/LSP/FASTA/CH-H-2299.fa.masked reference.fa

longstitch ntLink-arks draft="reference" reads="reads" G=200000000

````

## Scaffolding CH14 with CH2299 nanopore reads

````
#Request node to run in interactive mode
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

#Activate conda environment
conda activate TEST

#Navigate to working directory
cd scratch

#Create a directory for run
mkdir CH14
cd CH14

#Make symlinks for input files (required). Note that file input names are extremely strict including suffixes.

ln -s /scicore/home/ebertd/GROUP/oxford_nanopore_pascal/FI-OER-3-3_Hamiltosporidium_tvaerminnensis/oxford_nanopore.fastq.gz reads.fq.gz

ln -s /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta reference.fa

longstitch ntLink-arks draft="reference" reads="reads" G=200000000
````



```
pip install assemblyStatistics

assemblyStatistics reference.k32.w100.tigmint-ntLink.longstitch-scaffolds.fa
```

$Ragtag

````
#Request node to run in interactive mode
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

# install with conda
conda create --name RAGTAG
conda activate RAGTAG
conda install -c bioconda ragtag

QUERY=t2_17_3_4i_13.fasta
REF1=t3_12_3_1i_12.fasta
REF2=CH_H_2015_49.fasta

# scaffold with multiple references/maps
ragtag.py scaffold -t 8 -o out_1 $REF1 $QUERY
ragtag.py scaffold -t 8 -o out_2 $REF2 $QUERY
ragtag.py merge $QUERY out_*/*.agp


#Move XINB_new to folder
cp /scicore/home/ebertd/dexter0000/daphnia_ref/XINB_chromsome/D_magna_XINB3_chromosomal_scaffolds.fasta XINB3.fasta

ragtag.py scaffold -t 8 -o XINB XINB3.fasta $QUERY

cd XINB/

assemblyStatistics ragtag.scaffold.fasta


#NCBI
cp /scicore/home/ebertd/dexter0000/daphnia_ref/NCBI/NCBI_genome.fa NCBI.fasta

ragtag.py scaffold -t 8 -o NCBI NCBI.fasta $QUERY

cd NCBI/

assemblyStatistics ragtag.scaffold.fasta

ragtag.py merge $QUERY NCBI/ragtag.scaffold.agp XINB/ragtag.scaffold.agp out_*/*.agp

assemblyStatistics ragtag.merge.fasta
````



## Scaffolding with just chromosome-level genomes

CH14

````bash
#Request node to run in interactive mode
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

#Activate conda environment
conda activate RAGTAG

#Navigate to project folder
cd scratch/ragtag

#########################
#CH14
#########################
QUERY=t2_17_3_4i_13.fasta
REF1=XINB.fasta
REF2=NCBI.fasta

ragtag.py scaffold -t 8 -o CH14_XINB $REF1 $QUERY
ragtag.py scaffold -t 8 -o CH14_NCBI $REF2 $QUERY
ragtag.py merge $QUERY CH14_XINB/ragtag.scaffold.agp CH14_NCBI/ragtag.scaffold.agp -o CH14_BOTH

assemblyStatistics CH14_BOTH/ragtag.merge.fasta

#########################
#T1
#########################
mv /scicore/home/ebertd/dexter0000/LSP/genomes_original/t1_10_3_2.gz t1_10_3_2.fasta.gz
gunzip t1_10_3_2.fasta.gz

QUERY=t1_10_3_2.fasta
REF1=XINB.fasta
REF2=NCBI.fasta

ragtag.py scaffold -t 8 -o T1_XINB $REF1 $QUERY
ragtag.py scaffold -t 8 -o T1_NCBI $REF2 $QUERY
ragtag.py merge $QUERY T1_XINB/ragtag.scaffold.agp T1_NCBI/ragtag.scaffold.agp -o T1_BOTH

assemblyStatistics T1_BOTH/ragtag.merge.fasta

#########################
#CH434
#########################

mv /scicore/home/ebertd/dexter0000/LSP/genomes_original/CH_434-inb3-a-1.gz CH_434-inb3-a-1.fasta.gz
gunzip CH_434-inb3-a-1.fasta.gz

QUERY=CH_434-inb3-a-1.fasta
REF1=XINB.fasta
REF2=NCBI.fasta

ragtag.py scaffold -t 8 -o CH434_XINB $REF1 $QUERY
ragtag.py scaffold -t 8 -o CH434_NCBI $REF2 $QUERY
ragtag.py merge $QUERY CH434_XINB/ragtag.scaffold.agp CH434_NCBI/ragtag.scaffold.agp -o CH434_BOTH

assemblyStatistics CH434_BOTH/ragtag.merge.fasta

#####
#Move to genome folder and prepare for annotation pipeline
####
cd 
/scicore/home/ebertd/dexter0000/LSP/genomes_original

cp /scicore/home/ebertd/dexter0000/scratch/ragtag/CH14_BOTH/ragtag.merge.fasta ./CH14_scaffold
gzip CH14_scaffold

cp /scicore/home/ebertd/dexter0000/scratch/ragtag/T1_BOTH/ragtag.merge.fasta ./T1_scaffold
gzip T1_scaffold

cp /scicore/home/ebertd/dexter0000/scratch/ragtag/T1_BOTH/ragtag.merge.fasta ./CH434_scaffold
gzip CH434_scaffold

cp /scicore/home/ebertd/dexter0000/daphnia_ref/XINB_chromsome/D_magna_XINB3_chromosomal_scaffolds.fasta ./XINB_scaffold
gzip XINB_scaffold

cp /scicore/home/ebertd/dexter0000/daphnia_ref/NCBI/NCBI_genome.fa ./NCBI_scaffold
gzip NCBI_scaffold

#Make index for pipeline
#Format genomes
#Annotate genomes
#Run genespace

````

#Correction test

````bash
ragtag.py correct XINB.fasta t2_17_3_4i_13.fasta -o CH14_corr

QUERY=CH14_corr/ragtag.correct.fasta
REF1=XINB.fasta
REF2=NCBI.fasta

ragtag.py scaffold -t 8 -o CH14_corr_XINB $REF1 $QUERY
ragtag.py scaffold -t 8 -o CH14_corr_NCBI $REF2 $QUERY
ragtag.py merge $QUERY CH14_corr_XINB/ragtag.scaffold.agp CH14_corr_NCBI/ragtag.scaffold.agp -o CH14_corr_BOTH
````

