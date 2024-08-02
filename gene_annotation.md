# Gene annotation for the LSP

This script will use an Augustus gene model already trained on the D. magna genome to predict genes for the different LSP haplotypes, and will then functionally annotate the predicted genes. The repeats have already been identified using repeat modeler (only needs to be performed once), and the LSP reference haplotypes will be repeat masked before running Augustus. InterProScan then operates on the amino acids sequences predicted by augustus to make a functional annotation.

````bash
#Request interactive node
srun --nodes=1 --cpus-per-task=16 --mem=1G --pty bash

#Load required module
module load RepeatMasker

#Need to provide the output from repeat modeler to repeat masker
ln -s RM_120310.ThuApr201445402023/consensi.fa
RepeatMasker -pa 16 -gff -nolow -lib consensi.fa scratch/CH14_chr5_R.fa

#Use this once the run has finished. 
#[-dir] specifies the output directory
#[-nolow] is specific to gene prediction downstream usage
#[-small] soft-masks repetitive regions instead of replacement with N's
ln -s RM_120310.ThuApr201445402023/consensi.fa.classified
RepeatMasker -pa 16 -gff -xsmall -qq -lib consensi.fa.classified  -dir FASTA FASTA/CH14_LSP.fa

#Load required module
module purge
module load AUGUSTUS

#Required setup steps before initial run
#cp -r /scicore/soft/apps/AUGUSTUS/3.4.0-foss-2021a/config .

#Needs to be added at every session if not permanently added to bash.rc
export AUGUSTUS_CONFIG_PATH=/scicore/home/ebertd/dexter0000/LSP/config/

#Run Augustus
augustus --species=Dmagna_iso FASTA/CH14_LSP.fa.masked gff3=on > annotations/CH14_LSP.gff3
````

### Gene_predict.sh

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

RepeatMasker FASTA/"$GENOME"_LSP.fa -pa 16 -gff -xsmall -lib RM_120310.ThuApr201445402023/consensi.fa.classified  -dir FASTA 

if [ -s FASTA/"$GENOME"_LSP.fa.masked  ]
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

augustus --species=Dmagna_iso FASTA/"$GENOME"_LSP.fa.masked gff3=on > annotations/"$GENOME"_LSP.gff3

if [ -s annotations/"$GENOME"_LSP.gff ]
then 
        echo "Augustus gene prediction for "$GENOME" completed"
else 
        printf "Augustus gene prediction for "$GENOME" did not complete successfully"
fi

#Extract protein sequences from annotation (note that the order of the flags matters!)

echo "Extracting protein sequences from annotation"

scripts/protein_extract.sh -p annotations/"$GENOME".protein.faa -c annotations/"$GENOME".CDS.faa annotations/"$GENOME"_LSP.gff

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







## Predict GO with interproscan

````
#Request interactive node
srun --nodes=1 --cpus-per-task=8 --mem=8G --pty bash

#Load required module
module load InterProScan

#Navigate to project folder
cd LSP

#Run interproscan. The [-dp] flag disables pre-calculated match lookup service.
interproscan.sh -i test.faa -cpu 8 -dp -goterms
````

I need:

-The number of genes in each haplotype

-The gene density in each haplotype

-Some GO statistics for each haplotype

-A synteny plot for all 3 haplotypes

# A little trick to add a prefix to fasta headers

````
cat reference-genome-families.fa | seqkit fx2tab | awk '{ print "abcDef1_"$0 }' | seqkit tab2fx > reference-genome-families.prefix.fa
````

