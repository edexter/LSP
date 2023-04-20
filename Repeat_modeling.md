# Repeat modeling

Before repeats can be masked, they must be modeled. Modeling is the very slow part requiring  up to several days run with 32 threads. Since I'm not sure if it matters whether repeats are modeled from the whole genome, or just a subsample is OK I'm running scripts that produce models from just the 3 LSP haplotypes as well as the full CH14 genome assembly.



## Repeat modeling using the whole CH14 genome

This runs on a separate queue than the next script so that I can use more nodes.

````bash
#!/bin/bash

#SBATCH --job-name=repeat_model_CH14	#Job name
#SBATCH --cpus-per-task=32	        	#Number of cores reserved
#SBATCH --mem-per-cpu=1G            	#Memory reserved per core.
										#Total memory reserved: 32GB
#SBATCH --time=168:00:00	        	#Maximum time the job will run
#SBATCH --qos=1week           			#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/LSP/logs/repeat_model_Ch14.out

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/LSP/logs/repeat_model_Ch14.err

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

########################################################################
#Define variables
GENOME=/scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta

################################################################################
#Navigate to project folder
cd /scicore/home/ebertd/dexter0000/LSP

#Load required module
module load RepeatModeler

#Build database. [-pa] is the number of cores. This is very slow, so more is better.
BuildDatabase -name repeat_modeler/CH14 -engine ncbi "$GENOME"
RepeatModeler -database repeat_modeler/CH14 -engine ncbi -pa 32 -LTRStruct
````



## Repeat modeling for each LSP haplotype separately

This uses a job array to process each LSP haplotype from the same script, and they can run in one queue since the job should be be much faster than the whole genome.

````bash
#!/bin/bash

#SBATCH --job-name=repeat_model_LSP		#Job name
#SBATCH --cpus-per-task=32	        	#Number of cores reserved
#SBATCH --mem-per-cpu=1G            	#Memory reserved per core.
										#Total memory reserved: 32GB
#SBATCH --time=24:00:00	        		#Maximum time the job will run
#SBATCH --qos=1day           			#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/LSP/logs/repeat_model_%A_%a.out

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/LSP/logs/repeat_model_%A_%a.err

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
module load RepeatModeler

#Build database. [-pa] is the number of cores. This is very slow, so more is better.
BuildDatabase -name repeat_modeler/"$GENOME"_LSP -engine ncbi FASTA/"$GENOME"_LSP.fa
RepeatModeler -database repeat_modeler/"$GENOME"_LSP -engine ncbi -pa 16 -LTRStruct
````