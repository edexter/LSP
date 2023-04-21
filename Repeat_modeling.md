# Repeat modeling

Before repeats can be masked, they must be modeled. Modeling is the very slow part requiring several days with 32 threads. I'm using the CH14 genome to train the model because it's the most complete assembly.

### repeat_model.sh

````bash
#!/bin/bash

#SBATCH --job-name=repeat_model			#Job name
#SBATCH --cpus-per-task=32	        	#Number of cores reserved
#SBATCH --mem-per-cpu=1G            	#Memory reserved per core.
										#Total memory reserved: 32GB
#SBATCH --time=168:00:00	        	#Maximum time the job will run
#SBATCH --qos=1week           			#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/LSP/logs/repeat_model.out

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/LSP/logs/repeat_model.err

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
