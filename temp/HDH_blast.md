# HDH Haplotype BLAST

This script creates a BLAST database from a collection of D. magna genome assemblies and then finds the position of the chromosome 5 HRH region as defined by two 1Kb markers positioned where sequence homology resumes between assemblies.



### Request compute resources

````
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash
````



### Build the BLAST database


````
# First make sure that a combined fasta does not exist
rm blast/combined_assemblies.fasta

# The genome name has to be added to each contig and then all of the individual fasta concatenated together
for file in genomes_original/*.gz
do
    # Extract the filename without the .gz extension
    name=$(basename "$file" .gz)
    
    # Decompress, modify FASTA headers, and append to combined file
    zcat "$file" | awk -v name="$name" '/^>/{print ">" name "_" $0; next} {print}' >> blast/combined_assemblies.fasta
done

# Make a custom database from the concatenated fasta
module load BLAST
makeblastdb -in blast/combined_assemblies.fasta -dbtype nucl -out blast/daphnia_db
````



### Run local BLAST

````
blastn -query blast/LSPleft.fasta -db blast/daphnia_db -out blast/LSPleft_results.txt
blastn -query blast/LSPright.fasta -db blast/daphnia_db -out blast/LSPright_results.txt
````

