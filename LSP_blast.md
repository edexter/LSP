### Request compute resources

````
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash
````



### Unzip and concatenate the assemblies


````
# First make sure the file is empty
rm blast/combined_assemblies.fasta

#The filename has to be added to each contig
for file in genomes_original/*.gz
do
    # Extract the filename without the .gz extension
    name=$(basename "$file" .gz)
    
    # Decompress, modify FASTA headers, and append to combined file
    zcat "$file" | awk -v name="$name" '/^>/{print ">" name "_" $0; next} {print}' >> blast/combined_assemblies.fasta
done
````



### Build a local database from the assemblies

````
module load BLAST
makeblastdb -in blast/combined_assemblies.fasta -dbtype nucl -out blast/daphnia_db
````



### Run local blast

````
blastn -query blast/LSPleft.fasta -db blast/daphnia_db -out blast/LSPleft_results.txt
blastn -query blast/LSPright.fasta -db blast/daphnia_db -out blast/LSPright_results.txt
````

