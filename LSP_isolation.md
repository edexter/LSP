# Create LSP-only FASTA files

For each of the three LSP haplotypes we need to isolate the LSP region from the rest of the genome. Note that we have to flip the orientation in CH434 and T1 relative to CH14 due to the random orientation that the assemblers placed the containing contigs in. This is quick and can be run interactively.

````bash
#Request interactive node
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

#Load required module
module load SAMtools

#Navigate to project directory
cd LSP

#Extract the entire contig containing the LSP from the 3 genomes and flip if nessessary
#Also create a separate FASTA containing just the LSP

#CH14
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta "ptg000003l_1" > FASTA/CH14_contig3.fa

samtools faidx FASTA/CH14_contig3.fa "ptg000003l_1:1535355-7110581" > FASTA/CH14_LSP.fa

#CH434
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/alt2/CH_434-inb3-a-1.falcon.polish.mask.fasta "000012F" --reverse-complement > FASTA/CH434_contig12_rc.fa

samtools faidx FASTA/CH434_contig12_rc.fa "000012F/rc:1509248-4049114" > FASTA/CH434_LSP.fa

#T1
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000023l" --reverse-complement > FASTA/T1_contig23_rc.fa

samtools faidx FASTA/T1_contig23_rc.fa "ptg000023l/rc:2108160-5172826" > FASTA/T1_LSP.fa
````
