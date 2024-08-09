# Convert GFF to BED

This script converts the GFF output files from the previous step into BED files, which are required for GENESPACE. This script only takes a moment to run so it can just be performed interactively.

````
srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash
module load BEDOPS

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/t2_17_3_4i_13.gff3 > t2_17_3_4i_13.bed

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/CH_434-inb3-a-1.gff3 > CH_434-inb3-a-1.bed

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/t1_10_3_2.gff3 > t1_10_3_2.bed

gff2bed < /scicore/home/ebertd/dexter0000/LSP/annotations/t3_12_3_1i_12.gff3 > t3_12_3_1i_12.bed
````

### 