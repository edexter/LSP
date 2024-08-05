# Chromosome 5 scaffolding per genome

The following short code snippets allow for manual scaffolding of contigs into a single chromosome 5 in order to visualize synteny via dotplot alignment. The scaffolding is based on reciprocal Minimap alignment of all contigs against multiple chromosome 5 assemblies. 

### t2_17.3_4i_13 (Haplotype 1 reference)

````bash
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

module load SAMtools

cd LSP/scratch

# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta "ptg000019l_1" --reverse-complement > 19.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta "ptg000035l_1"> 35.fa

cat 19.fa 35.fa > CH14_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta "ptg000003l_1" > CH14_contig3_temp.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH14/t2_17.3_4i_13.v0.1.fasta "ptg000016l_1" --reverse-complement > CH14_contig16_temp.fa

cat CH14_contig16_temp.fa CH14_contig3_temp.fa > CH14_chr5_R.fa

# Assemble both arms
cat CH14_chr5_L.fa CH14_chr5_R.fa > CH14_chr5.fa
````



### CH_434-inb3-a-1 (Haplotype 2 reference)

````bash
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

module load SAMtools

cd LSP/scratch

# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/alt2/CH_434-inb3-a-1.falcon.polish.mask.fasta "000005F" > 5.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/alt2/CH_434-inb3-a-1.falcon.polish.mask.fasta "000019F" > 19.fa

cat 5.fa 19.fa > CH434_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/alt2/CH_434-inb3-a-1.falcon.polish.mask.fasta "000012F" --reverse-complement > CH434_contig12_reverse.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/alt2/CH_434-inb3-a-1.falcon.polish.mask.fasta "000037F" --reverse-complement > CH434_contig37_reverse.fa

cat CH434_contig37_reverse.fa CH434_contig12_reverse.fa > CH434_chr5_R.fa

# Assemble both arms
cat CH434_chr5_L.fa CH434_chr5_R.fa > CH434_chr5.fa
````



### t1_10.3_2 (Haplotype 3 reference)

````bash
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

module load SAMtools

cd LSP/scratch

# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000047l" > 47.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000082l" --reverse-complement > 82.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000035l" --reverse-complement > 35.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000046l" > 46.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000005l" --reverse-complement > 5.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000107l" --reverse-complement > 107.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000180l" --reverse-complement > 180.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000148l" > 148.fa

cat 47.fa 82.fa 35.fa 46.fa 5.fa 107.fa 180.fa 148.fa > T1_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T1/t1_10.3_2.hifiasm.fasta "ptg000023l" --reverse-complement > T1_chr5_R.fa

# Assemble both arms
cat T1_chr5_L.fa T1_chr5_R.fa > T1_chr5.fa
````



### t3_12.3_1i_12 (Haplotype 1 alternate reference)

````bash
# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000020l" > 20.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000087l" --reverse-complement > 87.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000033l" > 33.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000017l" --reverse-complement > 17.fa

cat 20.fa 87.fa 33.fa 17.fa > T3_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000305l" > 305.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000144l" > 144.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000105l" > 105.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000075l" --reverse-complement > 75.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/T3/t3_12.3_1i_12.31102021.fasta "ptg000023l" > 23.fa

cat 305.fa 144.fa 105.fa 75.fa 23.fa > T3_chr5_R.fa

# Assemble both arms
cat T3_chr5_L.fa T3_chr5_R.fa > T3_chr5.fa
````



### FI-SK-58-2-18-4 (Haplotype 2 alternate reference)

````bash
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

module load SAMtools

cd LSP/scratch

# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/FISK/FI-SK-58-2-18-4.falcon.polish.fasta "000004F|arrow" > 4.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/FISK/FI-SK-58-2-18-4.falcon.polish.fasta "000042F|arrow" --reverse-complement > 42R.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/FISK/FI-SK-58-2-18-4.falcon.polish.fasta "000039F|arrow" > 39.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/FISK/FI-SK-58-2-18-4.falcon.polish.fasta "000025F|arrow" --reverse-complement > 25R.fa

cat 4.fa 42R.fa 39.fa 25R.fa > FISK_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/FISK/FI-SK-58-2-18-4.falcon.polish.fasta "000008F|arrow" --reverse-complement > FISK_chr5_R.fa

# Assemble both arms
cat FISK_chr5_L.fa FISK_chr5_R.fa > FISK_chr5.fa
````



### CH-H-2299 (Haplotype 3 alternate reference)

````bash
# Assemble both arms
samtools faidx CH-H-2299.fa "000001F|arrow" --reverse-complement > 1R.fa

samtools faidx CH-H-2299.fa "000099F|arrow" > 99.fa

samtools faidx CH-H-2299.fa "000132F|arrow" --reverse-complement> 132R.fa

samtools faidx CH-H-2299.fa "000077F|arrow" > 77.fa

samtools faidx CH-H-2299.fa "000076F|arrow" > 76.fa

samtools faidx CH-H-2299.fa "000085F|arrow" --reverse-complement > 85R.fa

samtools faidx CH-H-2299.fa "000070F|arrow" > 70.fa

samtools faidx CH-H-2299.fa "000066F|arrow" > 66.fa

samtools faidx CH-H-2299.fa "000118F|arrow" > 118.fa

samtools faidx CH-H-2299.fa "000030F|arrow" > 30.fa

samtools faidx CH-H-2299.fa "000075F|arrow" > 75.fa

samtools faidx CH-H-2299.fa "000067F|arrow" > 67.fa

samtools faidx CH-H-2299.fa "000015F|arrow" > 15.fa

cat 1R.fa 99.fa 132R.fa 77.fa 76.fa 85R.fa 70.fa 66.fa 118.fa 30.fa 75.fa 67.fa 15.fa > CH2299_chr5.fa
````

