# Chromosome 5 scaffolding per genome

The following short code snippets allow for manual scaffolding of contigs into a single chromosome 5. The scaffolding is based on reciprocal minimap alignment of all contigs against multiple chromosome 5 assemblies. 

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



### RU-RM1-2

````bash
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

module load SAMtools

cd LSP/scratch

# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/RURM/RU-RM1-2.falcon.polish.fasta "000000F|arrow" --reverse-complement > 0.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/RURM/RU-RM1-2.falcon.polish.fasta "000085F|arrow" > 85.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/RURM/RU-RM1-2.falcon.polish.fasta "000036F|arrow" > 36.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/RURM/RU-RM1-2.falcon.polish.fasta "000087F|arrow"  > 87.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/RURM/RU-RM1-2.falcon.polish.fasta "000086F|arrow"  > 86.fa

cat 0.fa 85.fa 36.fa 87.fa 86.fa > RURM_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/RURM/RU-RM1-2.falcon.polish.fasta "000001F|arrow" --reverse-complement > RURM_chr5_R.fa

# Assemble both arms
cat RURM_chr5_L.fa RURM_chr5_R.fa > RURM_chr5.fa
````



### US-SP-221-1

````bash
srun --nodes=1 --cpus-per-task=1 --mem=8G --pty bash

module load SAMtools

cd LSP/scratch

# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/USSP/US-SP-221-1.falcon.polish.fasta "000006F|arrow" --reverse-complement > 6.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/USSP/US-SP-221-1.falcon.polish.fasta "000037F|arrow" --reverse-complement > 37.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/USSP/US-SP-221-1.falcon.polish.fasta "000020F|arrow" > 20.fa

cat 6.fa 37.fa 20.fa > USSP_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/USSP/US-SP-221-1.falcon.polish.fasta "000003F|arrow" > USSP_chr5_R.fa

# Assemble both arms
cat USSP_chr5_L.fa USSP_chr5_R.fa > USSP_chr5.fa
````



### CH-2016-H-34

````bash
# Assemble left arm

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2016/CH-2016-H-34.falcon.polish.fasta "000028F|arrow" --reverse-complement > CH2016_contig28.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2016/CH-2016-H-34.falcon.polish.fasta "000000F|arrow" > CH2016_contig0.fa

cat CH2016_contig28.fa CH2016_contig0.fa > CH2016_chr5_R.fa

# Assemble both arms
CH2016_chr5_L.fa CH2016_chr5_R.fa > CH2016_chr5.fa
````



### CH_H_2015_49

````bash
# Assemble left arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2015/CH_H_2015_49.hifiasm.fasta "ptg000076l" > 76.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2015/CH_H_2015_49.hifiasm.fasta "ptg000033l" --reverse-complement > 33.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2015/CH_H_2015_49.hifiasm.fasta "ptg000049l" --reverse-complement > 49.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2015/CH_H_2015_49.hifiasm.fasta "ptg000068l" > 68.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2015/CH_H_2015_49.hifiasm.fasta "ptg000108l" > 108.fa

samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2015/CH_H_2015_49.hifiasm.fasta "ptg000086l" > 86.fa

cat 76.fa 33.fa 49.fa 68.fa 108.fa 86.fa > CH2015_chr5_L.fa

# Assemble right arm
samtools faidx /scicore/home/ebertd/dexter0000/daphnia_ref/CH2015/CH_H_2015_49.hifiasm.fasta "ptg000035l" --reverse-complement > CH2015_chr5_R.fa

# Assemble both arms
cat CH2015_chr5_L.fa CH2015_chr5_R.fa > CH2015_chr5.fa
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
