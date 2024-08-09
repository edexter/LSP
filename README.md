This repository contains all of the code relevant to my investigation into the large structural polymorphism found on the right arm of chromosome 5 of the Daphnia magna genome. There are several modules:



## Module 1: Chromosome 5 synteny analysis

This module performs a genome synteny analysis across a set of reference Daphnia genome assemblies using the GENESPACE software, and then plots the results for chromosome 5 using custom R scripts. 

**Component script files:**

* genespace.md (contains multiple sub-scripts)
* repeat_modeling.sh
* protein_extract.sh

**Required input files:**

* A set of reference genome assemblies
* An index file ("genome_index.txt") of the genome assembly file names

**Optional input files:**

* Pre-trained gene models or RNA-seq data for Augustus gene prediction. The model used in this analysis is named "Dmagna_iso"



## Module 2: Chromosome 5 dotplot alignment

This module extracts all candidate contigs from individual genome assembly which appear to belong to chromosome 5, and concatenates them together into individual chromosome 5 FASTA files. All pairwise combinations are then reciprocally aligned together using an interactive DGENIES session (based on MiniMap). The alignment files are then visualized as a matrix of pairwise alignments using custom R code.

**Component script files:**

* chromosome_5_scaffolding.md (contains multiple sub-scripts)
* dotplot_matrix_main.R
* dotplot_matrix_utils.R

**Required input files:**

* A set of reference genome assemblies from which to extract the chromosome 5 contigs

* DGenies output folders for each pairwise chromosome 5 alignment



## Module 3: Linkage Disequilibrium (LD) calculations

This module calculates LD across chromosome 5 using a previously created VCF file from 258 D. magna short-reads mapped against a single reference genome.

**Component script files:**

* linkage_disequilibrium.md (contains multiple sub-scripts)

**Required input files:**

* An input VCF file "merged_daphnia_filtered.vcf"



## Module 4: GC bias calculations

This module calculates GC content across the right arm of chromosome 5 for multiple genome assemblies using a custom perl script and then plots the results using a custom R script.

 **Component script files:**

* GC_calculations.md (contains multiple sub-scripts)

**Required input files:**

* An input VCF file "merged_daphnia_filtered.vcf"

## Module 5: HDH orthology analysis for Lake Aegelsee

This module performs a gene orthology analysis on the three HDH haplotypes which were observed in Lake Aegelsee based on the genome annotation files which were previously created by AUGUSTUS in the GENESCAPE module.

**Component script files:**

* orthology_analysis.md (contains multiple sub-scripts)

**Required input files:**

* Genome annotation files from Augustus (one for each genome)
  * "t2_17_3_4i_13.gff3"
  * "CH_434-inb3-a-1.gff3"
  * "t1_10_3_2.gff3"
* Gene model predictions from Augustus (one for each genome)
  * "t2_17_3_4i_13_functional.tsv"
  * "CH_434-inb3-a-1_functional.tsv"
  * "t1_10_3_2_functional.tsv"
* Orthogroup master list from Augustus
  * "Orthogroups.tsv"



## Module 6: HDH coordinate BLAST