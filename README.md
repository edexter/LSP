This repository contains all of the code relevant to my investigation into the large structural polymorphism found on the right arm of chromosome 5 of the Daphnia magna genome. There are several facets of this:

### Genome annotation

The goal in this module is to identify and functionally annotate the genes present in the different LSP haplotypes. This requires several input files: a Reference genome assemblies containing different LSP haplotypes and Augustus gene prediction model parameters pre-trained for D. magna

* Step 1: Create fasta files containing just the LSP sequence from each of the reference assemblies
* Step 2: Create a catalog of repetitive genetic elements to be masked before gene prediction
* Step 3: Mask the repetitive elements and predict genes using a pre-trained Augustus model
* Step 4: Extract predicted protein sequences and functionally annotate using InterProScan 