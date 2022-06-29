# aln2pheno

**A bioinformatics tool for rapid screening of genetic features (nt or aa changes) potentially linked to specific phenotypes**


INTRODUCTION
------------------

The **aln2pheno** module screens a amino acid/nucleotide alignment against a "genotype-phenotype" database and reports the repertoire of mutations of interest per sequences and their potential impact on phenotype.


INPUT FILES
----------------

Table database in .tsv or .xlsx format 

**AND**        

Alignment (nucleotide or amino acid)   


(Mutation numbering will refer to the user-defined reference sequence included in the alignment).


USAGE
----------------

```bash

usage: algn2pheno.py [-h] [--db DB] [--sheet SHEET] [--table TABLE]
                     [--gencol GENCOL] [--phencol PHENCOL] -g GENE --algn ALGN
                     -r REFERENCE [--nucl] [--odir ODIR] [--output OUTPUT]
                     [--log LOG] [-f]

parse arguments

optional arguments:
  -h, --help            show this help message and exit
  --db DB               phenotype db. if excel, provide sheet name and columns
                        numbers for genotype and phenotype data. if not excel,
                        provide path to db, ensure 3 columns 'Mutation',
                        'Phenotype category' & 'Flag'
  --sheet SHEET, -s SHEET
                        Give sheet name (gene name or 'Lineages'). excel
                        input.
  --table TABLE, -t TABLE
                        table name in db. default: pokay. sqlite input.
  --gencol GENCOL       nr of column with genotype data. excel input.
  --phencol PHENCOL     nr of column with phenotype data. excel input.
  -g GENE, --gene GENE  Set gene or protein prefix
  --algn ALGN           Input alignment file
  -r REFERENCE, --reference REFERENCE
                        Give full reference sequence as in alignment file (use
                        quotation marks if this contains a space)
  --nucl                provide if nucleotide alignment instead of amino acid.
  --odir ODIR, -o ODIR  output directory
  --output OUTPUT       Set output file prefix
  --log LOG             logfile
  -f, --force           overwrite existing files


```


How to run (examples)
----------------------

1. database (.tsv) + amino acid alignment (SARS-CoV-2 Spike)

```bash
python algn2pheno.py --db database.tsv --algn alignment_aa_Spike.fasta -g S -r reference_header --odir output_folder --output output_prefix
```

2. database (.xlsx) + amino acid alignment (SARS-CoV-2 Spike)

```bash
python algn2pheno.py --db database.xlsx --algn alignment_aa_Spike.fasta --sheet S --gencol ["Mutation" column number] --phencol ["Phenotype category" column number] -g S -r reference_header --odir output_folder --output output_prefix
```

REQUIREMENTS
-----------------

The input database must contain mutations and associated phenotypes, listed in columns named `Mutation` and `Phenotype category` respectively. An additional column `Flag` may be provided for the ".tsv" input format.

In the case of using an excel file as input, an intermediary .tsv database will be created. A "Flag" column will be generated automatically with Flags "D" and "P" to differentiate mutation(s) Directly associated with the phenotype ("D") and mutation(s) Partially associated with the phenotype (i.e., the mutation is Part of a constellation associated with that phenotype) ("P"). Constellation of mutations should be included in the same row separated by commas (e.g., H69del,H70del,L452R,N501Y). 

Undefined amino acids should be included as "X" in the protein alignment.
Undefined nucleotides should be included as "N" in the gene alignment.

**Modules:**

> Main alignment to phenotype.
  - python=3.7
  - biopython
  - numpy 
  - pandas

> Further modules required to scrape Pokay database.
  - openpyxl
  - requests
  - pygithub
  - sqlalchemy
  - beautifulsoup4
  - lxml


INSTALLATION
-----------------

- git clone this repository.
- install and activate environment with required modules (see _pyenv.yml_).

CONFIGURATION
-----------------

This module does not require configuration.   



MAIN OUTPUTS
------------------------

> **_full_mutation_report.tsv**: provides the repertoire of "Flagged mutations" (i.e., database mutations detected in the alignment), the list of "phenotypes" supported by those mutations of interest and the list of "All mutations" detected in alignment for each sequence.

> **_flagged_mutation_report.tsv**: "Flagged mutation" binary matrix for all sequences and the "associated" phenotypes.


Other intermediate outputs:

> _all_mutation_matrix:  matrix of all mutations found in the alignment [mutations coded as 1 (presence), 0 (absence)]. At this stage, undefined amino acids ("X") or undefined nucleotides ("n") are not yet handled as no coverage positions.

> _all_mutation_report: matrix of all mutations found in the alignment [mutations coded as 1 (presence), 0 (absence) or nc (no coverage) positions] and associated phenotypes for the mutations in the database.

> _all_db_mutations_report: matrix of mutations in the database [mutations coded as 1 (presence), 0 (absence) or nc (no coverage) positions] and associated phenotypes, regardless if they were found or not in the alignment

> database.CLEAN.tsv: conversion of the ".xlsx" database into a clean ".tsv" format adding the gene "prefix" (-g argument) to each mutation (similar to when providing a .tsv database with the three headers 'Mutation', 'Phenotype category' and 'Flag')

MAINTAINERS
----------------

Current maintainers:

- Joao Santos (santosjgnd) 

- Bioinformatics Unit of the Portuguese National Institute of Health Dr. Ricardo Jorge (INSA) (insapathogenomics)



CITATION
----------

If you run algn2pheno, please cite:

João D. Santos,  Carlijn Bogaardt, Joana Isidro, João P. Gomes, Daniel Horton, Vítor Borges (2022). Algn2pheno: A bioinformatics tool for rapid screening of genetic features (nt or aa changes) potentially linked to specific phenotypes.


FUNDING
----------------

This work was supported by funding from the European Union’s Horizon 2020 Research and Innovation programme under grant agreement No 773830: One Health European Joint Programme.

