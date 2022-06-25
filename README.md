# Panta: a fast and progressive approach for bacterial pan-genome analysis
## Overview
Panta is a tool for bacterial pan-genome analysis. Our progressive algorithm enables fast-forward construction of an incrementally added-up collection without the need to re-run the whole pipeline again whenever new samples are available. Panta can handle large collections with thousands of samples, which requires less time and memory than other existing tools. Genome annotation step is also integrated into Panta, so we do not need to annotate every genome individually. Input data of Panta could be both genome annotation (GFF) and genome assembly (FASTA).

Panta is written in python, it includes the followings dependencies:
 * python (3.7)
 * BEDTools (2.30.0)
 * Prodigal (2.6.3)
 * CD-HIT (4.8.1)
 * blast (2.12.0)
 * DIAMOND (2.0.14)
 * HMMER (3.3.2)
 * MCL (14.137)
 * abPOA (1.4.0)
 * mafft (7.490)
 * parallel (20220222)


## Installation
The simplest method is installed via conda:

1. Download and install the appropriate conda, such as miniconda from [hear](https://docs.conda.io/en/latest/miniconda.html)
   
   
2. Create a conda environment with all the necessary dependencies: From the repository directory run

```bash

conda create -y -c conda-forge -c defaults --name panta python=3.7 mamba

source activate panta

mamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file requirements.txt

```

## Usage
```
usage: panta.py [initiate/add] [-h] [-g [GFF [GFF ...]]] [-b [FASTA [FASTA ...]]]
                                [-f TSV] -o OUTDIR [-d] [-i IDENTITY] [--LD LD]
                                [--AL AL] [--AS AS] [-e EVALUE] [-t THREADS]
                                [--table TABLE]

Initiate: run initial pan-genome analysis
Add: add samples into previous collection

optional arguments:
  -h, --help            show this help message and exit
  -g [GFF [GFF ...]], --gff [GFF [GFF ...]]
                        gff input files (default: None)
  -b [FASTA [FASTA ...]], --fasta [FASTA [FASTA ...]]
                        assembly input files (default: None)
  -f TSV, --tsv TSV     tsv input file (default: None)
  -o OUTDIR, --outdir OUTDIR
                        output directory (default: None)
  -d, --diamond         use Diamond for all-agaist-all alignment instead of
                        Blastp (default: False)
  -i IDENTITY, --identity IDENTITY
                        minimum percentage identity (default: 0.95)
  --LD LD               length difference cutoff between two sequences
                        (default: 0)
  --AL AL               alignment coverage for the longer sequence (default:
                        0)
  --AS AS               alignment coverage for the shorter sequence (default:
                        0)
  -e EVALUE, --evalue EVALUE
                        Blast evalue (default: 1e-06)
  -t THREADS, --threads THREADS
                        number of threads to use, 0 for all (default: 0)
  --table TABLE         codon table (default: 11)
```