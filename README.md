# Pan-genome pipeline
## Installation
Depencencies:
- bedtools
- CD-HIT
- BLAST
- DIAMOND
- MCL
- pandas
- SeqIO

The simplest method is installed via conda:

1. Download and install the appropriate conda, such as miniconda from [here](https://docs.conda.io/en/latest/miniconda.html)
   
   
2. Create a conda environment with all the necessary dependencies: From the repository directory run

```bash

conda create -y -c conda-forge -c defaults --name panta python=3.10 mamba

conda activate panta

mamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file requirements.txt

```

## Usage
Activate conda enviroment:
```
source activate panta
```
### Main pipeline: run pan-genome analysis for the first time
```
usage: pan-genome.py main [-h] [-g [GFF ...]] [-f TSV] -o OUTDIR [-s] [-b {diamond,blast}] [-i IDENTITY] [--LD LD] [--AL AL] [--AS AS] [-e EVALUE]
                          [-t THREADS] [--table TABLE] [-a [{nucleotide,protein} ...]]

Main pipeline: run pan-genome analysis for the first time

options:
  -h, --help            show this help message and exit
  -g [GFF ...], --gff [GFF ...]
                        gff input files (default: None)
  -f TSV, --tsv TSV     tsv input file (default: None)
  -o OUTDIR, --outdir OUTDIR
                        output directory (default: None)
  -s, --dont-split      dont split paralog clusters (default: False)
  -b {diamond,blast}, --blast {diamond,blast}
                        method for all-against-all alignment (default: diamond)
  -i IDENTITY, --identity IDENTITY
                        minimum percentage identity (default: 0.7)
  --LD LD               length difference cutoff between two sequences (default: 0)
  --AL AL               alignment coverage for the longer sequence (default: 0)
  --AS AS               alignment coverage for the shorter sequence (default: 0)
  -e EVALUE, --evalue EVALUE
                        Blast evalue (default: 1e-06)
  -t THREADS, --threads THREADS
                        number of threads to use, 0 for all (default: 0)
  --table TABLE         codon table (default: 11)
  -a [{nucleotide,protein} ...], --alignment [{nucleotide,protein} ...]
                        run alignment for each gene cluster (default: None)

```
### Add pipeline: add sample into previous collection
```
usage: pan-genome.py add [-h] [-g [GFF ...]] [-f TSV] -c COLLECTION_DIR [-s] [-b {diamond,blast}] [-i IDENTITY] [--LD LD] [--AL AL] [--AS AS]
                         [-e EVALUE] [-t THREADS] [--table TABLE] [-a [{nucleotide,protein} ...]]

Add pipeline: add sample into previous collection

options:
  -h, --help            show this help message and exit
  -g [GFF ...], --gff [GFF ...]
                        gff input files (default: None)
  -f TSV, --tsv TSV     tsv input file (default: None)
  -c COLLECTION_DIR, --collection-dir COLLECTION_DIR
                        previous collection directory (default: None)
  -s, --dont-split      dont split paralog clusters (default: False)
  -b {diamond,blast}, --blast {diamond,blast}
                        method for all-against-all alignment (default: diamond)
  -i IDENTITY, --identity IDENTITY
                        minimum percentage identity (default: 0.7)
  --LD LD               length difference cutoff between two sequences (default: 0)
  --AL AL               alignment coverage for the longer sequence (default: 0)
  --AS AS               alignment coverage for the shorter sequence (default: 0)
  -e EVALUE, --evalue EVALUE
                        Blast evalue (default: 1e-06)
  -t THREADS, --threads THREADS
                        number of threads to use, 0 for all (default: 0)
  --table TABLE         codon table (default: 11)
  -a [{nucleotide,protein} ...], --alignment [{nucleotide,protein} ...]
                        run alignment for each gene cluster (default: None)

```
## Example
Basic:
```
python pan-genome.py main -o examples/test/output -g examples/test/main/*.gff
python pan-genome.py add -c examples/test/output -g examples/test/add/*.gff
```
