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
## Usage
Activate conda enviroment:
```
source activate amromics-viz
```
### Main pipeline: run pan-genome analysis for the first time
```
usage: pan-genome.py main [-h] -c COLLECTION_DIR [-s] [-d] [-i IDENTITY] [-f] [-t THREADS] gff_files [gff_files ...]

Main pipeline: run pan-genome analysis for the first time

positional arguments:
  gff_files             a.gff b.gff ... (*.gff)

optional arguments:
  -h, --help            show this help message and exit
  -c COLLECTION_DIR, --collection-dir COLLECTION_DIR
                        collection directory (default: None)
  -s, --dont-split      dont split paralog clusters (default: False)
  -d, --diamond         use Diamond for all-agaist-all alignment instead of Blastp (default: False)
  -i IDENTITY, --identity IDENTITY
                        minimum percentage identity (default: 95)
  -f, --fasta           fasta files are seperated from gff files. (fasta file must have the same name, be in the same folder of coresponding gff
                        file, and have one of following extension: .fasta .fna .fnn) (default: False)
  -t THREADS, --threads THREADS
                        number of threads to use, 0 for all (default: 0)
```
### Add pipeline: add sample into previous collection
```
usage: pan-genome.py add [-h] -c COLLECTION_DIR [-s] [-d] [-i IDENTITY] [-f] [-t THREADS] gff_files [gff_files ...]

Add pipeline: add sample into previous collection

positional arguments:
  gff_files             a.gff b.gff ... (*.gff)

optional arguments:
  -h, --help            show this help message and exit
  -c COLLECTION_DIR, --collection-dir COLLECTION_DIR
                        previous collection directory (default: None)
  -s, --dont-split      dont split paralog clusters (default: False)
  -d, --diamond         use Diamond for all-agaist-all alignment instead of Blastp (default: False)
  -i IDENTITY, --identity IDENTITY
                        minimum percentage identity (default: 95)
  -f, --fasta           fasta files are seperated from gff files. (fasta file must have the same name, be in the same folder of coresponding gff
                        file, and have one of following extension: .fasta .fna .fnn) (default: False)
  -t THREADS, --threads THREADS
                        number of threads to use, 0 for all (default: 0)
```
## Example
Basic:
```
python3 pan-genome.py main -c tests/Kp26 -t 8 data/Kp26/main/*.gff
python3 pan-genome.py add -c tests/Kp26 -t 8 data/Kp26/add/*.gff
```