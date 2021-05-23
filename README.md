# Pan-genome pipeline
## Installation
Depencencies:
- bedtools
- CD-HIT
- BLAST
- DIAMOND
- MCL
- pandas
## Usage
Activate conda enviroment:
```
source activate amromics-viz
```
### Main pipeline: run pan-genome analysis for the first time
```
usage: pan-genome.py main [-options] -o OUT_DIR gff_files [gff_files ...]

positional arguments:
  gff files             a.gff b.gff ... (*.gff)

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory (default: None)
  -s, --dont-split      dont split paralog clusters (default: False)
  -d, --diamond         use Diamond for all-agaist-all alignment instead of
                        Blastp (default: False)
  -i IDENTITY, --identity IDENTITY
                        minimum percentage identity (default: 95)
  -f, --fasta           fasta files are seperated from gff files. (fasta file
                        must have the same name, be in the same folder of
                        coresponding gff file, and have one of following
                        extension: .fasta .fna .fnn) (default: False)
  -t THREADS, --threads THREADS
                        number of threads to use, 0 for all (default: 0)
  -w, --over-write      overwrite the previous results (default: False)
```
### Add pipeline: add sample into previous collection
```
usage: pan-genome.py add [-options] -c COLLECTION_DIR gff_files [gff_files ...]

positional arguments:
  gff_files             a.gff b.gff ... (*.gff)

optional arguments:
  -h, --help            show this help message and exit
  -c COLLECTION_DIR, --collection-dir COLLECTION_DIR
                        directory of previous collection (default: None)
  -s, --dont-split      dont split paralog clusters (default: False)
  -d, --diamond         use Diamond for all-agaist-all alignment instead of
                        Blastp (default: False)
  -i IDENTITY, --identity IDENTITY
                        minimum percentage identity (default: 95)
  -f, --fasta           fasta files are seperated from gff files. (fasta file
                        must have the same name, be in the same folder of
                        coresponding gff file, and have one of following
                        extensions: .fasta .fna .fnn) (default: False)
  -t THREADS, --threads THREADS
                        number of threads to use, 0 for all (default: 0)
```
## Example
Basic:
```
python pan-genome.py main -o tests/main -t 4 data/Kp26/main/*.gff
python pan-genome.py add -c tests/main -t 4 data/Kp26/add/*.gff
```