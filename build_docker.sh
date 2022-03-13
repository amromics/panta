#!/bin/bash

sudo docker build https://github.com/amromics/amromics.git#progressive:amromics/pan-genome -t panta


# docker run --rm -v /home/ted/pan-genome/fasta:/data panta sh -c 'python3 pan-genome.py main -b /data/main/*.fna -o /data/test -t 4'