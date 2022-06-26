#!/bin/bash

sudo docker build -t panta \
    --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) \
    https://github.com/amromics/amromics.git#progressive:amromics/pan-genome 


# docker run --rm -v /home/ted/pan-genome/fasta:/data panta \
#    sh -c 'python3 pan-genome.py main -b /data/main/*.fna -o /data/test -t 4'