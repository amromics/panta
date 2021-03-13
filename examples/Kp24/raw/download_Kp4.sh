#!/bin/bash

if [ ! -f GCF_000240185.1_ASM24018v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
fi

if [ ! -f SRR8607448_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607448; gzip SRR8607448_1.fastq;gzip SRR8607448_2.fastq
fi
if [ ! -f SRR8607451_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607451; gzip SRR8607451_1.fastq;gzip SRR8607451_2.fastq
fi
if [ ! -f SRR8607450_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607450; gzip SRR8607450_1.fastq;gzip SRR8607450_2.fastq
fi
if [ ! -f SRR8607449_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607449; gzip SRR8607449_1.fastq;gzip SRR8607449_2.fastq
fi
