#!/bin/bash

if [ ! -f GCF_000240185.1_ASM24018v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
fi

if [ ! -f GCF_000364385.2_ASM36438v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/385/GCF_000364385.2_ASM36438v2/GCF_000364385.2_ASM36438v2_genomic.fna.gz
fi


if [ ! -f SRR8607448_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607448; gzip SRR8607448_1.fastq;gzip SRR8607448_2.fastq
fi
if [ ! -f SRR8607464_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607464; gzip SRR8607464_1.fastq;gzip SRR8607464_2.fastq
fi
if [ ! -f SRR8607467_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607467; gzip SRR8607467_1.fastq;gzip SRR8607467_2.fastq
fi
if [ ! -f SRR8607459_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607459; gzip SRR8607459_1.fastq;gzip SRR8607459_2.fastq
fi
if [ ! -f SRR8607460_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607460; gzip SRR8607460_1.fastq;gzip SRR8607460_2.fastq
fi
if [ ! -f SRR8607456_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607456; gzip SRR8607456_1.fastq;gzip SRR8607456_2.fastq
fi
if [ ! -f SRR8607451_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607451; gzip SRR8607451_1.fastq;gzip SRR8607451_2.fastq
fi
if [ ! -f SRR8607450_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607450; gzip SRR8607450_1.fastq;gzip SRR8607450_2.fastq
fi
if [ ! -f SRR8607452_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607452; gzip SRR8607452_1.fastq;gzip SRR8607452_2.fastq
fi
if [ ! -f SRR8607471_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607471; gzip SRR8607471_1.fastq;gzip SRR8607471_2.fastq
fi
if [ ! -f SRR8607469_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607469; gzip SRR8607469_1.fastq;gzip SRR8607469_2.fastq
fi
if [ ! -f SRR8607462_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607462; gzip SRR8607462_1.fastq;gzip SRR8607462_2.fastq
fi
if [ ! -f SRR8607457_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607457; gzip SRR8607457_1.fastq;gzip SRR8607457_2.fastq
fi
if [ ! -f SRR8607455_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607455; gzip SRR8607455_1.fastq;gzip SRR8607455_2.fastq
fi
if [ ! -f SRR8607461_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607461; gzip SRR8607461_1.fastq;gzip SRR8607461_2.fastq
fi
if [ ! -f SRR8607458_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607458; gzip SRR8607458_1.fastq;gzip SRR8607458_2.fastq
fi
if [ ! -f SRR8607470_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607470; gzip SRR8607470_1.fastq;gzip SRR8607470_2.fastq
fi
if [ ! -f SRR8607449_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607449; gzip SRR8607449_1.fastq;gzip SRR8607449_2.fastq
fi
if [ ! -f SRR8607465_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607465; gzip SRR8607465_1.fastq;gzip SRR8607465_2.fastq
fi
if [ ! -f SRR8607454_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607454; gzip SRR8607454_1.fastq;gzip SRR8607454_2.fastq
fi
if [ ! -f SRR8607468_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607468; gzip SRR8607468_1.fastq;gzip SRR8607468_2.fastq
fi
if [ ! -f SRR8607463_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607463; gzip SRR8607463_1.fastq;gzip SRR8607463_2.fastq
fi
if [ ! -f SRR8607466_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607466; gzip SRR8607466_1.fastq;gzip SRR8607466_2.fastq
fi
if [ ! -f SRR8607453_1.fastq.gz ];then
  fastq-dump --split-3 SRR8607453; gzip SRR8607453_1.fastq;gzip SRR8607453_2.fastq
fi

