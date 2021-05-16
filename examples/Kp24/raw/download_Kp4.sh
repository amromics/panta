#!/bin/bash

if [ ! -f GCF_000240185.1_ASM24018v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
fi


for ac in SRR8607448 SRR8607451 SRR8607450 SRR8607449;do
    if [ ! -f ${ac}_2.fastq.gz ];then
        echo " Downloading Accession ${ac} ..."
        fasterq-dump --progress --split-3 ${ac} && gzip ${ac}_1.fastq && gzip ${ac}_2.fastq
    else
        echo " Accession ${ac} has been downloaded!"
    fi
done

