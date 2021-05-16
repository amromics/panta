#!/bin/bash

if [ ! -f GCF_000240185.1_ASM24018v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
fi

if [ ! -f GCF_000364385.2_ASM36438v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/385/GCF_000364385.2_ASM36438v2/GCF_000364385.2_ASM36438v2_genomic.fna.gz
fi



for ac in SRR8607448 SRR8607464 SRR8607467 SRR8607459 SRR8607460 SRR8607456 SRR8607451 SRR8607450 SRR8607452 SRR8607471 SRR8607469 SRR8607462 SRR8607457 SRR8607455 SRR8607461 SRR8607458 SRR8607470 SRR8607449 SRR8607465 SRR8607454 SRR8607468 SRR8607463 SRR8607466 SRR8607453;do
    if [ ! -f ${ac}_2.fastq.gz ];then
        echo " Downloading Accession ${ac} ..."
        fasterq-dump --progress --split-3 ${ac} && gzip ${ac}_1.fastq && gzip ${ac}_2.fastq
    else
        echo " Accession ${ac} has been downloaded!"
    fi
done




