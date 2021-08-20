

# /usr/bin/time -v satv_launch -d /home/ntanh1999/pan-genome/saturnv/Sp25 -c 8 -ann prokka -m lazy -a usearch -i 50 -ip 100 1>Sp25.log 2>&1

# /usr/bin/time -v satv_search-pangenome-lazy -g genomes_to_analyze.txt -c 8 -i 50 

import os
from glob import glob

gff_files = glob('/home/ntanh1999/pan-genome/amromics/amromics/pan-genome/data/Sp200/*.gff')

for gff_file in gff_files:
    sample_id = gff_file.split('/')[-1]
    sample_id = sample_id.split('.', 1)[0]
    # print(sample_id)
    # os.system(f'mkdir /home/ntanh1999/pan-genome/saturnv/Sp200/{sample_id}')
    os.system(f'cp {gff_file} /home/ntanh1999/pan-genome/saturnv/Sp200/')


./install -d /home/ntanh1999/pan-genome/saturnv -m core
