import os
from glob import glob

# gff_files = glob('/home/ntanh1999/pan-genome/amromics/amromics/pan-genome/data/Sp616/*.gff')

# /usr/bin/time -v PIRATE -i /home/ntanh1999/amromics/amromics/pan-genome/data/Sp616 -o Sp616_pirate_616 -s "95,96,97,98" -t 8 -z 2 1>> pirate.log 2>> pirate.log

# /usr/bin/time -v PIRATE -i /home/ntanh1999/amromics/amromics/pan-genome/data/Sp400 -o Sp616_pirate_400 -s "95,96,97,98" -t 8 -z 2 1>> pirate.log 2>> pirate.log

# for gff_file in gff_files[0:25]:
#     os.system(f'cp {gff_file} /home/ntanh1999/pan-genome/amromics/amromics/pan-genome/data/Sp25/')


gff_files = glob('/home/ntanh1999/data/Sp616/*.fasta')

# /usr/bin/time -v PIRATE -i /home/ntanh1999/pan-genome/amromics/amromics/pan-genome/data/Sp25 -o Sp25_pirate_1 -s "95,96,97,98" -t 8 -z 2 1> Sp25_pirate.log 2>&1

# /usr/bin/time -v PIRATE -i /home/ntanh1999/amromics/amromics/pan-genome/data/Sp400 -o Sp616_pirate_400 -s "95,96,97,98" -t 8 -z 2 1>> pirate.log 2>> pirate.log

for gff_file in gff_files[0:25]:
    os.system(f'cp {gff_file} /home/ntanh1999/pan-genome/saturnv/Sp25/')