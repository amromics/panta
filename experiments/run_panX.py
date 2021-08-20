import os
from glob import glob

os.chdir('/home/ntanh1999')

# for filename in os.listdir('/home/ntanh1999/pan-genome-analysis/data/Kp616'):
#     cells = filename.split('.')
#     new_filename = cells[0] + '.gbk'
#     os.system(f'mv /home/ntanh1999/pan-genome-analysis/data/Kp616/{filename} /home/ntanh1999/pan-genome-analysis/data/Kp616/{new_filename}')


# cmd = '/usr/bin/time -v ./panX.py -fn ./data/Kp616 -sl Streptococcus_pneumoniae -dmdc -sitr -t 8  1> Kp616.log 2> Kp616.log'

# /usr/bin/time -v ./panX.py -fn ./data/Sp25 -sl Streptococcus_pneumoniae -sitr -t 8  1> Sp25.log 2>&1


gff_files = glob('/home/ntanh1999/pan-genome/panX/data/Sp616/input_GenBank/*.gbk')

for gff_file in gff_files[0:25]:
    os.system(f'cp {gff_file} /home/ntanh1999/pan-genome/panX/data/Sp25')