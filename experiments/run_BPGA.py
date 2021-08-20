import os
import glob

a = os.listdir('/home/ntanh1999/pan-genome/BPGA/Sp200')

for fasta in a:
    sample_id = fasta.split('.')[0]
    os.system(f'mv /home/ntanh1999/pan-genome/BPGA/Sp200/{fasta} /home/ntanh1999/pan-genome/BPGA/Sp200/{sample_id}.fasta')