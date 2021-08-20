import os
import glob
from Bio import SeqIO
import json

a = os.listdir('/home/ntanh1999/pan-genome/amromics/amromics/pan-genome/tests/Sp616/616/samples')

gene_annotation = json.load(open('/home/ntanh1999/pan-genome/amromics/amromics/pan-genome/tests/Sp616/616/gene_annotation.json', 'r'))

with open('/home/ntanh1999/pan-genome/PanDelos/data/Sp25.faa', 'w') as fh:
    for sample_id in a[:25]:
        faa_file = os.path.join('/home/ntanh1999/pan-genome/amromics/amromics/pan-genome/tests/Sp616/616/samples',sample_id, sample_id + '.faa')
        for seq_record in SeqIO.parse(faa_file, "fasta"):
            seq_id = seq_record.id
            protein = seq_record.seq
            product = gene_annotation[seq_id]['product']
            fh.write(sample_id + '\t' + seq_id + '\t' + product + '\n')
            fh.write(str(protein) + '\n')



# /usr/bin/time -v bash pandelos.sh /home/ntanh1999/pan-genome/PanDelos/data/Sp25.faa /home/ntanh1999/pan-genome/PanDelos/data/Sp50  1>Sp25.log 2>&1

