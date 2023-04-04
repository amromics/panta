import os
import sys
import pandas as pd
import re
import csv
import json
from glob import glob
from datetime import datetime
import shutil
import gzip
fna_dir = sys.argv[1]
path_out_gff = sys.argv[2]
path_out_pangenome=path_out_gff+'/pan'
#fna_dir='/mnt/data/data/amromics/pantadata/Hoan/data'
threads=16
#path_out_gff='/mnt/data/data/amromics/pantadata/Hoan'
#path_out_pangenome='/mnt/data/data/amromics/pantadata/Hoan/pan'
# for fna_file in os.listdir(fna_dir):
#     sample_id=fna_file.rsplit('.', 1)[0]
#     if not os.path.exists(os.path.join(path_out_gff,sample_id+".gff")):
#         cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
#         cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=sample_id, path=path_out_gff)
#         cmd += ' ' + os.path.join(fna_dir,fna_file)
#         print(cmd)
#         os.system(cmd)
cmd=f'python pan-genome.py main -o {path_out_pangenome} -g {path_out_gff}/*.gff'
os.system(cmd)
#create samples.tsv
dict_samples=json.load(open(os.path.join(path_out_pangenome,'samples.json')))
#print(dict_samples)
map_stringid_to_numberic_id={}
for i in range(len(dict_samples)):
    #print(dict_samples[i]['id'])
    map_stringid_to_numberic_id[str(dict_samples[i]['id'])]=i
print(map_stringid_to_numberic_id)
dict_clusters=json.load(open(os.path.join(path_out_pangenome,'clusters.json')))
map_geneid_to_numberic_cid={}
cindex=0
for k in dict_clusters.keys():
    #print(k)
    #print(dict_clusters[k])

    map_geneid_to_numberic_cid[k] =cindex

    for g in dict_clusters[k]:
        map_geneid_to_numberic_cid[g]=cindex
    cindex=cindex+1
    #map_stringid_to_numberic_id[dict_samples[i]['id']]=i
#print(map_geneid_to_numberic_cid)
#map gene to strain
df_annotations= pd.read_csv(os.path.join(path_out_pangenome,'gene_annotation.csv'))
print(df_annotations.head())
map_geneid_to_info={}
for index, row  in df_annotations.iterrows():
    strand='1'
    if row['strand']=='-':
        strand='-1'
    map_geneid_to_info[row['gene_id']]={'sample_id':str(row['sample_id']),'seq_id':row['seq_id'],'length':row['length'],'strand':strand}
#print(map_geneid_to_info)
#samples.tsv
with open(os.path.join(path_out_pangenome,'samples.tsv'), 'w') as sampletsv:
    #sampletsv.write('Name\tSampleID\n')
    for k in map_stringid_to_numberic_id.keys():
        sampletsv.write(f'{k}\t{map_stringid_to_numberic_id[k]}\n')
#gene_info.tsv
with open(os.path.join(path_out_pangenome,'gene_info.tsv'), 'w') as genetsv:
    #genetsv.write('GeneName\tSampleID\tclusterID\n')
    for k in map_geneid_to_info.keys():
        print(map_geneid_to_info[k]['sample_id'])
        genetsv.write(k+'@'+map_geneid_to_info[k]['strand']+'\t'+str(map_stringid_to_numberic_id[str(map_geneid_to_info[k]['sample_id'])])+'\t'+str(map_geneid_to_numberic_cid[k])+'\n')
with open(os.path.join(path_out_pangenome,'gene_position.tsv'), 'w') as genepos_tsv,  open(os.path.join(path_out_pangenome,'gene_position.csv'), 'rt') as genepos_csv:
    lines = genepos_csv.readlines()

    #genepos_tsv.write('SampleID\tContigName\tGeneSequence\n')
    for line in lines:
        row=line[:-1].split(',')

        genepos_tsv.write(f'{map_stringid_to_numberic_id[row[0]]}\t{row[1]}\t')
        str_genes=''
        for i in range(2,len(row)):
            str_genes=str_genes+row[i]+'@'+map_geneid_to_info[row[i]]['strand']+';'
        genepos_tsv.write(str_genes[:-1]+'\n')
