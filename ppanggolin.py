# download dataset
# usage: python download.py [Klebsiella.tsv] [column] [output dir]
# column: 2 = assembly (fna.gz) 3 = annotation (gff.gz)  4 = protein (faa.gz)


import os
import sys
import pandas as pd
import re
import csv
from glob import glob
from datetime import datetime
import shutil
from Bio import SeqIO



def make_list_file(list_file,input_folder):

    f = open(list_file, "w")
    for gff in os.listdir(input_folder):
        file_path=os.path.join(input_folder,gff)
        if gff.endswith('.gff'):
            sample_id=gff.replace('.gff','')
            f.write(sample_id+'\t'+file_path+'\n')

    f.close()
def run_ppanggolin(file_list, out_dir_p,logfile):
    if not os.path.exists(out_dir_p):
        os.mkdir(out_dir_p)
    cmd=(f'/usr/bin/time -v ppanggolin workflow -f --anno {file_list} --output {out_dir_p} --cpu {threads}  1>> {logfile} 2>&1')
    os.system(cmd)
    return out_dir





out_dir = sys.argv[1]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)




threads=16

#running section

log_folder=os.path.join(out_dir,"logs")
if not os.path.exists(log_folder):
    os.mkdir(log_folder)


for Y in ['Sp600_batches/Sp600']:
    for Q in ['batch_gff_1','batch_gff_12','batch_gff_123','batch_gff_1234']:

        list_file=os.path.join(out_dir,"ppanggolin_"+Q+".list")
        if os.path.exists(list_file):
            os.remove(list_file)
        make_list_file(list_file,os.path.join(out_dir,Y+"/"+Q))

        run_dir_ppanggolin=os.path.join(out_dir,"ppanggolin_"+Q)
        if os.path.exists(run_dir_ppanggolin):
            shutil.rmtree(run_dir_ppanggolin)
        run_ppanggolin(list_file,run_dir_ppanggolin,os.path.join(log_folder,"ppanggolin_"+Q+".log"))
