import os
import shutil
import re
import json
import gzip
import csv
import logging
from Bio import SeqIO
import pandas as pd
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def run_cd_hit_2d(report, threads, timing_log):
    temp_dir = report['temp_dir']
    combined_faa_file = report['combined_faa_file']
    representative_fasta = report['representative_fasta']

    not_match_fasta_file = os.path.join(temp_dir, 'cd_hit_2d')
    cd_hit_cluster_file = os.path.join(temp_dir, 'cd_hit_2d.clstr')
    
    persent = 0.98
    cmd = f'cd-hit-2d -i {representative_fasta} -i2 {combined_faa_file} -o {not_match_fasta_file} -s {persent} -c {persent} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit')

    clusters = parse_cluster_file(cd_hit_cluster_file)

    report['not_match_fasta_file'] = not_match_fasta_file
    report['cd_hit_cluster_file'] = cd_hit_cluster_file
    report['cd_hit_2d_cluster'] = clusters
    return report

def filter_fasta(blast_result, fasta_file, out_dir):
    with open(blast_result, 'r') as fh:
        ls = []
        for line in fh:
            line = line.rstrip()
            cells = line.split('\t')
            ls.append(cells[0])

    blast_remain_fasta = os.path.join(out_dir, 'blast_remain_fasta')
    with open(blast_remain_fasta, 'w') as fh:
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            if seq_record.id in ls:
                continue
            SeqIO.write(seq_record, fh, 'fasta')

    return blast_remain_fasta

def reinflate_clusters(old_clusters, cd_hit_2d_cluster, blast_1_result_file, mcl_clusters):
    
    blast_dataframe = pd.read_csv(blast_1_result_file, sep='\t', header = None)
    dictionary = {}
    min_value = 100
    previous = None
    for i, row in blast_dataframe.iterrows():
        new = row[0]
        repre = row[1]
        value = row[10]
        if new != previous:
            min_value = 100
            previous = new
        if value < min_value:
            dictionary[new] = repre
            min_value = value
    for new in dictionary:
        repre = dictionary[new]
        cd_hit_2d_cluster[repre].append(new)

    for cluster in old_clusters:
        for g in cluster:
            if g in cd_hit_2d_cluster:
                cluster.extend(cd_hit_2d_cluster[g])
    
    with open(mcl_clusters, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            genes = line.split('\t')
            old_clusters.append(genes)

    return old_clusters
    