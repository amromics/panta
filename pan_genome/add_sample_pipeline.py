import os
import logging
from datetime import datetime
import pandas as pd
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def run_cd_hit_2d(database_1, database_2, out_dir, identity=95, threads=4, timing_log=None):
    starttime = datetime.now()

    not_match_fasta = os.path.join(out_dir, 'cd_hit_2d')
    cd_hit_cluster_file = os.path.join(out_dir, 'cd_hit_2d.clstr')
    
    persent = identity / 100
    cmd = f'cd-hit-2d -i {database_1} -i2 {database_2} -o {not_match_fasta} -s {persent} -c {persent} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit')

    clusters = parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT-2D with {identity}% identity -- time taken {str(elapsed)}')
    return not_match_fasta, clusters

def filter_fasta(blast_result, fasta_file, out_dir):
    starttime = datetime.now()

    with open(blast_result, 'r') as fh:
        ls = []
        for line in fh:
            line = line.rstrip()
            cells = line.split('\t')
            ls.append(cells[0])
        ls = set(ls)
    print(f'Number of filterd sequences: {len(ls)}')
    blast_remain_fasta = os.path.join(out_dir, 'blast_remain_fasta')
    create_fasta_exclude(fasta_file=fasta_file, exclude_list=ls, output_file=blast_remain_fasta)

    elapsed = datetime.now() - starttime
    logging.info(f'Filter fasta -- time taken {str(elapsed)}')
    return blast_remain_fasta

def reinflate_clusters(old_clusters, cd_hit_2d_clusters, blast_1_result_file, mcl_clusters):
    starttime = datetime.now()

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
        cd_hit_2d_clusters[repre].append(new)

    for cluster in old_clusters:
        for g in cluster:
            if g in cd_hit_2d_clusters:
                cluster.extend(cd_hit_2d_clusters[g])
    
    with open(mcl_clusters, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            genes = line.split('\t')
            old_clusters.append(genes)

    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate clusters -- time taken {str(elapsed)}')
    return old_clusters
    