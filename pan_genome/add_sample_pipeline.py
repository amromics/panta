import os
import logging
import copy
from datetime import datetime
import pandas as pd
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def run_cd_hit_2d(database_1, database_2, out_dir, threads=4):
    starttime = datetime.now()

    not_match_fasta = os.path.join(out_dir, 'cd-hit-2d.fasta')
    cd_hit_cluster_file = not_match_fasta + '.clstr'
    
    cmd = f'cd-hit-2d -i {database_1} -i2 {database_2} -o {not_match_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = os.system(cmd)
    if ret != 0:
        raise Exception('Error running cd-hit-2d')

    clusters = parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT-2D with 98% identity -- time taken {str(elapsed)}')
    return not_match_fasta, clusters


def combine_blast_results(blast_1, blast_2, blast_3, outdir):
    combined_blast_results = os.path.join(outdir, 'combined_blast_results')

    os.system(f'cat {blast_1} {blast_2} {blast_3} > {combined_blast_results}')

    return combined_blast_results


def combine_representative(new, old, out_dir):
    temp_file = os.path.join(out_dir, 'representative_temp')
    out_file = os.path.join(out_dir, 'representative.fasta')
    os.system(f'cat {old} {new} > {temp_file}')

    os.replace(temp_file, out_file)


def reinflate_clusters(old_clusters, cd_hit_2d_clusters, not_match_clusters, mcl_file):
    starttime = datetime.now()

    # clusters for next run
    new_clusters = copy.deepcopy(old_clusters)
    for gene in new_clusters:
        new_clusters[gene].extend(cd_hit_2d_clusters[gene])
    new_clusters.update(copy.deepcopy(not_match_clusters))

    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                if gene in old_clusters:
                    inflated_genes.extend(old_clusters[gene])
                    inflated_genes.extend(cd_hit_2d_clusters[gene])
                    del old_clusters[gene]
                    del cd_hit_2d_clusters[gene]
                if gene in not_match_clusters:
                    inflated_genes.extend(not_match_clusters[gene])
                    del not_match_clusters[gene]
            inflated_clusters.append(inflated_genes)
    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in old_clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(old_clusters[gene])
        inflated_genes.extend(cd_hit_2d_clusters[gene])
        inflated_clusters.append(inflated_genes)

    for gene in not_match_clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(not_match_clusters[gene])

    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate clusters -- time taken {str(elapsed)}')
    return inflated_clusters, new_clusters




def filter_fasta(blast_result, fasta_file, out_dir):
    starttime = datetime.now()

    with open(blast_result, 'r') as fh:
        ls = []
        for line in fh:
            line = line.rstrip()
            cells = line.split('\t')
            ls.append(cells[0])
        ls = set(ls)
    blast_remain_fasta = os.path.join(out_dir, 'blast_remain_fasta')
    create_fasta_exclude(fasta_file=fasta_file, exclude_list=ls, output_file=blast_remain_fasta)

    elapsed = datetime.now() - starttime
    logging.info(f'Filter fasta -- time taken {str(elapsed)}')
    return blast_remain_fasta



def reinflate_clusters_2(old_clusters, cd_hit_2d_clusters, blast_1_result, unmatched_clusters, mcl_file):
    starttime = datetime.now()
    
    gene_to_cluster_index = {gene:index for index,genes in enumerate(old_clusters) for gene in genes}

    for old in cd_hit_2d_clusters:
        cluster_index = gene_to_cluster_index[old]
        for gene in cd_hit_2d_clusters[old]:
            old_clusters[cluster_index].append(gene)

    blast_dataframe = pd.read_csv(blast_1_result, sep='\t', header = None)
    dictionary = {}
    min_value = 1000
    previous = None
    for i, row in blast_dataframe.iterrows():
        new = row[0]
        old = row[1]
        identity = row[2]
        value = row[10]
        if new != previous:
            min_value = 1000
            previous = new
        if value < min_value:
            dictionary[new] = old
            min_value = value
    for new in dictionary:
        old = dictionary[new]
        cluster_index = gene_to_cluster_index[old]
        old_clusters[cluster_index].append(new)
        old_clusters[cluster_index].extend(unmatched_clusters[new])

    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                inflated_genes.extend(unmatched_clusters[gene])
            old_clusters.append(inflated_genes)

    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate clusters -- time taken {str(elapsed)}')
    return old_clusters



def add_rep_fasta(old_rep_fasta, rep_clusters, new_seq_file, ourdir):
    rep_seq = set()
    for cluster in rep_clusters:
        for gene in cluster:
            rep_seq.add(gene)
    
    temp_fasta = os.path.join(ourdir, 'rep_temp.fasta')
    create_fasta_include(new_seq_file, rep_seq, temp_fasta)

    rep_fasta = rep_fasta = os.path.join(ourdir, 'rep.fasta')
    os.system(f'cat {old_rep_fasta} {temp_fasta} > {rep_fasta}')
    os.remove(temp_fasta)

    return temp_fasta