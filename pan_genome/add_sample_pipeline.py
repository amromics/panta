import os
import logging
from datetime import datetime
import pandas as pd
import pan_genome.utils as utils

logger = logging.getLogger(__name__)

def run_cd_hit_2d(database_1, database_2, out_dir, threads=4):
    starttime = datetime.now()

    not_match_fasta = os.path.join(out_dir, 'cd-hit-2d.fasta')
    cd_hit_cluster_file = not_match_fasta + '.clstr'
    
    cmd = f'cd-hit-2d -i {database_1} -i2 {database_2} -o {not_match_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = os.system(cmd)
    if ret != 0:
        raise Exception('Error running cd-hit-2d')

    clusters = utils.parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT-2D with 98% identity -- time taken {str(elapsed)}')
    return not_match_fasta, clusters


def add_gene_cd_hit_2d(old_clusters, cd_hit_2d_clusters):
    starttime = datetime.now()
    gene_to_cluster = {gene:index for index,genes in enumerate(old_clusters) for gene in genes}

    # add gene matched in CD-HIT-2D step
    for old in cd_hit_2d_clusters:
        cluster_index = gene_to_cluster[old]
        for gene in cd_hit_2d_clusters[old]:
            old_clusters[cluster_index].append(gene)
    
    elapsed = datetime.now() - starttime
    logging.info(f'Add new gene to clusters -- time taken {str(elapsed)}')
    return gene_to_cluster, old_clusters


def add_gene_blast(old_clusters, gene_to_cluster, unmatched_clusters, blast_result, fasta_file, out_dir, gene_annotation, identity, LD, AL, AS):
    starttime = datetime.now()
    
    blast_dataframe = pd.read_csv(blast_result, sep='\t', header = None)
    match_dict = {}
    min_evalue = 1000
    previous = None
    for i, row in blast_dataframe.iterrows():
        new = row[0]
        old = row[1]
        
        qlen = gene_annotation[new][2]
        slen = gene_annotation[old][2]
        
        pident = float(row[2]) / 100
        
        short_seq = min(qlen, slen)
        long_seq = max(qlen, slen)
        len_diff = short_seq / long_seq
        
        alignment_length = int(row[3]) * 3
        align_short = alignment_length / short_seq
        align_long = alignment_length / long_seq
            
        if pident <= identity or len_diff <= LD or align_short <= AS or align_long <= AL:
            continue
        
        evalue = row[10]
        if new != previous:
            min_evalue = 1000
            previous = new
        if evalue < min_evalue:
            match_dict[new] = old
            min_evalue = evalue
    
    for new in match_dict:
        old = match_dict[new]
        cluster_index = gene_to_cluster[old]
        old_clusters[cluster_index].append(new)
        old_clusters[cluster_index].extend(unmatched_clusters[new])
    

    blast_remain_fasta = os.path.join(out_dir, 'blast_remain_fasta')
    utils.create_fasta_exclude(fasta_file=fasta_file, exclude_list=match_dict.keys(), output_file=blast_remain_fasta)

    elapsed = datetime.now() - starttime
    logging.info(f'Add new gene to clusters -- time taken {str(elapsed)}')
    return blast_remain_fasta, old_clusters



def add_new_clusters(old_clusters, unmatched_clusters, mcl_file):
    starttime = datetime.now()

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
    logging.info(f'Add new clusters -- time taken {str(elapsed)}')
    return old_clusters
