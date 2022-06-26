import os
import logging
from datetime import datetime
import csv

import pan_genome.utils as utils

logger = logging.getLogger(__name__)

def run_cd_hit_2d(database_1, database_2, out_dir, threads, timing_log):
    starttime = datetime.now()

    not_match_fasta = os.path.join(out_dir, 'cd-hit-2d.fasta')
    cd_hit_cluster_file = not_match_fasta + '.clstr'
    
    cmd=(f'cd-hit-2d -i {database_1} -i2 {database_2} -o {not_match_fasta} '
         f'-s 0.98 -s2 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null')
    utils.run_command(cmd, timing_log)

    clusters = utils.parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(
        f'Run CD-HIT-2D with 98% identity -- time taken {str(elapsed)}')
    return not_match_fasta, clusters


def add_gene_cd_hit_2d(old_clusters, cd_hit_2d_clusters):
    starttime = datetime.now()

    for cluster_index in cd_hit_2d_clusters:
        new_seqs = cd_hit_2d_clusters[cluster_index]
        cluster_index = int(cluster_index)
        old_clusters[cluster_index].extend(new_seqs)

    elapsed = datetime.now() - starttime
    logging.info(f'Add new gene to clusters -- time taken {str(elapsed)}')
    return old_clusters


def add_gene_blast(old_clusters, unmatched_clusters, blast_result, 
                   fasta_file, out_dir, identity, LD, AL, AS):
    starttime = datetime.now()
    
    csv_reader = csv.reader(open(blast_result, 'r'), delimiter='\t')
    match_dict = {}
    min_evalue = 1000
    previous = None
    for row in csv_reader:
        new = row[0]
        cluster_index = row[1]
        
        qlen = int(row[12])
        slen = int(row[13])
        
        pident = float(row[2]) / 100
        
        short_seq = min(qlen, slen)
        long_seq = max(qlen, slen)
        len_diff = short_seq / long_seq
        
        alignment_length = int(row[3])
        align_short = alignment_length / short_seq
        align_long = alignment_length / long_seq
            
        if (pident <= identity 
            or len_diff <= LD 
            or align_short <= AS 
            or align_long <= AL):
            continue
        
        evalue = float(row[10])
        if new != previous:
            min_evalue = 1000
            previous = new
        if evalue < min_evalue:
            match_dict[new] = cluster_index
            min_evalue = evalue
    
    for new in match_dict:
        cluster_index = match_dict[new]
        cluster_index = int(cluster_index)
        old_clusters[cluster_index].append(new)
        old_clusters[cluster_index].extend(unmatched_clusters[new])
        del unmatched_clusters[new]
    
    blast_remain_fasta = os.path.join(out_dir, 'blast_remain_fasta')
    utils.create_fasta_exclude(
        fasta_file_list=[fasta_file], 
        exclude_list=match_dict.keys(), 
        output_file=blast_remain_fasta)

    elapsed = datetime.now() - starttime
    logging.info(f'Add new gene to clusters -- time taken {str(elapsed)}')
    return blast_remain_fasta, old_clusters



# def create_new_clusters(unmatched_clusters, mcl_file, gene_dictionary):
#     starttime = datetime.now()

#     new_clusters = []
#     new_represent_list = []
#     with open(mcl_file, 'r') as fh:
#         for line in fh:
#             inflated_genes = []
#             line = line.rstrip('\n')
#             genes = line.split('\t')
#             length_max = 0
#             representative = None
#             for gene in genes:
#                 length = gene_dictionary[gene][2]
#                 if length > length_max:
#                     representative = gene
#                     length_max = length
#                 inflated_genes.append(gene)
#                 inflated_genes.extend(unmatched_clusters[gene])
#             new_clusters.append(inflated_genes)
#             new_represent_list.append(representative)

#     elapsed = datetime.now() - starttime
#     logging.info(f'Add new clusters -- time taken {str(elapsed)}')
#     return new_clusters, new_represent_list
