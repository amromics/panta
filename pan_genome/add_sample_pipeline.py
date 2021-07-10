import os
import logging
from datetime import datetime
import pandas as pd
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def run_cd_hit_2d(database_1, database_2, out_dir, threads=4, timing_log=None):
    starttime = datetime.now()

    not_match_fasta = os.path.join(out_dir, 'cd-hit-2d.fasta')
    cd_hit_cluster_file = not_match_fasta + '.clstr'
    
    cmd = f'cd-hit-2d -i {database_1} -i2 {database_2} -o {not_match_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit-2d')

    clusters = parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT-2D with 98% identity -- time taken {str(elapsed)}')
    return not_match_fasta, clusters

# def filter_fasta(blast_result, fasta_file, out_dir):
#     starttime = datetime.now()

#     with open(blast_result, 'r') as fh:
#         ls = []
#         for line in fh:
#             line = line.rstrip()
#             cells = line.split('\t')
#             ls.append(cells[0])
#         ls = set(ls)
#     print(f'Number of filterd sequences: {len(ls)}')
#     blast_remain_fasta = os.path.join(out_dir, 'blast_remain_fasta')
#     create_fasta_exclude(fasta_file=fasta_file, exclude_list=ls, output_file=blast_remain_fasta)

#     elapsed = datetime.now() - starttime
#     logging.info(f'Filter fasta -- time taken {str(elapsed)}')
#     return blast_remain_fasta


def combine_blast_results(blast_1, blast_2, blast_3, outdir):
    combined_blast_results = os.path.join(outdir, 'combined_blast_results')

    os.system(f'cat {blast_1} {blast_2} {blast_3} > {combined_blast_results}')

    return combined_blast_results


def reinflate_clusters(old_clusters, cd_hit_2d_clusters, not_match_clusters, mcl_file):
    starttime = datetime.now()

    # clusters for next run
    new_clusters = {}
    new_clusters.update(old_clusters)
    for gene in new_clusters:
        new_clusters[gene].extend(cd_hit_2d_clusters[gene])
    new_clusters.update(not_match_clusters)

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