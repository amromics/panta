import os
import logging
from datetime import datetime
import csv

import pan_genome.utils as utils

logger = logging.getLogger(__name__)

def run_cd_hit_2d(database_1, database_2, out_dir, threads, timing_log):
    """
    Run CD-HIT-2D to compare new sequences from new samples with the 
    representative sequences of previous clusters. Then, parse the 
    resulting cluster file.

    Parameters
    ----------
    database_1 : path
        fasta file of representative sequences
    database 2 : path
        fasta file of new sequences
    out_dir : path
        directory of the output files (temp_dir)
    threads : int
        number of threads
    timing_log : path
        path of time.log
    
    Returns
    -------
    notmatch_faa : path
        fasta file containing not-matched sequences
    clusters : dict
        dictionary of CD-HIT-2D clusters
        {representative seq : [other sequences ids]}
    """
    starttime = datetime.now()
    
    # run CD-HIT-2D
    notmatch_faa = os.path.join(out_dir, 'cd-hit-2d.fasta')
    cd_hit_cluster_file = notmatch_faa + '.clstr'
    cmd=(f'cd-hit-2d -i {database_1} -i2 {database_2} -o {notmatch_faa} '
         f'-s 0.98 -s2 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null')
    utils.run_command(cmd, timing_log)

    # parse the cluster file
    clusters = utils.parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(
        f'Run CD-HIT-2D with 98% identity -- time taken {str(elapsed)}')
    return notmatch_faa, clusters


def add_gene_cd_hit_2d(previous_clusters, cd_hit_2d_clusters):
    """
    Add CD-HIT-2D matched sequences into coressponding previous clusters.

    Parameters
    ----------
    previous_clusters : list of []
        list of previous clusters
    cd_hit_2d_clusters : dict
        dictionary of CD-HIT-2D clusters
        {representative seq : [other sequences ids]}
    """
    starttime = datetime.now()

    for cluster_index in cd_hit_2d_clusters:
        new_seqs = cd_hit_2d_clusters[cluster_index]
        cluster_index = int(cluster_index)
        previous_clusters[cluster_index].extend(new_seqs)

    elapsed = datetime.now() - starttime
    logging.info(f'Add new gene to clusters -- time taken {str(elapsed)}')


def add_gene_blast(previous_clusters, cd_hit_clusters, blast_result, 
                   fasta_file, out_dir):
    """
    New sequences are added into previous clusters by its best hit 
    (lowest e-value), along with the corresponding CD-HIT sequences. 
    The not-matched sequences are write into a new fasta file.

    Parameters
    ----------
    previous_clusters : list
        list of previous clusters (with new sequences added)
    cd_hit_clusters : dict
        dictionary of CD-HIT clusters
        {representative seq : [other sequences ids]}
    blast_result : path
        blast 1 result file
    fasta_file : path
        fasta file, which is compared with previous clusters by BLAST
        It will be filtered out matched sequences
    out_dir : path
        directory of filtered fasta output
    
    Returns
    -------
    remain_fasta : path
        fasta file contains remainning sequences
    """
    starttime = datetime.now()
    
    # Filter blast result to find best matched sequences
    csv_reader = csv.reader(open(blast_result, 'r'), delimiter='\t')
    match_dict = {}
    min_evalue = 1000
    previous = None
    for row in csv_reader:
        new = row[0]
        cluster_index = row[1]
        evalue = float(row[10])
        if new != previous:
            min_evalue = 1000
            previous = new
        if evalue < min_evalue:
            match_dict[new] = cluster_index
            min_evalue = evalue
    
    # add matched sequences into previous clusters
    for new in match_dict:
        cluster_index = match_dict[new]
        cluster_index = int(cluster_index)
        previous_clusters[cluster_index].append(new)
        previous_clusters[cluster_index].extend(cd_hit_clusters[new])
        # deleted added sequences, so the remainning sequences
        # can be added later
        del cd_hit_clusters[new]
    
    # write the remaining sequences in to new fasta file
    remain_fasta = os.path.join(out_dir, 'remain_fasta')
    utils.create_fasta_exclude(
        fasta_file_list=[fasta_file], 
        exclude_list=match_dict.keys(), 
        output_file=remain_fasta)

    elapsed = datetime.now() - starttime
    logging.info(f'Add new gene to clusters -- time taken {str(elapsed)}')
    return remain_fasta
