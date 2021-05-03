import os
import shutil
import logging
import subprocess
from datetime import datetime
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def exclude_full_clusters(cd_hit_cluster_file, remain_faa_file, number_of_samples, excluded_cluster):
    
    clusters = parse_cluster_file(cd_hit_cluster_file)

    full_cluster_gene_names =[]
    for cluster_represent in clusters:
        other_genes = clusters[cluster_represent]
        if len(other_genes) >= number_of_samples -1:
            this_cluster = [cluster_represent] + other_genes
            full_cluster_gene_names.extend(this_cluster)
            excluded_cluster.append(this_cluster)
    cluster_filtered_faa_file = remain_faa_file + '.filtered'
    full_cluster_gene_names = set(full_cluster_gene_names)
    create_fasta_exclude(
        fasta_file=remain_faa_file, 
        exclude_list=full_cluster_gene_names, 
        output_file=cluster_filtered_faa_file
        )
    shutil.move(cluster_filtered_faa_file, remain_faa_file)


def run_cd_hit_iterative(combined_faa_file, samples, out_dir, threads=4, timing_log=None):
    statime = datetime.now()
    starttime = datetime.now()
    remain_faa_file = os.path.join(out_dir, 'remain.faa')
    shutil.copyfile(combined_faa_file, remain_faa_file)
    cd_hit_represent_fasta = os.path.join(out_dir, 'cluster')
    cd_hit_cluster_file = os.path.join(out_dir, 'cluster.clstr')
    excluded_cluster = []
    number_of_samples = len(samples)

    # run first time with persent = 100%
    cmd = f'cd-hit -i {remain_faa_file} -o {cd_hit_represent_fasta} -s 1 -c 1 -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit')
    exclude_full_clusters(cd_hit_cluster_file, remain_faa_file, number_of_samples, excluded_cluster)
    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 100% identity -- time taken {str(elapsed)}')
    
    # run iteratively
    lower = 0.98
    step = 0.005
    percent_match = 0.99
    while percent_match >= lower:
        starttime = datetime.now()
        cmd = f'cd-hit -i {remain_faa_file} -o {cd_hit_represent_fasta} -s {percent_match} -c {percent_match} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running cd-hit')
        exclude_full_clusters(cd_hit_cluster_file, remain_faa_file, number_of_samples, excluded_cluster)
        elapsed = datetime.now() - starttime
        logging.info(f'Run CD-HIT with {percent_match * 100}% identity -- time taken {str(elapsed)}')
        percent_match -= step

    # run last time without excluding
    starttime = datetime.now()
    cmd = f'cd-hit -i {remain_faa_file} -o {cd_hit_represent_fasta} -s {lower} -c {lower} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit')
    cd_hit_clusters = parse_cluster_file(cd_hit_cluster_file)
    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with {lower * 100}% identity -- time taken {str(elapsed)}')
    elapsed = datetime.now() - statime
    logging.info(f'Run CD-HIT iteratively -- time taken {str(elapsed)}')
    return remain_faa_file, cd_hit_represent_fasta, cd_hit_cluster_file, excluded_cluster, cd_hit_clusters


def all_against_all_blast(database_fasta, query_fasta, out_dir, threads=4, timing_log=None):
    starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    blast_db = os.path.join(out_dir, 'output_contigs')
    cmd = f"makeblastdb -in {database_fasta} -dbtype prot -out {blast_db} -logfile /dev/null"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running makeblastdb')
    
    # chunk fasta file
    chunk_dir = os.path.join(out_dir, 'chunk_files')
    chunked_file_list = chunk_fasta_file(query_fasta, chunk_dir)

    # run parallel all-against-all blast
    blast_cmds_file = os.path.join(out_dir,"blast_cmds.txt")
    blast_output_file_list = []
    with open(blast_cmds_file,'w') as fh:
        for chunked_file in chunked_file_list:
            blast_output_file = os.path.splitext(chunked_file)[0] + '.out'
            blast_output_file_list.append(blast_output_file)
            cmd = f"blastp -query {chunked_file} -db {blast_db} -evalue 1E-6 -num_threads 1 -outfmt 6 -max_target_seqs 2000 " + "| awk '{ if ($3 >= 95) print $0;}' 2> /dev/null 1> " + blast_output_file
            fh.write(cmd + '\n')
    cmd = f"parallel --bar -j {threads} -a {blast_cmds_file}"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parallel all-against-all blast')

    # combining blast results
    blast_result = os.path.join(out_dir, 'blast_results')
    if os.path.isfile(blast_result):
        os.remove(blast_result)
    for blast_output_file in blast_output_file_list:
        os.system(f'cat {blast_output_file} >> {blast_result}')

    elapsed = datetime.now() - starttime
    logging.info(f'All-against-all BLASTP -- time taken {str(elapsed)}')
    return blast_result


def run_diamond(database_fasta, query_fasta, out_dir, threads=4, timing_log=None):
    starttime = datetime.now()

    # make diamond database
    diamond_db = os.path.join(out_dir, 'diamond_db')
    cmd = f'diamond makedb --in {database_fasta} -d {diamond_db} -p {threads} --quiet'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running diamond makedb')
    
    # run diamond blastp
    diamond_result = os.path.join(out_dir, 'diamond.tsv')
    cmd = f"diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} --ultra-sensitive --outfmt 6 " + "| awk '{ if ($3 >= 95) print $0;}' 2> /dev/null 1> " + diamond_result
    subprocess.call(cmd, shell=True)

    elapsed = datetime.now() - starttime
    logging.info(f'Protein alignment with Diamond -- time taken {str(elapsed)}')
    return diamond_result


def cluster_with_mcl(blast_results, out_dir, threads=4, timing_log=None):
    starttime = datetime.now()
    mcl_file = os.path.join(out_dir, 'mcl_clusters')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_results} 2> /dev/null | mcl - --abc -I 1.5 -o {mcl_file} > /dev/null 2>&1"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    elapsed = datetime.now() - starttime
    logging.info(f'Cluster with MCL -- time taken {str(elapsed)}')
    return mcl_file


def reinflate_clusters(cd_hit_clusters, mcl_file, excluded_cluster):
    starttime = datetime.now()
    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                if gene in cd_hit_clusters:
                    inflated_genes.extend(cd_hit_clusters[gene])
                    del cd_hit_clusters[gene]
            inflated_clusters.append(inflated_genes)

    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in cd_hit_clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(cd_hit_clusters[gene])
        inflated_clusters.append(inflated_genes)

    # Add clusters which were excluded
    for cluster in excluded_cluster:
        inflated_clusters.append(cluster)

    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate clusters -- time taken {str(elapsed)}')
    return inflated_clusters