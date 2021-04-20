import os
import shutil
import re
import json
import gzip
import csv
import logging
from Bio import SeqIO
import pandas as pd
from pan_genome.utils import run_command

logger = logging.getLogger(__name__)


def parse_cluster_file(cd_hit_cluster_file):
    """
    Parse cd-hit .clstr file

    Parameters
    -------
    -------
    """
    
    clusters = {}
    with open(cd_hit_cluster_file, 'r') as fh:
        for line in fh:
            if re.match(r"^>", line) != None:
                cluster_name = re.findall(r'^>(.+)$', line)
                cluster_name = cluster_name[0]
                clusters[cluster_name] = {}
                clusters[cluster_name]['gene_names'] = []
            else:
                result = re.findall(r'[\d]+\t[\w]+, >(.+)\.\.\. (.+)$', line)
                if len(result) == 1:
                    gene_name = result[0][0]
                    identity = result[0][1]
                    if identity == '*':
                        clusters[cluster_name]['representative'] = gene_name
                    else:
                        clusters[cluster_name]['gene_names'].append(gene_name)
    # convert to a simple dictionary
    clusters_new ={}
    for cluster_name in clusters:
        clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']

    return clusters_new


def run_cd_hit_iterative(report, threads, timing_log):
    """
    Run CD-HIT iteratively

    Parameters
    -------
    -------
    """
    temp_dir = report['temp_dir']
    combined_faa_file = report['combined_faa_file']

    remain_faa_file = os.path.join(temp_dir, 'remain.faa')
    shutil.copyfile(combined_faa_file, remain_faa_file)
    report['remain_faa_file'] = remain_faa_file

    cd_hit_cluster_fasta = os.path.join(temp_dir, 'cluster')
    cd_hit_cluster_file = os.path.join(temp_dir, 'cluster.clstr')

    excluded_cluster = []

    greater_than_or_equal = True
    number_of_samples = len(report['samples'])
    
    lower = 0.98
    step = 0.005
    percent_match = 1
    while percent_match >= lower:
        cmd = f'cd-hit -i {remain_faa_file} -o {cd_hit_cluster_fasta} -s {percent_match} -c {percent_match} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running cd-hit')
        
        clusters = parse_cluster_file(cd_hit_cluster_file)

        if percent_match != 1:
            greater_than_or_equal = False
        full_cluster_gene_names =[]
        for cluster_represent in clusters:
            other_genes = clusters[cluster_represent]
            this_cluster = []
            if greater_than_or_equal == True and len(other_genes) >= number_of_samples -1:
                this_cluster.append(cluster_represent)
                this_cluster.extend(other_genes)
            if greater_than_or_equal == False and len(other_genes) == number_of_samples -1:
                this_cluster.append(cluster_represent)
                this_cluster.extend(other_genes)
            if len(this_cluster) != 0:
                full_cluster_gene_names.extend(this_cluster)
                excluded_cluster.append(this_cluster)

        cluster_filtered_faa_file = remain_faa_file + '.filtered'
        with open(cluster_filtered_faa_file, 'w') as fh:
            for seq_record in SeqIO.parse(remain_faa_file, "fasta"):
                if seq_record.id not in full_cluster_gene_names:
                    SeqIO.write(seq_record, fh, 'fasta')
        shutil.move(cluster_filtered_faa_file, remain_faa_file)
        percent_match -= step

    cmd = f'cd-hit -i {remain_faa_file} -o {cd_hit_cluster_fasta} -s {lower} -c {lower} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit')
    clusters = parse_cluster_file(cd_hit_cluster_file)

    report['cd_hit_cluster_fasta'] = cd_hit_cluster_fasta
    report['cd_hit_cluster_file'] = cd_hit_cluster_file
    report['excluded_cluster'] = excluded_cluster
    report['cd_hit_cluster'] = clusters
    return report


def chunk_fasta_file(fasta_file, out_dir):
    """
    Take a fasta file and chunk it up into smaller files with the length of 200000

    Parameters
    -------
    -------
    """
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)
    else:
        os.makedirs(out_dir)
    
    chunked_file_list = []
    chunk_number = 0
    current_chunk_length = 0
    chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
    chunked_fh = open(chunked_file, 'w')
    chunked_file_list.append(chunked_file)
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if current_chunk_length > 200000:
            chunked_fh.close()
            chunk_number += 1
            current_chunk_length = 0
            chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
            chunked_file_list.append(chunked_file)
            chunked_fh = open(chunked_file, 'w')
            SeqIO.write(seq_record, chunked_fh, 'fasta')
        chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
        SeqIO.write(seq_record, chunked_fh, 'fasta')
        current_chunk_length += len(seq_record.seq)
    chunked_fh.close()
    return chunked_file_list


def all_against_all_blast(report, threads, timing_log):
    """
    Run all against all blast in parallel

    Parameters
    -------
    -------
    """
    out_dir = os.path.join(report['temp_dir'], 'blast')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    fasta_file = report['cd_hit_cluster_fasta']
    blast_db = os.path.join(out_dir, 'output_contigs')
    cmd = f"makeblastdb -in {fasta_file} -dbtype prot -out {blast_db} -logfile /dev/null"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running makeblastdb')
    
    # chunk fasta file
    chunk_dir = os.path.join(out_dir, 'chunk_files')
    chunked_file_list = chunk_fasta_file(fasta_file, chunk_dir)

    # run parallel all-against-all blast
    blast_cmds_file = os.path.join(out_dir,"blast_cmds.txt")
    blast_output_file_list = []
    with open(blast_cmds_file,'w') as fh:
        for chunked_file in chunked_file_list:
            blast_output_file = os.path.splitext(chunked_file)[0] + '.out'
            blast_output_file_list.append(blast_output_file)
            cmd = f"blastp -query {chunked_file} -db {blast_db} -evalue 1E-6 -num_threads {threads} -outfmt 6 -max_target_seqs 2000 " + "| awk '{ if ($3 >= 98) print $0;}' 2> /dev/null 1> " + blast_output_file
            fh.write(cmd + '\n')
    cmd = f"parallel --bar -j {threads} -a {blast_cmds_file}"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parallel all-against-all blast')

    # combining blast results
    blast_result_file = os.path.join(report['temp_dir'], 'blast_results')
    if os.path.isfile(blast_result_file):
        os.remove(blast_result_file)
    for blast_output_file in blast_output_file_list:
        os.system(f'cat {blast_output_file} >> {blast_result_file}')

    report['blast_result_file'] = blast_result_file
    return report


def cluster_with_mcl(report, threads, timing_log):
    """
    Take blast results and outputs clustered results

    Parameters
    -------
    -------
    """
    out_dir = report['temp_dir']
    blast_results = report['blast_result_file']
    output_mcl_file = os.path.join(out_dir, 'uninflated_mcl_clusters')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_results} 2> /dev/null | mcl - --abc -I 1.5 -o {output_mcl_file} > /dev/null 2>&1"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    report['uninflated_mcl_clusters'] = output_mcl_file
    return report


def reinflate_clusters(report):
    """
    Take the clusters file from cd-hit and use it to reinflate the output of MCL

    Parameters
    -------
    -------
    """
    inflated_clusters = []
    clusters = report['cd_hit_cluster']

    # Inflate genes from cdhit which were sent to mcl
    mcl_file = report['uninflated_mcl_clusters']
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                if gene in clusters:
                    inflated_genes.extend(clusters[gene])
                    del clusters[gene]
            inflated_clusters.append(inflated_genes)

    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(clusters[gene])
        inflated_clusters.append(inflated_genes)

    # Add clusters which were excluded
    excluded_cluster = report['excluded_cluster']
    for cluster in excluded_cluster:
        inflated_clusters.append(cluster)

    report['inflated_unsplit_clusters'] = inflated_clusters
    return report

