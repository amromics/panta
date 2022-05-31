import os
import re
import logging
import shutil
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# import pyabpoa as pa
from datetime import datetime
import pan_genome.utils as utils
import multiprocessing
from functools import partial

def create_pro_file_for_each_cluster(samples, gene_to_cluster, collection_dir):
    starttime = datetime.now()
    clusters_dir = os.path.join(collection_dir, 'clusters')
    os.system(f'find {clusters_dir} -name "*.faa" -delete') # remove old sequence file

    for sample in samples:
        sample_id = sample['id']
        faa_file = os.path.join(collection_dir, 'samples', sample_id, sample_id + '.faa')
        with open(faa_file, 'r') as in_fh:
            last_seq_id = None
            seq_lines = []
            for line in in_fh:
                line = line.rstrip()
                result = re.match(r"^>(\S+)", line)
                if result != None:
                    seq_id = result.group(1)
                    if last_seq_id != None:
                        cluster_id = str(gene_to_cluster[last_seq_id])
                        cluster_dir = os.path.join(collection_dir, 'clusters', cluster_id)
                        if not os.path.exists(cluster_dir):
                            os.mkdir(cluster_dir)
                        with open(os.path.join(cluster_dir, cluster_id + '.faa'), 'a') as out_fh:
                            out_fh.write('>'+ last_seq_id + '\n')
                            out_fh.write('\n'.join(seq_lines) + '\n')
                    last_seq_id = seq_id
                    seq_lines = []
                else:
                    seq_lines.append(line)

            cluster_id = str(gene_to_cluster[last_seq_id])
            cluster_dir = os.path.join(collection_dir, 'clusters', cluster_id)
            if not os.path.exists(cluster_dir):
                os.mkdir(cluster_dir)
            with open(os.path.join(cluster_dir, cluster_id + '.faa'), 'a') as out_fh:
                out_fh.write('>'+ last_seq_id + '\n')
                out_fh.write('\n'.join(seq_lines) + '\n')

    elapsed = datetime.now() - starttime
    logging.info(f'Create protein sequence file for each gene cluster -- time taken {str(elapsed)}')


def create_poa(cluster_id, collection_dir, baseDir):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(collection_dir, 'clusters', cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + '.faa')

    matrix_file = os.path.join(baseDir, 'BLOSUM62.mtx')
    result_file = os.path.join(cluster_dir, cluster_id + '.result')
    cmd = f'abpoa {seq_file} -o {result_file} -r2 -t {matrix_file} -O 11,0 -E 1,0 -p -c 2> /dev/null'
    utils.run_command(cmd)
    cons_seq = ''
    msa_file = os.path.join(cluster_dir, cluster_id + '.msa')
    with open(result_file, 'r') as in_fh, open(msa_file, 'w') as out_fh:
        fasta_out = SeqIO.FastaIO.FastaWriter(out_fh, wrap=None)
        for seq_record in SeqIO.parse(in_fh, 'fasta'):
            if seq_record.id == 'Consensus_sequence':
                cons_seq = str(seq_record.seq)
            else:
                fasta_out.write_record(seq_record)
    # os.remove(result_file)

    return cons_seq

def create_poa_in_parallel(clusters_id_list, collection_dir, baseDir, out_dir, threads):
    starttime = datetime.now()

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(create_poa,collection_dir=collection_dir, baseDir=baseDir), clusters_id_list)

    representative_fasta = os.path.join(out_dir, 'reference_pangenome.fasta')
    with open(representative_fasta, 'w') as fh:
        for cluster_id, cons_seq in enumerate(results):
            cons_seq = re.sub('-', '', cons_seq) # remove dash character in consensus sequences
            new_record = SeqRecord(Seq(cons_seq), id = str(cluster_id), description = '')
            SeqIO.write(new_record, fh, 'fasta')

    elapsed = datetime.now() - starttime
    logging.info(f'Create multiple sequence alignment by abPOA -- time taken {str(elapsed)}')

    return representative_fasta


def add_poa(cluster_id, collection_dir, baseDir):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(collection_dir, 'clusters', cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + '.faa')
    msa_file = os.path.join(cluster_dir, cluster_id + '.msa')
    matrix_file = os.path.join(baseDir, 'BLOSUM62.mtx')
    result_file = os.path.join(cluster_dir, cluster_id + '.result')
    if os.path.isfile(seq_file):
        cmd = f'abpoa {seq_file} -i {msa_file} -o {result_file} -r2 -t {matrix_file} -O 11,0 -E 1,0 -p -c -m 1 2> /dev/null'
        utils.run_command(cmd)

    cons_seq = ''
    with open(result_file, 'r') as in_fh, open(msa_file, 'w') as out_fh:
        fasta_out = SeqIO.FastaIO.FastaWriter(out_fh, wrap=None)
        for seq_record in SeqIO.parse(in_fh, 'fasta'):
            if seq_record.id == 'Consensus_sequence':
                cons_seq = str(seq_record.seq)
            else:
                fasta_out.write_record(seq_record)
    # os.remove(result_file)

    return cons_seq

def add_poa_in_parallel(clusters_id_list, collection_dir, baseDir, threads):
    starttime = datetime.now()

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(add_poa,collection_dir=collection_dir, baseDir=baseDir), clusters_id_list)

    representative_fasta = os.path.join(collection_dir, 'reference_pangenome.fasta')
    with open(representative_fasta, 'w') as fh:
        for cluster_id, cons_seq in enumerate(results):
            cons_seq = re.sub('-', '', cons_seq) # remove dash character in consensus sequences
            new_record = SeqRecord(Seq(cons_seq), id = str(cluster_id), description = '')
            SeqIO.write(new_record, fh, 'fasta')

    elapsed = datetime.now() - starttime
    logging.info(f'Add new sequences to previous multiple sequence alignment by abPOA -- time taken {str(elapsed)}')

    return representative_fasta

def main_create_msa(clusters, samples, collection_dir, baseDir, threads):
    clusters_dir = os.path.join(collection_dir, 'clusters')
    if not os.path.exists(clusters_dir):
        os.mkdir(clusters_dir)
    else:
        shutil.rmtree(clusters_dir)
        os.mkdir(clusters_dir)

    gene_to_cluster_id = {gene:i for i, cluster in enumerate(clusters) for gene in cluster}
    clusters_id_list = range(0, len(clusters))

    create_pro_file_for_each_cluster(samples, gene_to_cluster_id, collection_dir)
    representative_fasta = create_poa_in_parallel(clusters_id_list, collection_dir, baseDir, collection_dir, threads)

    return representative_fasta

def add_create_msa(old_clusters, new_clusters, new_samples, collection_dir, baseDir, threads):
    gene_to_cluster_id = {} # find the cluster id for a specific gene
    for i, cluster in enumerate(old_clusters):
        for gene in cluster:
            gene_to_cluster_id[gene] = i
    for i, cluster in enumerate(new_clusters, len(old_clusters)): # continue the id of old clusters
        for gene in cluster:
            gene_to_cluster_id[gene] = i

    create_pro_file_for_each_cluster(new_samples, gene_to_cluster_id, collection_dir)

    # add new gene to existing clusters
    cluster_id_list = range(0, len(old_clusters))
    old_representative_fasta = add_poa_in_parallel(cluster_id_list, collection_dir, baseDir, threads)

    # create msa for the new clusters
    cluster_id_list = range(len(old_clusters) , len(old_clusters) + len(new_clusters))
    out_dir = os.path.join(collection_dir, 'temp')
    new_representative_fasta = create_poa_in_parallel(cluster_id_list, collection_dir, baseDir, out_dir, threads)

    os.system(f'cat {new_representative_fasta} >> {old_representative_fasta}')

    return new_representative_fasta
