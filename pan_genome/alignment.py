import os
import re
import logging
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pyabpoa as pa
from datetime import datetime
import pan_genome.utils as utils
import multiprocessing
from functools import partial

def create_pro_file_for_each_cluster(samples, gene_to_cluster, collection_dir):
    starttime = datetime.now()
    
    for sample in samples:
        sample_id = sample['id']
        faa_file = os.path.join(collection_dir, 'samples', sample_id, sample_id + '.faa')
        with open(faa_file, 'r') as in_fh:
            last_seq_id = None
            seq_lines = []
            for line in in_fh:
                line = line.rstrip()
                result = re.match(r"^>(.+)", line)
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


def run_poa(cluster_id, collection_dir):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(collection_dir, 'clusters', cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + '.faa')
    seqs = []
    id_list = []
    with open(seq_file, 'r') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            seq_id = record.id
            seq = str(record.seq)
            seqs.append(seq)
            id_list.append(seq_id)

    a = pa.msa_aligner(is_aa=True)
    res=a.msa(seqs, out_cons=True, out_msa=True)
    msa = res.msa_seq[:-1] # the last one is consensus sequence
    cons_seq = res.msa_seq[-1]
    msa_file = os.path.join(cluster_dir, cluster_id + '.msa')
    with open(msa_file, 'w') as fh:
        for seq_id, seq in zip(id_list, msa): 
            new_record = SeqRecord(Seq(seq), id = seq_id, description = '')
            SeqIO.write(new_record, fh, 'fasta')

    return cons_seq

def run_poa_in_parrallel(num_clusters, collection_dir, threads):
    starttime = datetime.now()

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(run_poa,collection_dir=collection_dir), range(num_clusters))

    cons_file = os.path.join(collection_dir, 'consensus_sequences.fasta')
    with open(cons_file, 'w') as fh:
        for cluster_id, cons_seq in enumerate(results):
            new_record = SeqRecord(Seq(cons_seq), id = str(cluster_id), description = '')
            SeqIO.write(new_record, fh, 'fasta')

    elapsed = datetime.now() - starttime
    logging.info(f'Create multiple sequence alignment by abPOA -- time taken {str(elapsed)}')


def create_msa(clusters, samples, collection_dir, threads):
    clusters_dir = os.path.join(collection_dir, 'clusters')
    if not os.path.exists(clusters_dir):
        os.mkdir(clusters_dir)

    gene_to_cluster_id = {gene:i for i, cluster in enumerate(clusters) for gene in cluster}
    num_clusters = len(clusters)

    create_pro_file_for_each_cluster(samples, gene_to_cluster_id, collection_dir)
    run_poa_in_parrallel(num_clusters, collection_dir, threads)


