import os
import subprocess
import multiprocessing
import re
import shutil
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def create_poa(cluster_id):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(clusters_dir, cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + '.faa')

    result_file = os.path.join(cluster_dir, cluster_id + '.result')
    cmd = f'abpoa {seq_file} -o {result_file} -r2 -t {matrix_file} -O 11,0 -E 1,0 -p -c 2> /dev/null'
    os.system(cmd)

def create_poa_in_parallel(clusters_id_list):

    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(create_poa, clusters_id_list)

def create_mafft(cluster_id):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(clusters_dir, cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + '.faa')
    msa_file = os.path.join(cluster_dir, cluster_id + '.mafft')
    cmd = f"mafft --auto --quiet --thread 1 {seq_file} > {msa_file}" # 9:43
    cmd = f"mafft --retree 2 --quiet --thread 1 {seq_file} > {msa_file}"  # 3:36
    cmd = f"mafft --retree 1 --maxiterate 0  --quiet --thread 1 {seq_file} > {msa_file}" # 3:21
    cmd = f"mafft --auto --maxiterate 10 --quiet --thread 1 {seq_file} > {msa_file}" # 10:01
    
    os.system(cmd)

def create_mafft_in_parallel(clusters_id_list):

    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(create_mafft, clusters_id_list)


def run_qscore(cluster_id):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(clusters_dir, cluster_id)
    poa_file = os.path.join(cluster_dir, cluster_id + '.aln.poa')
    mafft_file = os.path.join(cluster_dir, cluster_id + '.mafft')
    cmd = f"qscore -test {poa_file} -ref {mafft_file}"
    output = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result = re.match(r'Test=.+;Ref=.+;Q=([\d\.]+);TC=([\d\.]+)', output.stdout.rstrip())
    if result == None:
        return cluster_id, 0, 0
    else:
        SPS = result.group(1)
        CS = result.group(2)
        return cluster_id, SPS, CS
    

def get_gene_list(gene_id_file):
    gene_id_list = []
    for line in open(gene_id_file):
        gene_id = line.rstrip()
        gene_id_list.append(gene_id)

    return gene_id_list

def compare(clusters_id_list):
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(run_qscore, clusters_id_list)

    for cluster_id, SPS, CS in results:
        print(cluster_id, SPS, CS, sep='\t')

def rewrite_fasta_1_line(cluster_id):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(clusters_dir, cluster_id)
    msa_file = os.path.join(cluster_dir, cluster_id + '.mafft')
    temp_file = msa_file + '.temp'
    with open(msa_file, 'r') as in_fh, open(temp_file, 'w') as out_fh:
        fasta_out = SeqIO.FastaIO.FastaWriter(out_fh, wrap=None)
        for seq_record in SeqIO.parse(in_fh, 'fasta'):
            fasta_out.write_record(seq_record)

    shutil.move(temp_file, msa_file)

def rewrite_mafft(clusters_id_list):
    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(rewrite_fasta_1_line, clusters_id_list)



def call_snp(cluster_id):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(clusters_dir, cluster_id)
    vcf_file = os.path.join(cluster_dir, cluster_id + '.vcf')
    msa_file = os.path.join(cluster_dir, cluster_id + '.mafft')
    cmd = f'snp-sites -r -v -o {vcf_file} {msa_file}'
    os.system(cmd)

def call_snp_parallel(clusters_id_list):
    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(call_snp, clusters_id_list)

def split_by_length(seq_list, ratio):
    first_seqs = []
    second_seqs = []
    
    index_jump = int(1/ratio)

    for i, seq_tuple in enumerate(seq_list,0):
        if i % index_jump == 0:
            first_seqs.append(seq_tuple)
        else:
            second_seqs.append(seq_tuple)
    
    return first_seqs, second_seqs


def split_by_length_2(seq_list):
    prev_len = 1
    first_seqs = []
    second_seqs = []
    for seq_tuple in seq_list:
        seq_len = len(seq_tuple[1])
        if seq_len / prev_len > 1.1:
            first_seqs.append(seq_tuple)
            prev_len = seq_len
        else:
            second_seqs.append(seq_tuple)

    return first_seqs, second_seqs


def split_random(seq_list, ratio):
    random.seed(62)
    random.shuffle(seq_list)
    total_num = len(seq_list)
    first_seqs_len = int(total_num * ratio)
    first_seqs = seq_list[:first_seqs_len]
    second_seqs = seq_list[first_seqs_len:]
    return first_seqs, second_seqs

def split_seq_file(seq_file):
    seq_list = []
    with open(seq_file, 'r') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            seq_id = record.id
            seq = str(record.seq)
            seq_list.append((seq_id, seq))

    seq_list.sort(key= lambda x:len(x[1])) # sort sequences by length
    
    # first_seqs, second_seqs = split_random(seq_list, 0.2)
    first_seqs, second_seqs = split_by_length(seq_list, 0.2)

    mafft_seq = os.path.join(seq_file + '.1')
    with open(mafft_seq, 'w') as fh:
        for seq_id, seq in first_seqs: 
            new_record = SeqRecord(Seq(seq), id = seq_id, description = '')
            SeqIO.write(new_record, fh, 'fasta')
    
    poa_seq = os.path.join(seq_file + '.2')
    with open(poa_seq, 'w') as fh:
        for seq_id, seq in second_seqs: 
            new_record = SeqRecord(Seq(seq), id = seq_id, description = '')
            SeqIO.write(new_record, fh, 'fasta')

    return mafft_seq, poa_seq


def new_method(cluster_id):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(clusters_dir, cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + '.faa')

    mafft_seq, poa_seq = split_seq_file(seq_file)

    mafft_msa_file = os.path.join(cluster_dir, cluster_id + '.1.aln')
    cmd = f"mafft --localpair --maxiterate 1000 --quiet --thread 1 {mafft_seq} > {mafft_msa_file}"
    os.system(cmd)

    result_file = os.path.join(cluster_dir, cluster_id + '.new')
    cmd = f'abpoa {poa_seq} -i {mafft_msa_file} -o {result_file} -r1 -t {matrix_file} -O 11,0 -E 1,0 -p -c 2> /dev/null'
    os.system(cmd)

    os.remove(poa_seq)
    os.remove(mafft_seq)
    os.remove(mafft_msa_file)


def new_method_parallel(clusters_id_list):

    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(new_method, clusters_id_list)



if __name__ == "__main__":
	
    global threads
    threads = 8

    global matrix_file # Blosum matrix file for POA
    matrix_file = '/home/noideatt/TA/pan-genome/BLOSUM62.mtx'

    global clusters_dir
    clusters_dir = '/home/noideatt/TA/evaluate_msa/out/Pa100'

    clusters_id_file = '/home/noideatt/TA/evaluate_msa/clusters_id.txt'
    clusters_id_list = get_gene_list(clusters_id_file)
    # clusters_id_list =  range(0, 4671)



    # Compare POA and Mafft
    
    # create_poa_in_parallel(clusters_id_list)
    # create_mafft_in_parallel(clusters_id_list)
    # rewrite_mafft(clusters_id_list)
    # compare(clusters_id_list)
    # call_snp_parallel(clusters_id_list)

    # new_method: use mafft first to create high quality graph


    new_method_parallel(clusters_id_list)
    # compare(clusters_id_list)