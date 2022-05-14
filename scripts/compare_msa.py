import os
import subprocess
import multiprocessing
import re
import shutil
from Bio import SeqIO


def run_panta(panta_dir, input_dir, out_dir):
	os.chdir(panta_dir)
	if not os.path.isdir(out_dir):
		os.makedirs(out_dir)
	
	cmd = f'/usr/bin/time -v python pan-genome.py main -o {out_dir} -t {threads} -g {input_dir}/*.gff 1>> {out_dir}/panta.log 2>&1'
	os.system(cmd)

	return out_dir

def create_poa(cluster_id):
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(clusters_dir, cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + '.faa')

    matrix_file = os.path.join(baseDir, 'BLOSUM62.mtx')
    result_file = os.path.join(cluster_dir, cluster_id + '.result')
    cmd = f'abpoa {seq_file} -o {result_file} -r1 -p -c -m 1 2> /dev/null'
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
    poa_file = os.path.join(cluster_dir, cluster_id + '.result')
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



if __name__ == "__main__":
	
    global threads
    threads = 8

    panta_dir = '/home/noideatt/TA/pan-genome'
    input_dir = '/home/noideatt/TA/data/Kp26/main'
    out_dir = '/home/noideatt/TA/evaluate_msa/out/Kp20'
    # run_panta(panta_dir, input_dir, out_dir)

    # clusters_id_file = '/home/noideatt/TA/evaluate_msa/clusters_id.txt'
    # clusters_id_list = get_gene_list(clusters_id_file)
    clusters_id_list =  range(0, 9252)
    global clusters_dir
    clusters_dir = '/home/noideatt/TA/evaluate_msa/out/Kp20'
    global baseDir
    baseDir = panta_dir
    # create_poa_in_parallel(clusters_id_list)
    create_mafft_in_parallel(clusters_id_list)
    # rewrite_mafft(clusters_id_list)


    # compare(clusters_id_list)

    # call_snp_parallel(clusters_id_list)