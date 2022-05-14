from glob import glob
import os
import subprocess
import multiprocessing
import re
import shutil
from Bio import SeqIO


def create_poa(seq_file):
    matrix_file = os.path.join(baseDir, 'BLOSUM62.mtx')
    result_file = seq_file + '.poa'
    cmd = f'abpoa {seq_file} -o {result_file} -r1 -t {matrix_file} -O 11,0 -E 1,0 -p -c 2> /dev/null'
    os.system(cmd)

def create_poa_in_parallel(seq_files):

    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(create_poa, seq_files)

def create_mafft(seq_file):
    msa_file = seq_file + '.mafft'
    cmd = f"mafft --auto --quiet --thread 1 {seq_file} > {msa_file}" # 9:43
    # cmd = f"mafft --retree 2 --quiet --thread 1 {seq_file} > {msa_file}"  # 3:36
    # cmd = f"mafft --retree 1 --maxiterate 0  --quiet --thread 1 {seq_file} > {msa_file}" # 3:21
    # cmd = f"mafft --auto --maxiterate 10 --quiet --thread 1 {seq_file} > {msa_file}" # 10:01
    
    os.system(cmd)

def create_mafft_in_parallel(seq_files):

    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(create_mafft, seq_files)


def run_baliscore(seq_file):

    base_name = seq_file.rsplit('.', 1)[0]
    ref_file = base_name + '.xml'
    test_file = seq_file + ext


    cmd = f'bali-score -r {ref_file} -t {test_file}'

    output = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result = re.match(r'Test=.+;Ref=.+;Q=([\d\.]+);TC=([\d\.]+)', output.stdout.rstrip())
    result = re.match(r'AlignmentScores { sum_of_pairs: ([\d\.]+), column_score: ([\d\.]+) }', output.stdout.rstrip())
    if result == None:
        return base_name, 0, 0
    else:
        SPS = result.group(1)
        CS = result.group(2)
        return base_name, SPS, CS

def compare(seq_files):
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(run_baliscore, seq_files)

    sum_SPS = 0
    sum_CS = 0
    for cluster_id, SPS, CS in results:
        # print(cluster_id, SPS, CS, sep='\t')
        sum_SPS += float(SPS)
        sum_CS += float(CS)

    print(sum_SPS)
    print(sum_CS)

if __name__ == "__main__":
    global threads
    threads = 8
    
    global baseDir
    baseDir = '/home/noideatt/TA/pan-genome'

    global ext # extension (.mafft or .poa)
    

    seq_files = glob(f'/home/noideatt/TA/evaluate_msa/baliBase/bb3_release/RV50/*.tfa')

    create_poa_in_parallel(seq_files)
    create_mafft_in_parallel(seq_files)
    ext = '.poa'
    compare(seq_files)

    ext = '.mafft'
    compare(seq_files)



    