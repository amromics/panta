from email.mime import base
from glob import glob
import os
import subprocess
import multiprocessing
import re
import shutil
import random
from Bio.Align import substitution_matrices
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def run_command(cmd):
    cmd = f'/usr/bin/time --append -v -o {timing_log} bash -c "{cmd}"'
    os.system(cmd)

def rewrite_fasta(msa_file):
    temp_file = msa_file + '.temp'
    with open(msa_file, 'r') as in_fh, open(temp_file, 'w') as out_fh:
        fasta_out = SeqIO.FastaIO.FastaWriter(out_fh, wrap=None)
        for seq_record in SeqIO.parse(in_fh, 'fasta'):
            fasta_out.write_record(seq_record)
    shutil.move(temp_file, msa_file)


def create_poa(seq_file, msa_file):
    cmd = (f'abpoa {seq_file} -o {msa_file} -r1 -t {matrix_file} '
            '-O 11,0 -E 1,0 -p -c 2> /dev/null')
    run_command(cmd)

def create_poa_in_parallel(input_files, output_files):
    
    with multiprocessing.Pool(processes=threads) as pool:
        pool.starmap(create_poa, zip(input_files, output_files))

def create_mafft(seq_file, msa_file):
    # cmd = f"mafft --auto --quiet --thread 1 {seq_file} > {msa_file}" # 9:43
    # cmd = f"mafft --retree 2 --quiet --thread 1 {seq_file} > {msa_file}"  # 3:36
    # cmd = f"mafft --retree 1 --maxiterate 0  --quiet --thread 1 {seq_file} > {msa_file}" # 3:21
    # cmd = f"mafft --auto --maxiterate 1000 --quiet --thread 1 {seq_file} > {msa_file}" # 10:01
    cmd = f"mafft --localpair --maxiterate 1000 --quiet --thread 1 {seq_file} > {msa_file}"
    
    run_command(cmd)

    rewrite_fasta(msa_file)


def create_mafft_in_parallel(input_files, output_files):
    
    with multiprocessing.Pool(processes=threads) as pool:
        pool.starmap(create_mafft, zip(input_files, output_files))


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


# def split_by_length_2(seq_list):
#     prev_len = 1
#     first_seqs = []
#     second_seqs = []
#     for seq_tuple in seq_list:
#         seq_len = len(seq_tuple[1])
#         if seq_len / prev_len > 1.1:
#             first_seqs.append(seq_tuple)
#             prev_len = seq_len
#         else:
#             second_seqs.append(seq_tuple)

#     return first_seqs, second_seqs


def split_random(seq_list, ratio):
    random.seed(62)
    random.shuffle(seq_list)
    total_num = len(seq_list)
    first_seqs_len = int(total_num * ratio)
    first_seqs = seq_list[:first_seqs_len]
    second_seqs = seq_list[first_seqs_len:]
    return first_seqs, second_seqs

def split_num(seq_list, num):
    random.seed(62)
    random.shuffle(seq_list)
    first_seqs = seq_list[:num]
    second_seqs = seq_list[num:]
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
    # first_seqs, second_seqs = split_by_length(seq_list, 0.2)
    first_seqs, second_seqs = split_num(seq_list, 500)

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


def new_method(seq_file, msa_file):
    mafft_seq, poa_seq = split_seq_file(seq_file)

    tmp_msa_file = msa_file + '.tmp'
    cmd = ("mafft --localpair --maxiterate 1000 --quiet --thread 1 "
           f" {mafft_seq} > {tmp_msa_file}")
    run_command(cmd)

    cmd = (f'abpoa {poa_seq} -i {tmp_msa_file} -o {msa_file} -r1 '
           f'-t {matrix_file} -O 11,0 -E 1,0 -p -c 2> /dev/null')
    run_command(cmd)

    rewrite_fasta(msa_file)

    os.remove(poa_seq)
    os.remove(mafft_seq)
    os.remove(tmp_msa_file)


def new_method_parallel(input_files, output_files):
    
    with multiprocessing.Pool(processes=threads) as pool:
        pool.starmap(new_method, zip(input_files, output_files))


def run_qscore(test_file, ref_file):
    basename_1 = os.path.basename(test_file)
    basename_2 = os.path.basename(ref_file)
    cmd = f"qscore -test {test_file} -ref {ref_file}"
    output = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result = re.match(
        r'Test=.+;Ref=.+;Q=([\d\.]+);TC=([\d\.]+)', 
        output.stdout.rstrip())
    if result == None:
        return 0, 0
    else:
        SPS = result.group(1)
        CS = result.group(2)
        # print(basename_1, basename_2, SPS, CS, sep='\t')
        return SPS, CS


def qscore_parallel(test_list, reference_list):
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.starmap(run_qscore, zip(test_list, reference_list))
    
    sum_SP = 0
    sum_TC = 0
    for SP, TC in results:
        sum_SP += float(SP)
        sum_TC += float(TC)
    
    ave_SP = round(sum_SP / len(results), 4)
    ave_TC = round(sum_TC / len(results), 4)

    print('Average SP', ave_SP, sep=',')
    print('Average TC', ave_TC, sep=',')


def pairwise_score(S1, S2):
    score = 0
    prev_gap = False
    for A, B in zip(S1, S2):
        if (A == '-') and (B == '-'):
            continue
        gap = (A == '-') or (B == '-')
        if gap == False:
            score += matrix[A,B]
        else:
            if prev_gap == False:
                score += gop
            else:
                score += gep
        prev_gap = gap
    return score


def msa_score(msa_file):
    align = AlignIO.read(msa_file, "fasta")
    base_name = os.path.basename(msa_file)

    msa_score = 0
    arguments = []
    for i, S1 in enumerate(align):
        for j, S2 in enumerate(align):
            if i > j:
                arguments.append((S1,S2))

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.starmap(pairwise_score, arguments)
    
    for score in results:
        msa_score += score

    print(base_name, int(msa_score), sep='\t')

def run_iqtree_alisim():
    os.chdir('/home/ntanh1999/simulation')
    
    root_name = 'test_12'
    model = '-m LG '
    indel = '--indel 0.05,0.05 '
    random_tree = '-t RANDOM{yh/2000} -rlen 0.001 0.01 0.999 '
    seq_length = '--length 400 '
    number_msa = '--num-alignments 8'

    cmd = (f"iqtree2 --alisim {root_name} -af fasta " 
            + model + indel + random_tree 
            + seq_length + number_msa)

    os.system(cmd)
    return root_name

def pairwise_identity(S1, S2):
    total = 0
    same = 0
    for A, B in zip(S1, S2):
        if (A == '-') and (B == '-'):
            continue
        if A == B:
            same += 1
        total +=1
    identity = same / total
    # print(round(identity,2))
    return identity

def calucate_average_identity(msa_file):
    align = AlignIO.read(msa_file, "fasta")
    base_name = os.path.basename(msa_file)
    
    pair = []
    for i, S1 in enumerate(align):
        for j, S2 in enumerate(align):
            if i > j:
                pair.append((S1,S2))

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.starmap(pairwise_identity, pair)
    
    identity_score = sum(results) / len(results)

    # print(base_name, round(identity_score,2), sep='\t')
    return identity_score

def calucate_average_identity_2(msa_file):
    align = AlignIO.read(msa_file, "fasta")
    base_name = os.path.basename(msa_file)
    
    pair = []
    S1 = align[0]
    for S2 in align:
        pair.append((S1,S2))

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.starmap(pairwise_identity, pair)
    
    identity_score = sum(results) / len(results)

    # print(base_name, round(identity_score,2), sep='\t')
    return identity_score

def calculate_average_identity(ref_list):
    identity_ls = []
    for ref_msa in ref_list:
        identity = calucate_average_identity_2(ref_msa)
        identity_ls.append(identity)
    ave_identity = round(sum(identity_ls)/len(identity_ls), 2)
    print('Average Identity', ave_identity, sep=',')


def main_1():
    global matrix_file # Blosum matrix file for POA
    matrix_file = '/home/ntanh1999/amromics/amromics/pan-genome/BLOSUM62.mtx'
    global matrix
    matrix = substitution_matrices.load("BLOSUM62")
    global gop
    gop = -11
    global gep
    gep = -1

    root_name = run_iqtree_alisim()
    
    global timing_log
    timing_log = f'/home/ntanh1999/simulation/{root_name}.log'
     
    unaligned_seqs_list = glob(
        f'/home/ntanh1999/simulation/{root_name}_*.unaligned.fa')

    
    ref_list = []
    poa_msa_list = []
    mafft_msa_list = []
    new_method_list = []
    for seq_file in unaligned_seqs_list:
        dir_name = os.path.dirname(seq_file)
        base_name = os.path.basename(seq_file)
        root = base_name.split('.', 1)[0]
        
        ref = os.path.join(dir_name, root + '.fa')
        poa_msa = os.path.join(dir_name, root + '.poa')
        mafft_msa = os.path.join(dir_name, root + '.mafft')
        new_method_msa = os.path.join(dir_name, root + '.new')

        ref_list.append(ref)
        poa_msa_list.append(poa_msa)
        mafft_msa_list.append(mafft_msa)
        new_method_list.append(new_method_msa)
    
    calculate_average_identity(ref_list)

    create_poa_in_parallel(unaligned_seqs_list, poa_msa_list)
    # create_mafft_in_parallel(unaligned_seqs_list, mafft_msa_list)
    new_method_parallel(unaligned_seqs_list, new_method_list)

    qscore_parallel(poa_msa_list, ref_list)
    # qscore_parallel(mafft_msa_list, ref_list)
    qscore_parallel(new_method_list, ref_list)


    # for msa in mafft_msa_list:
    #     msa_score(msa)
    # for msa in poa_msa_list:
    #     msa_score(msa)
    # for msa in new_method_list:
    #     msa_score(msa)
    # for msa in ref_list:
    #     msa_score(msa)


def main_2():
    for i in range(0, 5001):
        msa_path = f'/home/ntanh1999/evaluate_msa/Pa100/{str(i)}/{str(i)}.mafft'
        calucate_average_identity_2(msa_path)

if __name__ == "__main__":
	
    global threads
    threads = 8

    main_1()
    # main_2()