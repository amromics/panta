from Bio.Align import substitution_matrices
from Bio import AlignIO
import multiprocessing
from glob import glob
import os


def pairwise_score(S1S2):
    S1, S2 = S1S2
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



def multiple_alignment_score(msa_file):
    align = AlignIO.read(msa_file, "fasta")

    msa_score = 0
    compare_list = []
    for i, S1 in enumerate(align):
        for j, S2 in enumerate(align):
            if i > j:
                compare_list.append((S1,S2))

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(pairwise_score, compare_list)
    
    for score in results:
        msa_score += score

    return int(msa_score)


def get_gene_list(gene_id_file):
    gene_id_list = []
    for line in open(gene_id_file):
        gene_id = line.rstrip()
        gene_id_list.append(gene_id)

    return gene_id_list


def compare_collection():
    clusters_id_file = '/home/noideatt/TA/evaluate_msa/clusters_id.txt'
    clusters_id_list = get_gene_list(clusters_id_file)
    clusters_dir = '/home/noideatt/TA/evaluate_msa/out/Pa100'
    for ID in clusters_id_list:
        msa_file = f"{clusters_dir}/{str(ID)}/{str(ID)}.poa"
        poa = multiple_alignment_score(msa_file)

        msa_file = f"{clusters_dir}/{str(ID)}/{str(ID)}.mafft"
        mafft = multiple_alignment_score(msa_file)

        msa_file = f"{clusters_dir}/{str(ID)}/{str(ID)}.new"
        new = multiple_alignment_score(msa_file)

        print(poa, mafft, new, sep = '\t')





if __name__ == "__main__":
    global threads
    threads = 8
    global matrix
    matrix = substitution_matrices.load("BLOSUM62")
    global gop
    gop = -11
    global gep
    gep = -1

    compare_collection()
