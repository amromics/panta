import csv
import gzip
import json
import re
import multiprocessing
from glob import glob
from sklearn.metrics.cluster import rand_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import normalized_mutual_info_score

def get_cluster_panta(input_file):
    gene_to_cluster_id = {}
    reader = csv.reader(gzip.open(input_file, 'rt'), delimiter=',')
    next(reader)
    for cluster_id, row in enumerate(reader):
        for cell in row[1:]: # exclude the ID column
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                gene_to_cluster_id[gene] = cluster_id      
    return gene_to_cluster_id

def get_cluster_roary(input_file):
    gene_to_cluster_id = {}
    reader = csv.reader(open(input_file, 'r'), delimiter=',')
    next(reader)
    for cluster_id, row in enumerate(reader):
        for cell in row[14:]: # exclude the cluster info columns
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                gene_to_cluster_id[gene] = cluster_id      
    return gene_to_cluster_id

def get_cluster_panaroo(input_file):
    gene_to_cluster_id = {}
    reader = csv.reader(open(input_file, 'r'), delimiter=',')
    next(reader)
    for cluster_id, row in enumerate(reader):
        for cell in row[14:]: # exclude the cluster info columns
            genes = cell.split(';')
            for gene in genes:
                if gene == '':
                    continue
                gene_to_cluster_id[gene] = cluster_id      
    return gene_to_cluster_id


def get_cluster_pirate(input_file, prirate_prev_gene_id):
    gene_to_cluster_id = {}
    reader = csv.reader(open(input_file, 'r'), delimiter='\t')
    next(reader)
    for cluster_id, row in enumerate(reader):
        for cell in row[20:]: # exclude the cluster info columns
            cell = re.sub(r'[\(\)]', '', cell)
            cell = re.sub(r';', ':', cell)
            genes = cell.split(':')
            for gene in genes:
                if gene == '':
                    continue
                prev_gene_id = prirate_prev_gene_id[gene]
                gene_to_cluster_id[prev_gene_id] = cluster_id      
    return gene_to_cluster_id



def get_cluster_panx(input_file):
    gene_to_cluster_id = {}
    reader = csv.reader(open(input_file, 'r'), delimiter='\t')
    for cluster_id, row in enumerate(reader):
        for cell in row:
            elements = cell.split('|')
            sample_id = elements[0]
            gene_id = elements[1]
            id_count = re.sub(r'[A-Z_]+', '', gene_id)
            new_gene_id = sample_id + "_" + id_count
            
            gene_to_cluster_id[new_gene_id] = cluster_id      
    return gene_to_cluster_id


def get_prev_gene_id(gff_file):
    prirate_prev_gene_id = {}
    found_fasta = False
    with open(gff_file,'r') as in_fh:
        for line in in_fh:
            if found_fasta == True:
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
                continue
            if re.match(r"^#", line) != None:
                continue
            line = line.rstrip('\n')
            cells = line.split('\t')
            if cells[2] != 'CDS':
                continue
            tags = cells[8].split(';')
            gene_id = None
            prev_gene_id = None
            for tag in tags:
                ID = re.match(r"^ID=(.+)", tag)
                if ID != None:
                    gene_id = ID.group(1)
                    continue
                
                prev_ID = re.match(r"^prev_ID=(.+)", tag)
                if prev_ID != None:
                    prev_gene_id = prev_ID.group(1)
                    continue
            
            if gene_id == None or prev_gene_id == None:
                print('Error')
                continue
            else:
                prirate_prev_gene_id[gene_id] = prev_gene_id
        
    
    return prirate_prev_gene_id

def get_prev_gene_id_pirate(modified_gff_files):
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(get_prev_gene_id, modified_gff_files)

    prirate_prev_gene_id = {}
    for result in results:
        prirate_prev_gene_id.update(result)
    return prirate_prev_gene_id


def compute_similarity_score(gene_to_cluster_id_1, gene_to_cluster_id_2):
    cluster_id_ls_1 = []
    cluster_id_ls_2 = []
    for gene in gene_to_cluster_id_1:
        cluster_id_1 = gene_to_cluster_id_1[gene]
        if gene in gene_to_cluster_id_2:
            cluster_id_2 = gene_to_cluster_id_2[gene]
            cluster_id_ls_1.append(cluster_id_1)
            cluster_id_ls_2.append(cluster_id_2)

    
    # rand_index = rand_score(cluster_id_ls_1, cluster_id_ls_2)
    adjusted_rand_index = adjusted_rand_score(cluster_id_ls_1, cluster_id_ls_2)
    # adjusted_mutual_info = adjusted_mutual_info_score(cluster_id_ls_1, cluster_id_ls_2)
    # normalized_mutual_info = normalized_mutual_info_score(cluster_id_ls_1, cluster_id_ls_2)


    # number_of_shared_gene = len(cluster_id_ls_1)
    # print('Number of gene:', len(gene_to_cluster_id_1), len(gene_to_cluster_id_2), sep=' ')
    # print('Number of shared gene: ', number_of_shared_gene)
    # print('Rand Index: ', round(rand_index, 4) )
    print('Adjusted Rand Index: ', round(adjusted_rand_index, 4) )
    # print('Adjusted mutual info score: ', round(adjusted_mutual_info, 4) )
    # print('Normalized mutual info score: ', round(normalized_mutual_info, 4) )


def compare_tools():
    panta_file = f'{baseDir}/{collection}/out/panta/gene_presence_absence.csv.gz'
    roary_file = f'{baseDir}/{collection}/out/roary_nosplit/gene_presence_absence.csv'
    panaroo_file = f'{baseDir}/{collection}/out/panaroo/gene_presence_absence_roary.csv'
    pirate_file = f'{baseDir}/{collection}/out/PIRATE/PIRATE.gene_families.tsv'
    panx_file = f'{baseDir}/{collection}/gbk/allclusters_final.tsv'

    pirate_modified_gffs = glob(f'{baseDir}/{collection}/out/PIRATE/modified_gffs/*.gff')
    prirate_prev_gene_id = get_prev_gene_id_pirate(pirate_modified_gffs)


    gene_to_cluster_id_1 = get_cluster_panta(panta_file)
    gene_to_cluster_id_2 = get_cluster_roary(roary_file)
    gene_to_cluster_id_3 = get_cluster_panaroo(panaroo_file)
    gene_to_cluster_id_4 = get_cluster_pirate(pirate_file, prirate_prev_gene_id)
    gene_to_cluster_id_5 = get_cluster_panx(panx_file)


    print('Panta vs Roary')
    compute_similarity_score(gene_to_cluster_id_1, gene_to_cluster_id_2)
    print('Panta vs Panaroo')
    compute_similarity_score(gene_to_cluster_id_1, gene_to_cluster_id_3)
    print('Panta vs PIRATE')
    compute_similarity_score(gene_to_cluster_id_1, gene_to_cluster_id_4)
    print('Panta vs PanX')
    compute_similarity_score(gene_to_cluster_id_1, gene_to_cluster_id_5)
    print('Roary vs Panaroo')
    compute_similarity_score(gene_to_cluster_id_2, gene_to_cluster_id_3)
    print('Roary vs PIRATE')
    compute_similarity_score(gene_to_cluster_id_2, gene_to_cluster_id_4)
    print('Roary vs PanX')
    compute_similarity_score(gene_to_cluster_id_2, gene_to_cluster_id_5)
    print('Panaroo vs PIRATE')
    compute_similarity_score(gene_to_cluster_id_3, gene_to_cluster_id_4)
    print('Panaroo vs PanX')
    compute_similarity_score(gene_to_cluster_id_3, gene_to_cluster_id_5)
    print('PIRATE vs PanX')
    compute_similarity_score(gene_to_cluster_id_4, gene_to_cluster_id_5)


def compare_add_pipeline():
    panta_file = f'{baseDir}/{collection}/out/panta/gene_presence_absence.csv.gz'
    gene_to_cluster_id_1 = get_cluster_panta(panta_file)

    for n in [2,25,50,75,99]:
        panta_add_file = f'{baseDir}/{collection}/out/panta_add/{str(n)}/gene_presence_absence.csv.gz'
        gene_to_cluster_id_add = get_cluster_panta(panta_add_file)
        print(n)
        compute_similarity_score(gene_to_cluster_id_1, gene_to_cluster_id_add)

def compare_add_pipeline_other_tools():
    panta_file = f'{baseDir}/{collection}/out/panta/gene_presence_absence.csv.gz'
    roary_file = f'{baseDir}/{collection}/out/roary_nosplit/gene_presence_absence.csv'
    panaroo_file = f'{baseDir}/{collection}/out/panaroo/gene_presence_absence_roary.csv'
    pirate_file = f'{baseDir}/{collection}/out/PIRATE/PIRATE.gene_families.tsv'
    panx_file = f'{baseDir}/{collection}/gbk/allclusters_final.tsv'

    pirate_modified_gffs = glob(f'{baseDir}/{collection}/out/PIRATE/modified_gffs/*.gff')
    prirate_prev_gene_id = get_prev_gene_id_pirate(pirate_modified_gffs)


    gene_to_cluster_id_1 = get_cluster_panta(panta_file)
    gene_to_cluster_id_2 = get_cluster_roary(roary_file)
    gene_to_cluster_id_3 = get_cluster_panaroo(panaroo_file)
    gene_to_cluster_id_4 = get_cluster_pirate(pirate_file, prirate_prev_gene_id)
    gene_to_cluster_id_5 = get_cluster_panx(panx_file)

    for n in [2,25,50,75,99]:
        panta_add_file = f'{baseDir}/{collection}/out/panta_add/{str(n)}/gene_presence_absence.csv.gz'
        gene_to_cluster_id_add = get_cluster_panta(panta_add_file)
        print(n)
        compute_similarity_score(gene_to_cluster_id_1, gene_to_cluster_id_add)
        compute_similarity_score(gene_to_cluster_id_2, gene_to_cluster_id_add)
        compute_similarity_score(gene_to_cluster_id_3, gene_to_cluster_id_add)
        compute_similarity_score(gene_to_cluster_id_4, gene_to_cluster_id_add)
        compute_similarity_score(gene_to_cluster_id_5, gene_to_cluster_id_add)



if __name__ == '__main__':
    global baseDir
    baseDir = '/home/ntanh1999'
    global collection
    collection = 'Pa100'
    global threads
    threads = 8

    # compare_tools()

    # compare_add_pipeline()
    compare_add_pipeline_other_tools()
