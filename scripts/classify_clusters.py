import csv
from functools import total_ordering
import gzip
import re

def get_cluster_sample_size(cluster):
    samples = set()
    for gene in cluster:
        cells = gene.rsplit('_', 1)
        sample = cells[0]
        samples.add(sample)
    
    sample_size = len(samples)
    return sample_size



def get_cluster_panta(input_file):
    clusters = []
    
    reader = csv.reader(gzip.open(input_file, 'rt'), delimiter=',')
    next(reader)
    for row in reader:
        this_cluster = set()
        for cell in row[1:]: # exclude the ID column
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                this_cluster.add(gene)
        sample_count = get_cluster_sample_size(this_cluster)
        clusters.append(sample_count)

                
    return clusters

def get_cluster_roary(input_file):
    clusters = []
    
    reader = csv.reader(open(input_file, 'r'), delimiter=',')
    next(reader)
    for row in reader:
        this_cluster = set()
        for cell in row[14:]: # exclude the cluster info columns
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue  
                this_cluster.add(gene)
        sample_count = get_cluster_sample_size(this_cluster)
        clusters.append(sample_count)
    return clusters

def get_cluster_panaroo(input_file):
    clusters = []
    
    reader = csv.reader(open(input_file, 'r'), delimiter=',')
    next(reader)
    for row in reader:
        this_cluster = set()
        for cell in row[14:]: # exclude the cluster info columns
            genes = cell.split(';')
            for gene in genes:
                if gene == '':
                    continue
                this_cluster.add(gene)
        sample_count = get_cluster_sample_size(this_cluster)
        clusters.append(sample_count)
    return clusters


def get_cluster_pirate(input_file):
    clusters = []
    
    reader = csv.reader(open(input_file, 'r'), delimiter='\t')
    next(reader)
    for row in reader:
        this_cluster = set()
        for cell in row[20:]: # exclude the cluster info columns
            cell = re.sub(r'[\(\)]', '', cell)
            cell = re.sub(r';', ':', cell)
            genes = cell.split(':')
            for gene in genes:
                if gene == '':
                    continue
                this_cluster.add(gene)
        sample_count = get_cluster_sample_size(this_cluster)
        clusters.append(sample_count)
    return clusters

def get_cluster_panx(input_file):
    clusters = []
    reader = csv.reader(open(input_file, 'r'), delimiter='\t')
    for row in reader:
        this_cluster = set()
        for cell in row:
            elements = cell.split('|')
            sample_id = elements[0]
            sample_id = re.sub(r'\.', '_', sample_id)
            gene_id = elements[1]
            id_count = re.sub(r'[A-Z_]+', '', gene_id)
            new_gene_id = sample_id + "_" + id_count
            this_cluster.add(new_gene_id)
        sample_count = get_cluster_sample_size(this_cluster)
        clusters.append(sample_count)
    return clusters


def classify_cluster(sample_size_ls):
    # 4 value corespond to count of core, softcore, shell, cloud gene
    count_list = [0] * 4 
    for num_sample in sample_size_ls:
        percent = num_sample / total
        if percent >= 0.99:
            count_list[0] += 1 # core
        elif percent >= 0.95:
            count_list[1] += 1 # softcore
        elif percent >= 0.15:
            count_list[2] += 1 # shell
        else:
            count_list[3] += 1 # cloud

    print(count_list)
    print(sum(count_list))


if __name__ == '__main__':
    
    baseDir = '/home/noideatt/TA'
    collection = 'Sp100'
    
    global total
    total = 100

    panta_file = f'{baseDir}/{collection}/out/panta/gene_presence_absence.csv.gz'
    roary_file = f'{baseDir}/{collection}/out/roary_nosplit/gene_presence_absence.csv'
    panaroo_file = f'{baseDir}/{collection}/out/panaroo/gene_presence_absence_roary.csv'
    pirate_file = f'{baseDir}/{collection}/out/PIRATE/PIRATE.gene_families.tsv'
    panx_file = f'{baseDir}/{collection}/out/panx/allclusters_final.tsv'

    # clusters = get_cluster_panta(panta_file)
    # classify_cluster(clusters)
    # clusters = get_cluster_roary(roary_file)
    # classify_cluster(clusters)
    # clusters = get_cluster_panaroo(panaroo_file)
    # classify_cluster(clusters)
    # clusters = get_cluster_pirate(pirate_file)
    # classify_cluster(clusters)
    clusters = get_cluster_panx(panx_file)
    # classify_cluster(clusters)

    

    for i in clusters:
        print(i)