import csv
import gzip
import re

def get_cluster_panta(input_file):
    clusters = []
    gene_set = set()
    reader = csv.reader(gzip.open(input_file, 'rt'), delimiter=',')
    next(reader)
    for row in reader:
        this_cluster = set()
        clusters.append(this_cluster)
        for cell in row[1:]: # exclude the ID column
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                gene_set.add(gene)
                this_cluster.add(gene)
                
    return clusters, gene_set

def get_cluster_roary(input_file):
    clusters = []
    gene_set = set()
    reader = csv.reader(open(input_file, 'r'), delimiter=',')
    next(reader)
    for row in reader:
        this_cluster = set()
        clusters.append(this_cluster)
        for cell in row[14:]: # exclude the cluster info columns
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                gene_set.add(gene)
                this_cluster.add(gene)
                
    return clusters, gene_set

def get_cluster_panaroo(input_file):
    clusters = []
    gene_set = set()
    reader = csv.reader(open(input_file, 'r'), delimiter=',')
    next(reader)
    for row in reader:
        this_cluster = set()
        clusters.append(this_cluster)
        for cell in row[14:]: # exclude the cluster info columns
            genes = cell.split(';')
            for gene in genes:
                if gene == '':
                    continue
                gene_set.add(gene)
                this_cluster.add(gene)
                
    return clusters, gene_set


def get_cluster_pirate(input_file):
    clusters = []
    gene_set = set()
    reader = csv.reader(open(input_file, 'r'), delimiter='\t')
    next(reader)
    for row in reader:
        this_cluster = set()
        clusters.append(this_cluster)
        for cell in row[20:]: # exclude the cluster info columns
            cell = re.sub(r'[\(\)]', '', cell)
            cell = re.sub(r';', ':', cell)
            genes = cell.split(':')
            for gene in genes:
                if gene == '':
                    continue
                gene_set.add(gene)
                this_cluster.add(gene)
                
    return clusters, gene_set


def remove_differnt_gene(clusters_1, gene_set_1, clusters_2, gene_set_2):
    new_cluster_1 = []
    for cluster in clusters_1:
        new_cluster = set()
        for gene in cluster:
            if gene in gene_set_2:
                new_cluster.add(gene)
        if len(new_cluster) != 0:
            new_cluster_1.append(new_cluster)
    new_cluster_2 = []
    for cluster in clusters_2:
        new_cluster = set()
        for gene in cluster:
            if gene in gene_set_1:
                new_cluster.add(gene)
        if len(new_cluster) != 0:
            new_cluster_2.append(new_cluster)

    return new_cluster_1, new_cluster_2

def count_different(clusters_1, clusters_2):
    same = 0
    diff = 0
    total = 0
    for i in clusters_2:
        found = False
        for j in clusters_1:
            if i == j:
                same += 1
                found = True
                break
        if found == False:
            diff += 1
        total += 1

    print(str(same), str(diff), str(total), sep="\n")



if __name__ == '__main__':
    
    panta_file = '/home/noideatt/TA/Sp100/out/panta/gene_presence_absence.csv.gz'
    roary_file = '/home/noideatt/TA/check_consistency/Sp616_roary_nosplit.csv'
    panaroo_file = '/home/noideatt/TA/check_consistency/Sp616_panaroo.csv'
    pirate_file = '/home/noideatt/TA/Sp100/out/PIRATE/PIRATE.gene_families.tsv'

    clusters_1, gene_set_1 = get_cluster_panta(panta_file)
    # clusters_2, gene_set_2 = get_cluster_roary(roary_file)
    # clusters_2, gene_set_2 = get_cluster_panaroo(panaroo_file)
    clusters_2, gene_set_2 = get_cluster_pirate(pirate_file)

    # panta_98 = '/home/noideatt/TA/out/Sp100/50/gene_presence_absence.csv.gz'
    # panta_95 = '/home/noideatt/TA/out/Sp100_new/50/gene_presence_absence.csv.gz'

    # clusters_1, gene_set_1 = get_cluster_panta(panta_98)
    # clusters_2, gene_set_2 = get_cluster_panta(panta_95)

    clusters_1, clusters_2 = remove_differnt_gene(clusters_1, gene_set_1, clusters_2, gene_set_2)
    count_different(clusters_1, clusters_2)
    count_different(clusters_2, clusters_1)