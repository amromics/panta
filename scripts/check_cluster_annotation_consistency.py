import json
import csv
import gzip
import re

def extract_panta(input_file, product_dictionary):
    product_dictionary = {}
    clusters_list = []
    with gzip.open(input_file, 'rt') as fh:
        reader = csv.reader(fh, delimiter=',')
        header = next(reader)
        for cluster_id, row in enumerate(reader):
            cluster_product = set()
            clusters_list.append(cluster_product)
            for cell in row[1:]:
                genes = cell.split('\t')
                for gene in genes:
                    if gene == '':
                        continue
                    product = product_dictionary[gene]
                    if (product != 'hypothetical protein' and 
                        product != "putative protein"):
                        cluster_product.add(product)
                        product_dictionary.setdefault(
                            product, set()).add(cluster_id)
            
    return product_dictionary, clusters_list
            


def extract_panaroo(input_file, product_dictionary):
    product_dictionary = {}
    clusters_list = []
    with open(input_file, 'r') as fh:
        reader = csv.reader(fh, delimiter=',')
        header = next(reader)
        excluded_gene = []
        for cluster_id, row in enumerate(reader):
            cluster_product = set()
            clusters_list.append(cluster_product)
            for cell in row[14:]:
                genes = cell.split(';')
                for gene in genes:
                    if gene == '':
                        continue
                    if gene not in product_dictionary:
                        excluded_gene.append(gene)
                        continue
                    product = product_dictionary[gene]
                    if (product != 'hypothetical protein' and 
                        product != "putative protein"):
                        cluster_product.add(product)
                        product_dictionary.setdefault(
                            product, set()).add(cluster_id)
    print('Excluded gene: ', len(excluded_gene), sep = '')
    return product_dictionary, clusters_list

def extract_roary(input_file, product_dictionary):
    product_dictionary = {}
    clusters_list = []
    with open(input_file, 'r') as fh:
        reader = csv.reader(fh, delimiter=',')
        header = next(reader)
        excluded_gene = []
        for cluster_id, row in enumerate(reader):
            cluster_product = set()
            clusters_list.append(cluster_product)
            for cell in row[14:]:
                genes = cell.split('\t')
                for gene in genes:
                    if gene == '':
                        continue
                    if gene not in product_dictionary:
                        excluded_gene.append(gene)
                        continue
                    product = product_dictionary[gene]
                    if (product != 'hypothetical protein' and 
                        product != "putative protein"):
                        cluster_product.add(product)
                        product_dictionary.setdefault(
                            product, set()).add(cluster_id)
    print('Excluded gene: ', len(excluded_gene), sep = '')
    return product_dictionary, clusters_list

def extract_roary_2(input_file, product_dictionary):
    product_dictionary = {}
    clusters_list = []
    with open(input_file, 'r') as fh:
        reader = csv.reader(fh, delimiter='\t')
        excluded_gene = []
        for cluster_id, row in enumerate(reader):
            cluster_product = set()
            clusters_list.append(cluster_product)
            for cell in row:
                genes = cell.split(' ')
                for gene in genes:
                    if gene == '':
                        continue
                    if gene not in product_dictionary:
                        print(gene)
                        excluded_gene.append(gene)
                        continue
                    product = product_dictionary[gene]
                    if (product != 'hypothetical protein' and 
                        product != "putative protein"):
                        cluster_product.add(product)
                        product_dictionary.setdefault(
                            product, set()).add(cluster_id)
    print('Excluded gene: ', len(excluded_gene), sep = '')
    return product_dictionary, clusters_list

def extract_pirate(input_file, product_dictionary):
    product_dictionary = {}
    clusters_list = []
    with open(input_file, 'r') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        excluded_gene = []
        for cluster_id, row in enumerate(reader):
            cluster_product = set()
            clusters_list.append(cluster_product)
            for cell in row[20:]:
                cell = re.sub(r'[\(\)]', '', cell)
                cell = re.sub(r';', ':', cell)
                genes = cell.split(':')
                for gene in genes:
                    if gene == '':
                        continue
                    if gene not in product_dictionary:
                        excluded_gene.append(gene)
                        continue
                    product = product_dictionary[gene]
                    if (product != 'hypothetical protein' and 
                        product != "putative protein"):
                        cluster_product.add(product)
                        product_dictionary.setdefault(
                            product, set()).add(cluster_id)
    print('Excluded gene: ', len(excluded_gene), sep = '')
    return product_dictionary, clusters_list

def get_product_dictionary(gene_dictionary_file):

    gene_dictionary = json.load(open(gene_dictionary_file, 'r'))

    product_dict = {}
    for gene in gene_dictionary:
        new_gene = re.sub(r'\.', '_', gene)
        product_dict[new_gene] = gene_dictionary[gene][4]

    json.dump(product_dict, open('Sp616.json', 'w'), indent = 4) 



def count_merged_cluster(clusters_list):
    merged = 0
    consistent = 0
    unknown = 0
    for cluster_product in clusters_list:
        if len(cluster_product) == 0:
            unknown += 1
        elif len(cluster_product) > 1:
            merged += 1
        else:
            consistent += 1

    print(f"Merged cluster: {merged}")
    print(f"Consistent cluster: {consistent}")
    print(f"Unknown cluster: {unknown}")

def count_split_cluster(product_dictionary):
    split_product = 0
    for product in product_dictionary:
        # print(product, product_dictionary[product], sep = ': ')
        if len(product_dictionary[product]) > 1:
            split_product += 1 

    print(f"Split product: {split_product}")
    print(f"Total product: {len(product_dictionary)}")



if __name__ == '__main__':
    # raw_gene_dictionary_file = '/home/noideatt/TA/check_consistency/Sp616_raw.json'
    # get_product_dictionary(raw_gene_dictionary_file)
    product_dictionary_file = '/home/noideatt/TA/check_consistency/Sp616.json'
    product_dictionary = json.load(open(product_dictionary_file, 'r'))
    
    # panta
    input_file = '/home/noideatt/TA/out/Sp100/100/gene_presence_absence.csv.gz'
    product_dictionary, clusters_list = extract_panta(
        input_file, product_dictionary)
    
    # pararoo
    # input_file = '/home/noideatt/TA/check_consistency/Pa800_panaroo.csv'
    # product_dictionary, clusters_list = extract_panaroo(
    # input_file, product_dictionary)
    
    # roary
    # input_file = '/home/noideatt/TA/check_consistency/Pa800_roary_nosplit.csv'
    # product_dictionary, clusters_list = extract_roary_2(
    # input_file, product_dictionary)


    # pirate
    # input_file = '/home/noideatt/TA/check_consistency/Pa800_pirate.tsv'
    # product_dictionary, clusters_list = extract_pirate(
    # input_file, product_dictionary)


    count_merged_cluster(clusters_list)
    count_split_cluster(product_dictionary)