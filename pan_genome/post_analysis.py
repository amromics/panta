import os
import logging
from datetime import datetime
from pan_genome.utils import run_command

logger = logging.getLogger(__name__)


def find_paralogs(cluster, gene_annotation):
    samples = {}
    for gene_id in cluster:
        sample_id = gene_annotation[gene_id]['sample_id']
        samples.setdefault(sample_id, []).append(gene_id)
    
    # pick paralogs with the smallest number of genes
    smallest_number = 1000000
    paralog_genes = None
    for sample_id in samples:
        genes = samples[sample_id]
        count = len(genes)
        if count > 1 and count < smallest_number:
            paralog_genes = genes
            smallest_number = count

    return paralog_genes


def create_orthologs(cluster, paralog_genes, gene_annotation, gene_to_cluster_index):
    # find cluster indices of all the neighbour genes of each paralog gene
    cluster_indices_around_paralogs = []
    for p in paralog_genes:
        neighbours_of_p = gene_annotation[p]['neighbour_genes']
        cluster_indices_around_p = set()
        for neighbour_gene in neighbours_of_p:
            try:
                cluster_index = gene_to_cluster_index[neighbour_gene]
            except:
                continue
            cluster_indices_around_p.add(cluster_index)
        cluster_indices_around_paralogs.append(cluster_indices_around_p)

    # create data structure to hold new clusters
    new_clusters = [[p] for p in paralog_genes]
    new_clusters.append([]) # extra "leftovers" list to gather genes that don't share CGN with any paralog gene

    # add other members of the cluster to their closest match
    for g in cluster:
        if g in paralog_genes:
            continue

        neighbour_genes_of_g = gene_annotation[g]['neighbour_genes']
        if len(neighbour_genes_of_g) == 0:
            new_clusters[-1].append(g)
            continue

        # find paralog gene which is closest match with g
        best_score = 0
        best_score_index = -1  # -1 is the index of "leftovers" list
        for p,_ in enumerate(paralog_genes):
            cluster_indices_around_p = cluster_indices_around_paralogs[p]
            score_of_p = 0
            for neighbour_gene in neighbour_genes_of_g:
                try:
                    cluster_index = gene_to_cluster_index[neighbour_gene]
                except:
                    continue
                if cluster_index in cluster_indices_around_p:
                    score_of_p += 1
            score_of_p = score_of_p / len(neighbour_genes_of_g)
            if score_of_p > best_score:
                best_score = score_of_p
                best_score_index = p

        new_clusters[best_score_index].append(g)

    # check for "leftovers", remove if absent
    if len(new_clusters[-1]) == 0:
        del new_clusters[-1]
    
    return new_clusters

def split_paralogs(gene_annotation, unsplit_clusters, dontsplit):
    if dontsplit == True:
        return unsplit_clusters

    starttime = datetime.now()
    
    clusters_not_paralogs = set()
    # run iteratively
    out_clusters = unsplit_clusters
    for i in range(50):
        stime = datetime.now()
        in_clusters = out_clusters
        out_clusters = []
        any_paralogs = 0
        # convert in_clusters so we can find the cluster index of gene
        gene_to_cluster_index = {gene:index for index, genes in enumerate(in_clusters) for gene in genes}
        
        for cluster in in_clusters:
            if len(cluster) == 1:
                out_clusters.append(cluster)
                continue
            first_gene = cluster[0]
            if first_gene in clusters_not_paralogs:
                out_clusters.append(cluster)
                continue

            # check paralogs
            paralog_genes = find_paralogs(cluster, gene_annotation)

            if paralog_genes == None:
                clusters_not_paralogs.add(first_gene)
                out_clusters.append(cluster)
                continue

            # split paralogs
            orthologs_clusters = create_orthologs(cluster, paralog_genes, gene_annotation, gene_to_cluster_index)
            out_clusters.extend(orthologs_clusters)
            any_paralogs = 1

        # check if next iteration is required
        if any_paralogs == 0:
            break
        elapsed = datetime.now() - stime
        logging.info(f'Split paralogs iterate {i}-- time taken {str(elapsed)}')
    split_clusters = out_clusters

    elapsed = datetime.now() - starttime
    logging.info(f'Split paralogs -- time taken {str(elapsed)}')
    return split_clusters


def annotate_cluster(unlabeled_clusters, gene_annotation):
    starttime = datetime.now()

    clusters = {'groups_' + str(i) : cluster for i, cluster in enumerate(unlabeled_clusters)}

    annotated_clusters = {}
    suffix = 1
    for cluster_name in clusters:
        cluster_new_name = cluster_name
        cluster_product = None
        gene_name_count = {}
        max_number = 0
        gene_id_list = clusters[cluster_name]
        for gene_id in gene_id_list:
            this_gene = gene_annotation[gene_id]
            if 'name' in this_gene:
                gene_name = this_gene['name']
                gene_name_count[gene_name] = gene_name_count.get(gene_name, 0) + 1
                if gene_name_count[gene_name] > max_number:
                    cluster_new_name = gene_name
                    max_number = gene_name_count[gene_name]
                    if 'product' in this_gene:
                        cluster_product = this_gene['product']
        if cluster_product == None:
            cluster_product =[]
            for gene_id in gene_id_list:
                this_gene = gene_annotation[gene_id]
                if 'product' in this_gene:
                    gene_product = this_gene['product']
                    if gene_product not in cluster_product:
                        cluster_product.append(gene_product)
            if len(cluster_product) > 0:
                cluster_product = ', '.join(cluster_product)
            else:
                cluster_product = 'unknown'
        # check if cluster_new_name is already exist
        if cluster_new_name in annotated_clusters:
            cluster_new_name += '_{:05d}'.format(suffix)
            suffix += 1
        annotated_clusters[cluster_new_name] = {'gene_id':gene_id_list, 'product':cluster_product}
    
    elapsed = datetime.now() - starttime
    logging.info(f'Annotate clusters -- time taken {str(elapsed)}')
    return annotated_clusters