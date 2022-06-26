from Bio import SeqIO


def read_database(db_file):
    clusters_annotation = []
    clusters = []
    gene_to_cluster = {}
    index = 0
    for seq_record in SeqIO.parse(open(db_file), "fasta"):
        desc = seq_record.description.split('~~~')
        desc[0] = desc[0].split(' ')[1] # remove sequence ID
        clusters_annotation.append(desc)
        clusters.append([])
        gene_to_cluster[seq_record.id] = index
        index += 1
    
    return clusters, clusters_annotation, gene_to_cluster


def combine_result(old_clusters, old_clusters_annotation, 
                   new_clusters, new_clusters_annotation):
    out_clusters = []
    out_clusters_annotation = []

    for old_cluster, annotation in zip(old_clusters, old_clusters_annotation):
        if len(old_cluster) == 0:
            continue
        else:
            out_clusters.append(old_cluster)
            out_clusters_annotation.append(annotation)
    
    for new_cluster, annotation in zip(new_clusters, new_clusters_annotation):
        out_clusters.append(new_cluster)
        out_clusters_annotation.append(annotation)

    return out_clusters, out_clusters_annotation

def combine_clusters(old_clusters, new_clusters):
    out_clusters = []
    for old_cluster in old_clusters:
        if len(old_cluster) == 0:
            continue
        else:
            out_clusters.append(old_cluster)
    
    out_clusters.extend(new_clusters)

    return out_clusters

