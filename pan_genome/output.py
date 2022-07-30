import os
import csv
import sys
import gzip
import logging
import shutil
from datetime import datetime

logger = logging.getLogger(__name__)

def classify_cluster(num_sample, total, count):
    """
    Clasify a cluster into core, softcore, shell and cloud gene.

    Parameters
    ----------
    num_sample : int
        number of sample in this cluster
    total : int
        total number of sample in collection
    count : list
        a data structure to keep track the number of 
        core, softcore, shell and cloud gene.
    """
    percent = num_sample / total
    if percent >= 0.99:
        count[0] += 1 # core
    elif percent >= 0.95:
        count[1] += 1 # softcore
    elif percent >= 0.15:
        count[2] += 1 # shell
    else:
        count[3] += 1 # cloud

def write_summary(count, out_dir, total_samples):
    """
    Write the summary file.

    Parameters
    ----------
    count : list
        a data structure to keep track the number of 
        core, softcore, shell and cloud gene.
    out_dir : path
        directory of output file. 
    total_samples : int
        Number of samples in collection.   
    """
    total = sum(count)
    summary_file = os.path.join(out_dir, 'summary_statistics.txt')
    with open(summary_file, 'w') as fh:
        fh.write('Core genes' + '\t' + '(99% <= strains <= 100%)' 
                 + '\t'+ str(count[0]) + '\n')
        fh.write('Soft core genes' + '\t' + '(95% <= strains < 99%)' 
                 + '\t'+ str(count[1]) + '\n')
        fh.write('Shell genes' + '\t' + '(15% <= strains < 95%)' 
                 + '\t' + str(count[2]) + '\n')
        fh.write('Cloud genes' + '\t' + '(0% <= strains < 15%)' 
                 + '\t'+ str(count[3]) + '\n')
        fh.write('Total genes' + '\t' + '(0% <= strains <= 100%)' 
                 + '\t'+ str(total) + '\n')
        fh.write('Total number of samples:' + '\t' + str(total_samples))

def get_number_of_samples(summary_file):
    """
    Get total number of samples in previous collection.

    Parameters:
    -----------
    summary_file: path
        path to summary_statistics.txt
    
    Returns:
    --------
    int
        number of samples of previous collection
    """
    with open(summary_file, 'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[0] == 'Total number of samples:':
                total_samples = line[1]
                return int(total_samples)


def output_cluster_info(clusters, clusters_annotation, 
                  gene_dictionary, samples, out_dir):
    """
    Create 2 output file:
        - cluster_info.csv
        - summary_statistics.txt
    
    Parameters
    ----------
    clusters : list of list
        list of sequence IDs of each cluster
    clusters_annotation : list of list
        list of annotation information of each cluster
    gene dictionary : dict
        contain information of each gene
        {gene_id: (sample_id, contig, length, gene_name, gene_product)}
    samples : list of dict
        list of samples information {id: , gff_file: , assembly: }
    out_dir : path
        output directory
    """
    starttime = datetime.now()
    # output files
    cluster_info_file = os.path.join(out_dir, 'cluster_info.csv')
    total_sample = len(samples)
    
    # count: a data structure to keep track the number of 
    # core, softcore, shell and cloud gene.
    count = [0] * 4 

    with open(cluster_info_file, 'w') as fh_1:
        writer_1 = csv.writer(
            fh_1, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)

        # write header
        header_1 = ['ID','Gene', 'Annotation', 'No. isolates', 'No. sequences', 
                    'Avg sequences per isolate', 'Min group size nuc', 
                    'Max group size nuc', 'Avg group size nuc' ]
        writer_1.writerow(header_1)

        # write row
        cluster_id = 0
        for cluster, cluster_annotation in zip(clusters, clusters_annotation):
            row_1 = []
            sample_dict = {}
            length_list = []
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
                length = gene_dictionary[gene][2]
                length_list.append(length)
            
            # ID
            row_1.append(cluster_id)
            # Gene
            row_1.append(cluster_annotation[0])
            # Annotation
            row_1.append(cluster_annotation[1])
            # No. isolates
            num_sample = len(sample_dict)
            row_1.append(num_sample)
            classify_cluster(num_sample, total_sample, count)
            # No. sequences
            row_1.append(len(cluster))
            # Avg sequences per isolate
            avg_seq = len(cluster) / len(sample_dict)
            row_1.append(round(avg_seq,2))
            # Min group size nuc
            row_1.append(min(length_list))
            # Max group size nuc
            row_1.append(max(length_list))
            # Avg group size nuc
            nuc_size = sum(length_list) / len(length_list)
            row_1.append(round(nuc_size,0))

            writer_1.writerow(row_1)

            cluster_id += 1
    # write summary file
    write_summary(count, out_dir, total_sample)

    elapsed = datetime.now() - starttime
    logging.info(f'Output cluster info -- time taken {str(elapsed)}')


def update_cluster_info(
        previous_clusters, new_clusters, new_clusters_annotation, 
        gene_dictionary, new_samples, temp_dir, collection_dir):
    """
    Update cluster_info.csv
    Write a new summary_statistics.txt
    
    Parameters
    ----------
    previous_clusters : list of list
        list of sequence IDs of previous clusters
    new_clusters : list of list
        list of sequence IDs of new clusters
    new_clusters_annotation : list of list
        list of annotation information of new clusters
    gene dictionary : dict
        contain information of each gene
        {gene_id: (sample_id, contig, length, gene_name, gene_product)}
    new_samples : list of dict
        list of new samples {id: , gff_file: , assembly: }
    temp_dir : path
        temporary directory
    collection_dir : path
        collection directory
    """    
    starttime = datetime.now()
    
    new_cluster_info_file = os.path.join(
        temp_dir, 'cluster_info.csv')
    old_cluster_info_file = os.path.join(
        collection_dir, 'cluster_info.csv')
    summary_file = os.path.join(collection_dir, 'summary_statistics.txt')
    total_sample = get_number_of_samples(summary_file) + len(new_samples)

    # count: a data structure to keep track the number of 
    # core, softcore, shell and cloud gene.
    count = [0] * 4 

    with open(old_cluster_info_file, 'r') as in_fh_1, \
        open(new_cluster_info_file, 'w') as out_fh_1:
        csv.field_size_limit(sys.maxsize)
        writer_1 = csv.writer(
            out_fh_1, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        reader_1 = csv.reader(in_fh_1, delimiter=',')
        
        # update header
        header_1 = next(reader_1)
        writer_1.writerow(header_1)
        
        # update row of previous clusters
        cluster_id = 0
        for row_1, cluster in zip(reader_1, previous_clusters):
            num_seq = int(row_1[4])
            num_iso = int(row_1[3])
            min_len = int(row_1[6])
            max_len = int(row_1[7])
            new_sample_dict = {}
            new_length_list = []
            for new_gene in cluster:
                sample_id = gene_dictionary[new_gene][0]
                new_sample_dict.setdefault(sample_id, []).append(new_gene)
                length = gene_dictionary[new_gene][2]
                if length > max_len:
                    max_len = length
                if length < min_len:
                    min_len = length
                new_length_list.append(length)
            
            # No. isolates
            num_sample = num_iso + len(new_sample_dict)
            row_1[3] = num_sample
            classify_cluster(num_sample, total_sample, count)
            # No. sequences
            row_1[4] = num_seq + len(cluster)
            # Avg sequences per isolate
            avg_seq = row_1[4] / row_1[3]
            row_1[5] = round(avg_seq,2)
            # Min group size nuc
            row_1[6] = min_len
            # Max group size nuc
            row_1[7] = max_len
            # Avg group size nuc
            nuc_size = (
                (float(row_1[8]) * num_seq + len(new_length_list)) 
                / row_1[4])
            row_1[8] = round(nuc_size,0)
            
            writer_1.writerow(row_1)
            cluster_id += 1


        # write new row of new clusters
        for cluster, cluster_annotation in zip(new_clusters, 
                                                new_clusters_annotation):
            row_1 = []
            sample_dict = {}
            length_list = []
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
                length = gene_dictionary[gene][2]
                length_list.append(length)
            
            # ID
            row_1.append(cluster_id)
            cluster_id += 1
            # Gene
            row_1.append(cluster_annotation[0])
            # Annotation
            row_1.append(cluster_annotation[1])
            # No. isolates
            num_sample = len(sample_dict)
            row_1.append(num_sample)
            classify_cluster(num_sample, total_sample, count)
            # No. sequences
            row_1.append(len(cluster))
            # Avg sequences per isolate
            avg_seq = len(cluster) / len(sample_dict)
            row_1.append(round(avg_seq,2))
            # Min group size nuc
            row_1.append(min(length_list))
            # Max group size nuc
            row_1.append(max(length_list))
            # Avg group size nuc
            nuc_size = sum(length_list) / len(length_list)
            row_1.append(round(nuc_size,0))
            
            writer_1.writerow(row_1)
    
    shutil.move(new_cluster_info_file, old_cluster_info_file)
    
    # write new summary file
    write_summary(count, collection_dir, total_sample)

    elapsed = datetime.now() - starttime
    logging.info(f'Update cluster info output -- time taken {str(elapsed)}')


def output_gene_info(clusters, gene_dictionary, out_dir):
    """
        Output gene_info.tsv file
        
        Parameters
        ----------
        clusters : list of list
            list of sequence IDs of each cluster
        gene dictionary : dict
            contain information of each gene
            {gene_id: (sample_id, contig, length, gene_name, gene_product)}
        out_dir : path
            output directory
        """    
    starttime = datetime.now()
    gene_file = os.path.join(out_dir, 'gene_info.tsv')
    with  open(gene_file, 'a') as fh:
        writer = csv.writer(fh, delimiter='\t')
        for i, cluster in enumerate(clusters):
            for gene_id in cluster:
                sample_id = gene_dictionary[gene_id][0]
                row = []
                row.append(gene_id)
                row.append(sample_id)
                row.append(i)
                writer.writerow(row)
    elapsed = datetime.now() - starttime
    logging.info(f'Output gene info -- time taken {str(elapsed)}')


def create_rtab(clusters, gene_dictionary, samples, out_dir):
    starttime = datetime.now()
    rtab_file = os.path.join(out_dir, 'gene_presence_absence.Rtab.gz')
    with gzip.open(rtab_file, 'wt') as fh:
        writer = csv.writer(fh, delimiter='\t')

        # write header
        header = ['ID']
        for sample in samples:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        cluster_id = 0
        for cluster in clusters:
            row = []
            # ID
            row.append(cluster_id)
            cluster_id += 1
            # Samples
            sample_dict = {}
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
            for sample in samples:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)
    elapsed = datetime.now() - starttime
    logging.info(f'Create Rtab -- time taken {str(elapsed)}')
    return rtab_file

def update_rtab(old_file, previous_clusters, new_clusters, 
                gene_dictionary, new_samples, temp_dir):
    starttime = datetime.now()
    new_file = os.path.join(temp_dir, 'gene_presence_absence.Rtab.gz')
    with gzip.open(new_file, 'wt') as out_fh, \
         gzip.open(old_file, 'rt') as in_fh:
        csv.field_size_limit(sys.maxsize)
        writer = csv.writer(out_fh, delimiter='\t')
        reader = csv.reader(in_fh, delimiter='\t')
        
        # update header
        header = next(reader)
        num_old_samples = len(header[1:])
        for sample in new_samples:
            header.append(sample['id'])
        writer.writerow(header)
        
        # update old row
        cluster_id = 0
        for row, cluster in zip(reader, previous_clusters):
            cluster_id += 1
            new_sample_dict = {}
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                new_sample_dict.setdefault(sample_id, []).append(gene)
    
            for sample in new_samples:
                sample_id = sample['id']
                gene_list = new_sample_dict.get(sample_id, [])
                row.append(len(gene_list))

            writer.writerow(row)

        # write new row
        for cluster in new_clusters:
            row = []
            # ID
            row.append(cluster_id)
            cluster_id += 1
            # Samples
            sample_dict = {}
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
            for i in range(0, num_old_samples):
                row.append(0)
            for sample in new_samples:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)

    shutil.move(new_file, old_file)
    elapsed = datetime.now() - starttime
    logging.info(f'Update Rtab -- time taken {str(elapsed)}')
    return old_file


def write_gene_dictionary(gene_dictionary, out_dir, mode='w'):
    # starttime = datetime.now()
    
    with open(os.path.join(out_dir, 'gene_dictionary.tsv'),mode) as fh:
        writer = csv.writer(fh, delimiter='\t')
        for gene in gene_dictionary:
            row = []
            row.append(gene)
            row.extend(gene_dictionary[gene])
            writer.writerow(row)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Export gene annotation -- time taken {str(elapsed)}')

def read_gene_dictionary(annotation_file):
    # starttime = datetime.now()
    
    gene_dictionary = {}
    with open(annotation_file,'r') as fh:
        csv_reader = csv.reader(fh, delimiter='\t')
        for row in csv_reader:
            row[3] = int(row[3])
            gene_dictionary[row[0]] = tuple(row[1:])
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Import gene annotation -- time taken {str(elapsed)}')
    return gene_dictionary

def write_gene_position(gene_position, out_dir):
    # starttime = datetime.now()
    
    with open(os.path.join(out_dir, 'gene_position.tsv'), 'a') as fh:
        writer = csv.writer(fh, delimiter='\t')
        for sample in gene_position:
            for seq in gene_position[sample]:
                row = []
                row.append(sample)
                row.append(seq)
                row.append(';'.join(gene_position[sample][seq]))
                writer.writerow(row)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Export gene annotation -- time taken {str(elapsed)}')

def read_gene_position(position_file):
    # starttime = datetime.now()
    
    gene_position = {}
    with open(position_file,'r') as fh:
        csv_reader = csv.reader(fh, delimiter='\t')
        for row in csv_reader:
            ls = row[1].split(';')
            gene_position[row[0]] = ls
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Import gene position -- time taken {str(elapsed)}')
    return gene_position