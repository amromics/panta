import os
import csv
import logging
import shutil
from datetime import datetime
import pandas as pd
import pan_genome.utils as utils

logger = logging.getLogger(__name__)


def create_spreadsheet(clusters, clusters_annotation, gene_dictionary, samples, out_dir):
    starttime = datetime.now()
    spreadsheet_file = os.path.join(out_dir, 'gene_presence_absence.csv')
    with open(spreadsheet_file, 'w') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)

        # write header
        header = ['Gene', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc' ]
        for sample in samples:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        for cluster, cluster_annotation in zip(clusters, clusters_annotation):
            row = []
            sample_dict = {}
            length_list = []
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
                length = gene_dictionary[gene][2]
                length_list.append(length)
            
            # Gene
            row.append(cluster_annotation[0])
            # Annotation
            row.append(cluster_annotation[1])
            # No. isolates
            row.append(len(sample_dict))
            # No. sequences
            row.append(len(cluster))
            # Avg sequences per isolate
            avg_seq = len(cluster) / len(sample_dict)
            row.append(round(avg_seq,2))
            # Min group size nuc
            row.append(min(length_list))
            # Max group size nuc
            row.append(max(length_list))
            # Avg group size nuc
            nuc_size = sum(length_list) / len(length_list)
            row.append(round(nuc_size,0))

            # sample columns
            for sample in samples:
                sample_id = sample['id']
                if sample_id in sample_dict:
                    gene_list = sample_dict[sample_id]
                    row.append('\t'.join(gene_list))
                else:
                    row.append('')
            writer.writerow(row)
    elapsed = datetime.now() - starttime
    logging.info(f'Create spreadsheet -- time taken {str(elapsed)}')
    return spreadsheet_file


def update_spreadsheet(old_file, old_clusters, new_clusters, new_clusters_annotation, gene_dictionary, new_samples, all_samples, temp_dir):
    starttime = datetime.now()
    new_file = os.path.join(temp_dir, 'gene_presence_absence.csv')
    with open(new_file, 'w') as out_fh, open(old_file, 'r') as in_fh:
        writer = csv.writer(out_fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        reader = csv.reader(in_fh, delimiter=',')
        
        # update header
        header = next(reader)
        for sample in new_samples:
            header.append(sample['id'])
        writer.writerow(header)
        
        for row, cluster in zip(reader, old_clusters):
            num_seq = int(row[3])
            num_iso = int(row[2])
            min_len = int(row[5])
            max_len = int(row[6])
            new_sample_dict = {}
            new_length_list = []
            for i,gene in enumerate(cluster):
                if i < num_seq:
                    continue
                sample_id = gene_dictionary[gene][0]
                new_sample_dict.setdefault(sample_id, []).append(gene)
                length = gene_dictionary[gene][2]
                if length > max_len:
                    max_len = length
                if length < min_len:
                    min_len = length
                new_length_list.append(length)
            
            
            # No. isolates
            row[2] = num_iso + len(new_sample_dict)
            # No. sequences
            row[3] = len(cluster)
            # Avg sequences per isolate
            avg_seq = row[3] / row[2]
            row[4] = round(avg_seq,2)
            # Min group size nuc
            row[5] =min_len
            # Max group size nuc
            row[6] = max_len
            # Avg group size nuc
            nuc_size = (float(row[7]) * num_seq + len(new_length_list)) / len(cluster)
            row[7] = round(nuc_size,0)

            for sample in new_samples:
                sample_id = sample['id']
                if sample_id in new_sample_dict:
                    gene_list = new_sample_dict[sample_id]
                    row.append('\t'.join(gene_list))
                else:
                    row.append('')
            
            writer.writerow(row)


        # write new row
        for cluster, cluster_annotation in zip(new_clusters, new_clusters_annotation):
            row = []
            sample_dict = {}
            length_list = []
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
                length = gene_dictionary[gene][2]
                length_list.append(length)
            
            # Gene
            row.append(cluster_annotation[0])
            # Annotation
            row.append(cluster_annotation[1])
            # No. isolates
            row.append(len(sample_dict))
            # No. sequences
            row.append(len(cluster))
            # Avg sequences per isolate
            avg_seq = len(cluster) / len(sample_dict)
            row.append(round(avg_seq,2))
            # Min group size nuc
            row.append(min(length_list))
            # Max group size nuc
            row.append(max(length_list))
            # Avg group size nuc
            nuc_size = sum(length_list) / len(length_list)
            row.append(round(nuc_size,0))

            # sample columns
            for sample in all_samples:
                sample_id = sample['id']
                if sample_id in sample_dict:
                    gene_list = sample_dict[sample_id]
                    row.append('\t'.join(gene_list))
                else:
                    row.append('')
            writer.writerow(row)
    
    shutil.move(new_file, old_file)
    elapsed = datetime.now() - starttime
    logging.info(f'Create spreadsheet -- time taken {str(elapsed)}')


def create_rtab(clusters, clusters_annotation, gene_dictionary, samples, out_dir):
    starttime = datetime.now()
    rtab_file = os.path.join(out_dir, 'gene_presence_absence.Rtab')
    with open(rtab_file, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')

        # write header
        header = ['Gene']
        for sample in samples:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        for cluster, cluster_annotation in zip(clusters, clusters_annotation):
            row = []
            # Gene
            row.append(cluster_annotation[0])
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

def update_rtab(old_file, old_clusters, new_clusters, new_clusters_annotation, gene_dictionary, new_samples, all_samples, temp_dir):
    starttime = datetime.now()
    new_file = os.path.join(temp_dir, 'gene_presence_absence.Rtab')
    with open(new_file, 'w') as out_fh, open(old_file, 'r') as in_fh:
        writer = csv.writer(out_fh, delimiter='\t')
        reader = csv.reader(in_fh, delimiter='\t')
        
        # update header
        header = next(reader)
        for sample in new_samples:
            header.append(sample['id'])
        writer.writerow(header)
  
        for row, cluster in zip(reader, old_clusters):
            num_seq = 0
            for cell in row[1:]:
                num_seq += int(cell)
            new_sample_dict = {}
            for i,gene in enumerate(cluster):
                if i < num_seq:
                    continue
                sample_id = gene_dictionary[gene][0]
                new_sample_dict.setdefault(sample_id, []).append(gene)
    
            for sample in new_samples:
                sample_id = sample['id']
                gene_list = new_sample_dict.get(sample_id, [])
                row.append(len(gene_list))

            writer.writerow(row)

        # write row
        for cluster, new_cluster_annotation in zip(new_clusters,new_clusters_annotation):
            row = []
            # Gene
            row.append(new_cluster_annotation[0])
            # Samples
            sample_dict = {}
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
            for sample in all_samples:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)

    shutil.move(new_file, old_file)
    elapsed = datetime.now() - starttime
    logging.info(f'Create Rtab -- time taken {str(elapsed)}')
    return old_file


def create_summary(rtab_file, out_dir):
    starttime = datetime.now()
    num_core = 0
    num_soft_core = 0
    num_shell = 0
    num_cloud = 0
    with open(rtab_file, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            cells = line.split('\t')
            if cells[0] == 'Gene':
                num_sample = len(cells) - 1
                continue
            num_zero = 0
            for cell in cells:
                if cell == '0':
                    num_zero += 1
            percent = (num_sample - num_zero) / num_sample
            if percent >= 0.99:
                num_core += 1
            elif percent >= 0.95:
                num_soft_core += 1
            elif percent >= 0.15:
                num_shell += 1
            else:
                num_cloud += 1
    total = num_core + num_soft_core + num_shell + num_cloud

    summary_file = os.path.join(out_dir, 'summary_statistics.txt')
    with open(summary_file, 'w') as fh:
        fh.write('Core genes' + '\t' + '(99% <= strains <= 100%)' + '\t'+ str(num_core) + '\n')
        fh.write('Soft core genes' + '\t' + '(95% <= strains < 99%)' + '\t'+ str(num_soft_core) + '\n')
        fh.write('Shell genes' + '\t' + '(15% <= strains < 95%)' + '\t' + str(num_shell) + '\n')
        fh.write('Cloud genes' + '\t' + '(0% <= strains < 15%)' + '\t'+ str(num_cloud) + '\n')
        fh.write('Total genes' + '\t' + '(0% <= strains <= 100%)' + '\t'+ str(total))
    elapsed = datetime.now() - starttime
    logging.info(f'Create summary -- time taken {str(elapsed)}')
    return summary_file


# def create_representative_fasta(clusters, gene_dictionary, fasta_list, out_dir):
#     starttime = datetime.now()
#     representative_fasta = os.path.join(out_dir, 'representative.fasta')
#     representative_list = set()
#     for cluster in clusters:
#         length_max = 0
#         representative = None
#         for gene_id in cluster:
#             length = gene_dictionary[gene_id][2]
#             if length > length_max:
#                 representative = gene_id
#                 length_max = length
#         representative_list.add(representative)
#     utils.create_fasta_include(
#         fasta_file_list=fasta_list, 
#         include_list=representative_list, 
#         output_file=representative_fasta
#         )
#     elapsed = datetime.now() - starttime
#     logging.info(f'Create representative fasta -- time taken {str(elapsed)}')
    
#     return representative_fasta


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

def write_gene_position(gene_position, out_dir, mode='w'):
    # starttime = datetime.now()
    
    with open(os.path.join(out_dir, 'gene_position.tsv'),mode) as fh:
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