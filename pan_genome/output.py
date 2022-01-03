import os
import csv
import logging
from datetime import datetime
import pandas as pd
import pan_genome.utils as utils

logger = logging.getLogger(__name__)


def create_spreadsheet(annotated_clusters, gene_annotation, samples, out_dir):
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
        for cluster in annotated_clusters:
            row = []
            sample_dict = {}
            length_list = []
            this_cluster = annotated_clusters[cluster]
            for gene in this_cluster['gene_id']:
                sample_id = gene_annotation[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
                length = gene_annotation[gene][2]
                length_list.append(length)
            
            # Gene
            row.append(cluster)
            # Annotation
            row.append(this_cluster['product'])
            # No. isolates
            row.append(len(sample_dict))
            # No. sequences
            row.append(len(this_cluster['gene_id']))
            # Avg sequences per isolate
            avg_seq = len(this_cluster['gene_id']) / len(sample_dict)
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


def create_rtab(annotated_clusters, gene_annotation, samples, out_dir):
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
        for cluster in annotated_clusters:
            row = []
            # Gene
            row.append(cluster)
            # Samples
            sample_dict = {}
            for gene in annotated_clusters[cluster]['gene_id']:
                sample_id = gene_annotation[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
            for sample in samples:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)
    elapsed = datetime.now() - starttime
    logging.info(f'Create Rtab -- time taken {str(elapsed)}')
    return rtab_file


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


def create_representative_fasta(clusters, gene_annotation, faa_fasta, out_dir):
    starttime = datetime.now()
    representative_fasta = os.path.join(out_dir, 'representative.fasta')
    representative_list = set()
    for cluster in clusters:
        length_max = 0
        representative = None
        for gene_id in cluster:
            length = gene_annotation[gene_id][2]
            if length > length_max:
                representative = gene_id
                length_max = length
        representative_list.add(representative)
    utils.create_fasta_include(
        fasta_file=faa_fasta, 
        include_list=representative_list, 
        output_file=representative_fasta
        )
    elapsed = datetime.now() - starttime
    logging.info(f'Create representative fasta -- time taken {str(elapsed)}')
    
    return representative_fasta


def export_gene_annotation(gene_annotation, out_dir):
    # starttime = datetime.now()
    
    with open(os.path.join(out_dir, 'gene_annotation.tsv'),'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        for gene in gene_annotation:
            row = []
            row.append(gene)
            row.extend(gene_annotation[gene])
            writer.writerow(row)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Export gene annotation -- time taken {str(elapsed)}')

def import_gene_annotation(annotation_file):
    # starttime = datetime.now()
    
    gene_annotation = {}
    with open(annotation_file,'r') as fh:
        csv_reader = csv.reader(fh, delimiter='\t')
        for row in csv_reader:
            row[3] = int(row[3])
            gene_annotation[row[0]] = tuple(row[1:])
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Import gene annotation -- time taken {str(elapsed)}')
    return gene_annotation