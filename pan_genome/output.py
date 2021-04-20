import os
import shutil
import re
import json
import gzip
import csv
import logging
from Bio import SeqIO
import pandas as pd
from pan_genome.utils import run_command

logger = logging.getLogger(__name__)


def create_spreadsheet(report):
    spreadsheet_file = os.path.join(report['pan_genome'], 'gene_presence_absence.csv')
    annotated_clusters = report['annotated_clusters']
    gene_annotation = report['gene_annotation']
    with open(spreadsheet_file, 'w') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # write header
        header = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment','Order within Fragment', 'Accessory Fragment','Accessory Order with Fragment', 'QC','Min group size nuc', 'Max group size nuc', 'Avg group size nuc' ]
        for sample in report['samples']:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        for cluster in annotated_clusters:
            sample_dict = {}
            for gene in annotated_clusters[cluster]['gene_id']:
                sample_id = gene_annotation[gene]['sample_id']
                if sample_id not in sample_dict:
                    sample_dict[sample_id] = []
                sample_dict[sample_id].append(gene)
            
            row = []
            # Gene
            row.append(cluster)
            # Non-unique Gene name
            row.append("")
            # Annotation
            row.append(annotated_clusters[cluster]['product'])
            # No. isolates
            row.append(len(sample_dict))
            # No. sequences
            row.append(len(annotated_clusters[cluster]['gene_id']))
            # Avg sequences per isolate
            row.append("")
            # Genome Fragment
            row.append("")
            # Order within Fragment
            row.append("")
            # Accessory Fragment
            row.append("")
            # Accessory Order with Fragment
            row.append("")
            # QC
            row.append("")
            # Min group size nuc
            row.append("")
            # Max group size nuc
            row.append("")
            # Avg group size nuc
            row.append("")
            # sample columns
            for sample in report['samples']:
                sample_id = sample['id']
                if sample_id in sample_dict:
                    gene_list = sample_dict[sample_id]
                    row.append('\t'.join(gene_list))
                else:
                    row.append('')
            writer.writerow(row)
    report['spreadsheet'] = spreadsheet_file
    return report


def create_rtab(report):
    rtab_file = os.path.join(report['pan_genome'], 'gene_presence_absence.Rtab')
    annotated_clusters = report['annotated_clusters']
    gene_annotation = report['gene_annotation']
    with open(rtab_file, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')

        # write header
        header = ['Gene']
        for sample in report['samples']:
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
                sample_id = gene_annotation[gene]['sample_id']
                if sample_id not in sample_dict:
                    sample_dict[sample_id] = []
                sample_dict[sample_id].append(gene)
            for sample in report['samples']:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)
    report['rtab'] = rtab_file
    return report


def create_summary(report):
    rtab_file = report['rtab']
    cluster_df = pd.read_csv(rtab_file, sep='\t', index_col='Gene')

    num_core = 0
    num_soft_core = 0
    num_shell = 0
    num_cloud = 0
    num_sample = len(report['samples'])
    for cluster, row in cluster_df.iterrows():
        absent = len(row[row == 0])
        percent = (num_sample - absent) / num_sample
        if percent >= 0.99:
            num_core += 1
        elif percent >= 0.95:
            num_soft_core += 1
        elif percent >= 0.15:
            num_shell += 1
        else:
            num_cloud += 1
    total = num_core + num_soft_core + num_shell + num_cloud

    summary_file = os.path.join(report['pan_genome'], 'summary_statistics.txt')
    with open(summary_file, 'w') as fh:
        fh.write('Core genes' + '\t' + '(99% <= strains <= 100%)' + '\t'+ str(num_core) + '\n')
        fh.write('Soft core genes' + '\t' + '(95% <= strains < 99%)' + '\t'+ str(num_soft_core) + '\n')
        fh.write('Shell genes' + '\t' + '(15% <= strains < 95%)' + '\t' + str(num_shell) + '\n')
        fh.write('Cloud genes' + '\t' + '(0% <= strains < 15%)' + '\t'+ str(num_cloud) + '\n')
        fh.write('Total genes' + '\t' + '(0% <= strains <= 100%)' + '\t'+ str(total))
    report['summary'] = summary_file
    return report


def create_representative_fasta(report):
    unsplit_clusters = report['inflated_unsplit_clusters']
    gene_annotation = report['gene_annotation']
    combined_fasta = report['combined_faa_file']
    representative_fasta = os.path.join(report['pan_genome'], 'representative.fasta')

    representative_list = []
    for cluster in unsplit_clusters:
        length_max = 0
        representative = None
        for gene_id in cluster:
            length = gene_annotation[gene_id]['length']
            if length > length_max:
                representative = gene_id
                length_max = length
        representative_list.append(representative)
    
    with open(representative_fasta, 'w') as fh:
        for seq_record in SeqIO.parse(combined_fasta, 'fasta'):
            if seq_record.id in representative_list:
                SeqIO.write(seq_record, fh, 'fasta')
    
    report['representative_fasta'] = representative_fasta
    return report