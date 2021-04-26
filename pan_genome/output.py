import os
import csv
import logging
from datetime import datetime
from pan_genome.utils import run_command, include_fasta

logger = logging.getLogger(__name__)


def create_spreadsheet(annotated_clusters, gene_annotation, samples, out_dir):
    starttime = datetime.now()
    spreadsheet_file = os.path.join(out_dir, 'gene_presence_absence.csv')
    with open(spreadsheet_file, 'w') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # write header
        header = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment','Order within Fragment', 'Accessory Fragment','Accessory Order with Fragment', 'QC','Min group size nuc', 'Max group size nuc', 'Avg group size nuc' ]
        for sample in samples:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        for cluster in annotated_clusters:
            this_cluster = annotated_clusters[cluster]
            sample_dict = {}
            for gene in this_cluster['gene_id']:
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
            row.append(this_cluster['product'])
            # No. isolates
            row.append(len(sample_dict))
            # No. sequences
            row.append(len(this_cluster['gene_id']))
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
                sample_id = gene_annotation[gene]['sample_id']
                if sample_id not in sample_dict:
                    sample_dict[sample_id] = []
                sample_dict[sample_id].append(gene)
            for sample in samples:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)
    elapsed = datetime.now() - starttime
    logging.info(f'Create Rtab -- time taken {str(elapsed)}')
    return rtab_file


def create_summary(split_clusters, out_dir, samples):
    starttime = datetime.now()
    num_core = 0
    num_soft_core = 0
    num_shell = 0
    num_cloud = 0
    num_sample = len(samples)
    for cluster in split_clusters:
        num = len(cluster)
        percent = num / num_sample
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
    representative_list = []
    for cluster in clusters:
        length_max = 0
        representative = None
        for gene_id in cluster:
            length = gene_annotation[gene_id]['length']
            if length > length_max:
                representative = gene_id
                length_max = length
        representative_list.append(representative)
    representative_list=set(representative_list)
    include_fasta(
        fasta_file=faa_fasta, 
        include_list=representative_list, 
        output_file=representative_fasta
        )
    elapsed = datetime.now() - starttime
    logging.info(f'Create representative fasta -- time taken {str(elapsed)}')
    return representative_fasta