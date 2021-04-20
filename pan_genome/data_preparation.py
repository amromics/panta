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

def parse_gff_file(ggf_file, bed_out_file, fasta_out_file, sample_id, dictionary):
    """
    Parse gff file, filter out small genes(<18 base). Create a bed file, 
    fasta file and a dictionary contain gene annotation. 

    Parameters
    -------
    -------
    """
    found_fasta = 0
    last_seq_id = None
    count = 1
    with open(ggf_file,'r') as in_fh, open(bed_out_file, 'w') as bed_fh, open(fasta_out_file, 'w') as fasta_fh:
        for line in in_fh:
            if found_fasta == 1:
                fasta_fh.write(line)
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = 1
                continue
            if re.match(r"^##", line) != None:
                continue
            line = line.rstrip('\n')
            cells = line.split('\t')
            if cells[2] != 'CDS':
                continue
                
            #filter out small genes (<18 base)
            start = int(cells[3])
            end = int(cells[4])
            length = end - start
            if length < 18:
                continue

            seq_id = cells[0]
            
            #if seq_id != last_seq_id:
            #    count = 1
            #    last_seq_id = seq_id
            #gene_id = sample_id + '_' + seq_id + '_' + str(count)
            #count = count + 1
            
            gene_id = re.findall(r"ID=(.+?);",cells[8])
            gene_id = gene_id[0]
            #if gene_id in dictionary:
            #    gene_id = gene_id + datetime.now().strftime("%M%S%f")
            
            trand = cells[6]
            row = [seq_id, str(start -1), str(end), gene_id, '1', trand]
            bed_fh.write('\t'.join(row)+ '\n')

            # create annotation dictionary
            dictionary[gene_id] = {}
            dictionary[gene_id]['sample_id'] = sample_id
            dictionary[gene_id]['length'] = length
            gene_name = re.findall(r"Name=(.+?);",cells[8])
            if len(gene_name) != 0:
                gene_name = gene_name[0]
                dictionary[gene_id]['name'] = gene_name
            gene_product = re.findall(r"product=(.+?)$",cells[8])
            if len(gene_product) != 0:
                gene_product = gene_product[0]
                dictionary[gene_id]['product'] = gene_product
            dictionary[gene_id]['contig'] = seq_id


def extract_proteins(report, gene_annotation, timing_log=None):
    """
    Take in GFF file and create protein sequences in FASTA format.
    Create json file contain annotation information of each gene_id

    Parameters
    -------
    -------
    """
    for sample in report['samples']:
        sample_id = sample['id']
        sample_dir = os.path.join(report['pan_genome'], 'samples', sample_id)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        # parse gff file
        ggf_file = sample['input_file']
        fna_file = os.path.join(sample_dir, sample_id + '.fna')
        bed_file = os.path.join(sample_dir, sample_id + '.bed')
        parse_gff_file(ggf_file, bed_file, fna_file, sample_id, gene_annotation)
        # extract nucleotide region
        extracted_fna_file = os.path.join(sample_dir, sample_id +'.extracted.fna')
        cmd = f"bedtools getfasta -s -fi {fna_file} -bed {bed_file} -fo {extracted_fna_file} -name > /dev/null 2>&1"
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running bedtools')
        # translate nucleotide to protein
        faa_file = os.path.join(sample_dir, sample_id +'.faa')
        with open(faa_file, 'w') as faa_hd:
            for seq_record in SeqIO.parse(extracted_fna_file, "fasta"):
                headers = seq_record.id.split(':')
                seq_record.id = headers[0]
                seq_record.name = ''
                seq_record.description = ''
                seq_record.seq = seq_record.seq.translate(table=11)
                SeqIO.write(seq_record, faa_hd, 'fasta')
        sample['bed'] = bed_file
        sample['extracted_fna_file'] = extracted_fna_file
        sample['faa_file'] = faa_file
    
    report['gene_annotation'] = gene_annotation
    return report


def combine_proteins(report, timing_log=None):
    """
    Take in multiple FASTA sequences containing proteomes and concat them together 
    and output a FASTA file

    Parameters
    -------
    -------
    """
    temp_dir = report['temp_dir']
    combined_faa_file = os.path.join(temp_dir, 'combined.faa')
    faa_file_list =[]
    for sample in report['samples']:
        faa_file = sample['faa_file']
        faa_file_list.append(faa_file)

    cmd = "cat {} > {}".format(" ".join(faa_file_list),combined_faa_file)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error combining protein sequences')
    
    report['combined_faa_file'] = combined_faa_file
    return report