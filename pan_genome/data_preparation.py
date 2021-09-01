import os
import re
import logging
import multiprocessing
from functools import partial
from datetime import datetime
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def parse_gff_file(ggf_file, sample_dir, sample_id):
    bed_file = os.path.join(sample_dir, sample_id + '.bed')
    fna_file = os.path.join(sample_dir, sample_id + '.fna')
    gene_annotation = {}
    gene_position = {}
    found_fasta = False
    # suffix = 1
    with open(ggf_file,'r') as in_fh, open(bed_file, 'w') as bed_fh, open(fna_file, 'w') as fna_fh:
        for line in in_fh:
            if found_fasta == True:
                fna_fh.write(line)
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
                continue
            if re.match(r"^#", line) != None:
                continue
            line = line.rstrip('\n')
            cells = line.split('\t')
            if cells[2] != 'CDS':
                continue
            
            start = int(cells[3])
            end = int(cells[4])
            length = end - start
            if length < 120:
                continue
            seq_id = cells[0]
            trand = cells[6]
            tags = cells[8].split(';')
            gene_id = None
            gene_name = ''
            gene_product = ''
            for tag in tags:
                ID = re.match(r"^ID=(.+)", tag)
                if ID != None:
                    gene_id = ID.group(1)
                    # gene_id = re.sub(r'\W', '_', gene_id)
                    continue

                gene = re.match(r"^gene=(.+)", tag)
                if gene != None:
                    gene_name = gene.group(1)
                    # gene_name = re.sub(r'\W', '_', gene_name)
                    continue
                
                product = re.match(r"^product=(.+)", tag)
                if product != None:
                    gene_product = product.group(1)
            if gene_id == None:
                continue
            
            # if re.match(sample_id, gene_id) == None:
            #     gene_id = sample_id + '_' + gene_id
            if gene_id in gene_annotation:
                raise Exception(f'{gene_id} of {sample_id} appear the second time. Please fix gff files')
                # logging.info(f'{gene_id} already exists -- add suffix')
                # gene_id += '_{:05d}'.format(suffix)
                # suffix += 1
            
            # create bed file
            row = [seq_id, str(start-1), str(end), gene_id, '1', trand]
            bed_fh.write('\t'.join(row)+ '\n')
            # add to gene_annotation           
            gene_annotation[gene_id] = (sample_id, seq_id, length, gene_name, gene_product)
            # add to gene_position
            gene_position.setdefault(seq_id, []).append(gene_id)

    return bed_file, fna_file, gene_annotation, gene_position
    

def process_single_sample(sample, out_dir, table):
    # starttime = datetime.now()
    
    sample_id = sample['id']
    sample_dir = os.path.join(out_dir, 'samples', sample_id)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    
    # parse gff file
    bed_file, fna_file, gene_annotation, gene_position = parse_gff_file(
        ggf_file = sample['gff_file'], 
        sample_dir = sample_dir, 
        sample_id = sample_id
        )
    
    # extract nucleotide region
    extracted_fna_file = os.path.join(sample_dir, sample_id +'.extracted.fna')
    os.system(f"bedtools getfasta -s -fi {fna_file} -bed {bed_file} -fo {extracted_fna_file} -name > /dev/null 2>&1")
    
    # translate nucleotide to protein
    faa_file = os.path.join(sample_dir, sample_id +'.faa')
    translate_protein(nu_fasta=extracted_fna_file, pro_fasta=faa_file, table=table)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Extract protein of {sample_id} -- time taken {str(elapsed)}')

    return gene_annotation, faa_file, sample_dir, gene_position, extracted_fna_file


def extract_proteins(samples, out_dir, gene_annotation, gene_position, table, threads):
    starttime = datetime.now()
    
    if threads == 0:
        threads = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(process_single_sample, out_dir=out_dir, table=table), samples)
    
    for sample, result in zip(samples, results):
        gene_annotation.update(result[0])
        sample['faa_file'] = result[1]
        sample['sample_dir'] = result[2]
        gene_position[sample['id']] = result[3]
        sample['fna_file'] = result[4]
    
    elapsed = datetime.now() - starttime
    logging.info(f'Extract protein -- time taken {str(elapsed)}')


def combine_proteins(out_dir, samples, timing_log=None):
    # starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    combined_faa_file = os.path.join(out_dir, 'combined.faa')
    faa_file_list =[]
    for sample in samples:
        faa_file = sample['faa_file']
        faa_file_list.append(faa_file)

    cmd = "cat {} > {}".format(" ".join(faa_file_list),combined_faa_file)
    os.system(cmd)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file