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
    assembly_file = os.path.join(sample_dir, sample_id + '.fasta')
    gene_annotation = {}
    gene_position = {}
    found_fasta = False
    suffix = 1
    with open(ggf_file,'r') as in_fh, open(bed_file, 'w') as bed_fh, open(assembly_file, 'w') as fna_fh:
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
            length = end - start + 1
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
                    gene_name = re.sub(r'\W', '_', gene_name)
                    continue
                
                product = re.match(r"^product=(.+)", tag)
                if product != None:
                    gene_product = product.group(1)
            if gene_id == None:
                continue
            
            # if re.match(sample_id, gene_id) == None:
            #     gene_id = sample_id + '_' + gene_id
            if gene_id in gene_annotation:
                # raise Exception(f'{gene_id} of {sample_id} appear the second time. Please fix gff files')
                logging.info(f'{gene_id} already exists -- add suffix')
                gene_id += '_{:05d}'.format(suffix)
                suffix += 1
            
            # create bed file
            row = [seq_id, str(start-1), str(end), gene_id, '1', trand]
            bed_fh.write('\t'.join(row)+ '\n')
            # add to gene_annotation           
            gene_annotation[gene_id] = (sample_id, seq_id, length, gene_name, gene_product)
            # add to gene_position
            gene_position.setdefault(seq_id, []).append(gene_id)

    return bed_file, assembly_file, gene_annotation, gene_position
    

def process_single_sample(sample, out_dir, table):
    # starttime = datetime.now()
    
    sample_id = sample['id']
    sample_dir = os.path.join(out_dir, 'samples', sample_id)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    
    # parse gff file
    bed_file, assembly_file, gene_annotation, gene_position = parse_gff_file(
        ggf_file = sample['gff_file'], 
        sample_dir = sample_dir, 
        sample_id = sample_id
        )
    if sample['assembly'] != None:
        assembly_file = sample['assembly']
    
    # extract nucleotide region
    fna_file = os.path.join(sample_dir, sample_id +'.fna')
    os.system(f"bedtools getfasta -s -fi {assembly_file} -bed {bed_file} -fo {fna_file} -name > /dev/null 2>&1")
    
    # translate nucleotide to protein
    faa_file = os.path.join(sample_dir, sample_id +'.faa')
    translate_protein(nu_fasta=fna_file, pro_fasta=faa_file, table=table)
    
    if sample['assembly'] == None:
        os.remove(assembly_file)
    os.remove(bed_file)
    os.remove(assembly_file + '.fai')

    # elapsed = datetime.now() - starttime
    # logging.info(f'Extract protein of {sample_id} -- time taken {str(elapsed)}')

    return gene_annotation, gene_position


def extract_proteins(samples, out_dir, gene_annotation, gene_position, table, threads):
    starttime = datetime.now()
    
    if threads == 0:
        threads = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(process_single_sample, out_dir=out_dir, table=table), samples)
    
    for sample, result in zip(samples, results):
        # gene_annotation.update(result[0])
        for k, v in result[0].items():
            if k in gene_annotation:
                logging.info(f'{k} already exists -- add prefix')
                k = sample['id'] + '_' + k
            gene_annotation[k] = v
        gene_position[sample['id']] = result[1]
        
    elapsed = datetime.now() - starttime
    logging.info(f'Extract protein -- time taken {str(elapsed)}')


def combine_proteins(collection_dir, out_dir, samples):
    # starttime = datetime.now()

    combined_faa_file = os.path.join(out_dir, 'combined.faa')
    protein_files = os.path.join(out_dir, 'protein.txt')
    with open(protein_files, 'w') as fh:
        for sample in samples:
            sample_id = sample['id']
            faa_file = os.path.join(collection_dir, 'samples', sample_id, sample_id + '.faa')
            if os.path.isfile(faa_file):
                fh.write(faa_file + '\n')
            else:
                raise Exception(f'{faa_file} does not exist')

    cmd = f"cat {protein_files} | xargs cat > {combined_faa_file}"
    os.system(cmd)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file