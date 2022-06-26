import os
import re
import logging
import multiprocessing
from functools import partial
from datetime import datetime
import gzip

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pan_genome.utils as utils

logger = logging.getLogger(__name__)

def parse_gff_file(ggf_file, sample_dir, sample_id):
    bed_file = os.path.join(sample_dir, sample_id + '.bed')
    assembly_file = os.path.join(sample_dir, sample_id + '.fasta')
    gene_dictionary = {}
    found_fasta = False
    suffix = 1

    if ggf_file.endswith('.gz'):
        in_fh = gzip.open(ggf_file,'rt')
    else:
        in_fh = open(ggf_file,'r')
    
    with open(bed_file, 'w') as bed_fh, open(assembly_file, 'w') as fna_fh:
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
                    gene_name = re.sub(r'_\d+$', '', gene_name)
                    continue
                
                product = re.match(r"^product=(.+)", tag)
                if product != None:
                    gene_product = product.group(1)
            if gene_id == None:
                continue
            
            # if re.match(sample_id, gene_id) == None:
            #     gene_id = sample_id + '_' + gene_id
            if gene_id in gene_dictionary:
                logging.info(f'{gene_id} already exists -- add suffix')
                gene_id += '_{:05d}'.format(suffix)
                suffix += 1
            
            # create bed file
            row = [seq_id, str(start-1), str(end), gene_id, '1', trand]
            bed_fh.write('\t'.join(row)+ '\n')
            # add to gene_dictionary           
            gene_dictionary[gene_id] = (
                sample_id, seq_id, length, gene_name, gene_product)
    
    in_fh.close()
    
    return bed_file, assembly_file, gene_dictionary
    

def process_single_sample_gff(sample, out_dir, table):
    sample_id = sample['id']
    sample_dir = os.path.join(out_dir, 'samples', sample_id)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    
    # parse gff file
    bed_file, assembly_file, gene_dictionary = parse_gff_file(
        ggf_file = sample['gff_file'], 
        sample_dir = sample_dir, 
        sample_id = sample_id
        )
    if sample['assembly'] != None:
        assembly_file = sample['assembly']
    
    # extract nucleotide region
    fna_file = os.path.join(sample_dir, sample_id +'.fna')
    cmd = (f"bedtools getfasta -s -fi {assembly_file} -bed {bed_file} "
           f"-fo {fna_file} -name > /dev/null 2>&1")
    utils.run_command(cmd)

    # translate nucleotide to protein
    faa_file = os.path.join(sample_dir, sample_id +'.faa')
    utils.translate_protein(nu_fasta=fna_file, pro_fasta=faa_file, table=table)
    
    if sample['assembly'] == None: # do not remove input assembly
        os.remove(assembly_file)
    os.remove(bed_file)
    os.remove(assembly_file + '.fai')
    os.remove(fna_file)

    return gene_dictionary


def process_single_sample_fasta(sample, out_dir, table):

    sample_id = sample['id']
    sample_dir = os.path.join(out_dir, 'samples', sample_id)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)

    # gene prediction
    assembly_file = sample['assembly']
    faa_file = os.path.join(sample_dir, sample_id +'.original.faa')
    
    if assembly_file.endswith('.gz'):
        cmd = (f'zcat {assembly_file} | '
               f'prodigal -a {faa_file} -g {table} -c -m -q -o /dev/null')
    else:
        cmd = (f'prodigal -i {assembly_file} -a {faa_file} -g {table} '
               f'-c -m -q -o /dev/null')
    
    if not os.path.isfile(faa_file):
        utils.run_command(cmd)
                    

    # change gene id and extract coordinates
    gene_dictionary = {}
    rewrite_faa_file = os.path.join(sample_dir, sample_id +'.faa')
    count = 1
    passed_genes = set()
    with open(faa_file, 'r') as in_fh, open(rewrite_faa_file, 'w') as out_fh:
        for record in SeqIO.parse(in_fh, "fasta"):
            contig = record.id.rsplit('_', 1)[0] 
            cells = record.description.split(' # ')
            pro = str(record.seq)
            
            # filter short genes
            start = int(cells[1])
            end = int(cells[2])
            trand = cells[3]
            length = end - start + 1
            if length < 120:
                # logger.info('Short gene')
                continue

            # filter seq with premature codon
            results = re.findall(r'\*', pro)
            if len(results) > 1:
                # logger.info('Have premature codon')
                continue

            # filter seq lacking start and stop codon
            if pro[0] != 'M' and pro[-1] != '*':
                # logger.info('Lack both start and stop codon')
                continue
            
            # filter seq which has more than 5% of unknown
            results = re.findall(r'X', pro)
            if len(results) / len (pro) > 0.05:
                # logger.info('Too many unknowns')
                continue
            
            gene_id = sample_id + '_{:05d}'.format(count)
            count += 1
            desc = '{}~~~{}~~~{}~~~{}'.format(contig, start, end, trand)
            new_record = SeqRecord(
                record.seq, id = gene_id, description = desc)
            SeqIO.write(new_record, out_fh, 'fasta')
            passed_genes.add(gene_id)
            # add to gene_dictionary           
            gene_dictionary[gene_id] = (sample_id, contig, length)

    return gene_dictionary


def extract_proteins(samples, out_dir, args, timing_log):
    starttime = datetime.now()
    
    if args.threads == 0:
        threads = multiprocessing.cpu_count()
    else:
        threads = args.threads

    with multiprocessing.Pool(processes=threads) as pool:
        if args.fasta == None:
            results = pool.map(
                partial(process_single_sample_gff, out_dir=out_dir, 
                        table=args.table, timing_log=timing_log), 
                samples)
        else:
            results = pool.map(
                partial(process_single_sample_fasta, out_dir=out_dir, 
                        table=args.table, timing_log=timing_log), 
                samples)
    
    gene_dictionary = {}
    for sample, result in zip(samples, results):
        # gene_dictionary.update(result[0])
        for k, v in result.items():
            if k in gene_dictionary:
                logging.info(f'{k} already exists -- add prefix')
                k = sample['id'] + '_' + k
            gene_dictionary[k] = v
        
    elapsed = datetime.now() - starttime
    logging.info(f'Extract protein -- time taken {str(elapsed)}')

    return gene_dictionary


def combine_proteins(collection_dir, out_dir, samples, timing_log):
    # starttime = datetime.now()

    combined_faa_file = os.path.join(out_dir, 'combined.faa')
    protein_files = os.path.join(out_dir, 'protein.txt')
    with open(protein_files, 'w') as fh:
        for sample in samples:
            sample_id = sample['id']
            faa_file = os.path.join(
                collection_dir, 'samples', sample_id, sample_id + '.faa')
            if os.path.isfile(faa_file):
                fh.write(faa_file + '\n')
            else:
                raise Exception(f'{faa_file} does not exist')

    cmd = f"cat {protein_files} | xargs cat > {combined_faa_file}"
    utils.run_command(cmd, timing_log)
            
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file