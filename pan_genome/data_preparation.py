import os
import re
import logging
from datetime import datetime
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def parse_gff_file(ggf_file, sample_dir, sample_id, gene_annotation, gene_position):
    bed_file = os.path.join(sample_dir, sample_id + '.bed')
    fna_file = os.path.join(sample_dir, sample_id + '.fna')
    found_fasta = False
    sample_dict = {}
    gene_position[sample_id] = sample_dict
    with open(ggf_file,'r') as in_fh, open(bed_file, 'w') as bed_fh, open(fna_file, 'w') as fna_fh:
        for line in in_fh:
            if found_fasta == True:
                fna_fh.write(line)
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
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
            gene_id = re.findall(r"ID=(.+?);",cells[8])
            gene_id = gene_id[0]
            if gene_id in gene_annotation:
                gene_id = gene_id + datetime.now().strftime("%f")
            trand = cells[6]

            # create bed file
            row = [seq_id, str(start-1), str(end), gene_id, '1', trand]
            bed_fh.write('\t'.join(row)+ '\n')

            # add to gene_annotation           
            gene_dict = {}
            gene_dict['sample_id'] = sample_id
            gene_dict['contig'] = seq_id
            gene_dict['length'] = length
            gene_name = re.findall(r"Name=(.+?);",cells[8])
            if len(gene_name) != 0:
                gene_name = gene_name[0]
                gene_dict['name'] = gene_name
            gene_product = re.findall(r"product=(.+?)$",cells[8])
            if len(gene_product) != 0:
                gene_product = gene_product[0]
                gene_dict['product'] = gene_product
            gene_annotation[gene_id] = gene_dict

            # add to gene_position
            sample_dict.setdefault(seq_id, []).append(gene_id)
    return bed_file, fna_file

def extract_proteins(samples, out_dir, gene_annotation, gene_position, timing_log=None):
    statime = datetime.now()
    for sample in samples:
        starttime = datetime.now()
        sample_id = sample['id']
        sample_dir = os.path.join(out_dir, 'samples', sample_id)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        # parse gff file
        [bed_file, fna_file] = parse_gff_file(
            ggf_file = sample['input_file'], 
            sample_dir = sample_dir, 
            sample_id = sample_id, 
            gene_annotation = gene_annotation,
            gene_position = gene_position
            )
        # extract nucleotide region
        extracted_fna_file = os.path.join(sample_dir, sample_id +'.extracted.fna')
        cmd = f"bedtools getfasta -s -fi {fna_file} -bed {bed_file} -fo {extracted_fna_file} -name > /dev/null 2>&1"
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running bedtools')
        # translate nucleotide to protein
        faa_file = os.path.join(sample_dir, sample_id +'.faa')
        translate_protein(nu_fasta=extracted_fna_file, pro_fasta=faa_file)
        elapsed = datetime.now() - starttime
        logging.info(f'Extract protein of {sample_id} -- time taken {str(elapsed)}')
        
        sample['bed'] = bed_file
        sample['extracted_fna_file'] = extracted_fna_file
        sample['faa_file'] = faa_file
        sample['sample_dir'] = sample_dir
    elapsed = datetime.now() - statime
    logging.info(f'Extract protein -- time taken {str(elapsed)}')
    return gene_annotation


def combine_proteins(out_dir, samples, timing_log=None):
    starttime = datetime.now()

    combined_faa_file = os.path.join(out_dir, 'combined.faa')
    faa_file_list =[]
    for sample in samples:
        faa_file = sample['faa_file']
        faa_file_list.append(faa_file)

    cmd = "cat {} > {}".format(" ".join(faa_file_list),combined_faa_file)
    os.system(cmd)
    
    elapsed = datetime.now() - starttime
    logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file