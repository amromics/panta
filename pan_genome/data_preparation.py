import os
import re
import logging
import multiprocessing
from functools import partial
from datetime import datetime
import gzip
from Bio import SeqIO, Seq
from pan_genome.utils import *
from collections import OrderedDict

logger = logging.getLogger(__name__)

def parse_gff_file(ggf_file, sample_dir, sample_id):
    assembly_file = os.path.join(sample_dir, sample_id + '.fasta')
    gene_annotation = OrderedDict()
    gene_position = OrderedDict()
    found_fasta = False
    suffix = 1
    bed_records = []
    if ggf_file.endswith('gff.gz'):
        in_fh = gzip.open(ggf_file,'rt')
    else: #end with .gff
        in_fh = open(ggf_file)

    with open(assembly_file, 'w') as fna_fh:
        gene_index = 0
        seq_id = None
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
            #print(line)
            cells = line.split('\t')
            if cells[2] != 'CDS':
                continue

            start = int(cells[3])
            end = int(cells[4])
            length = end - start + 1            
            # if length < 120:
            #     continue
            if length % 3 != 0:
                continue
            cells[0] = cells[0].replace('-','_') #make sure seq_id has no -
            if seq_id != cells[0]:
                seq_id = cells[0]
                gene_index = 0

            strand = cells[6]
            tags = cells[8].split(';')
            gene_id = None
            gene_name = ''
            gene_product = ''
            for tag in tags:
                ID = re.match(r"^ID=(.+)", tag)
                if ID != None:
                    gene_id = ID.group(1)                    
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

            # Ensure gene_id is in the format of sample_id-seq_id-gene_tag
            if not gene_id.startswith(sample_id + '-'):
                gene_id = sample_id + '-' + gene_id

            if not gene_id.startswith(sample_id + '-' + seq_id + '-'):
                gene_id = sample_id + '-' + seq_id + '-' + gene_id[len(sample_id)+1:]
            
            suffix += 1            
            bed_records.append((seq_id, start - 1, end, gene_id, strand))            
            # add to gene_annotation
            gene_annotation[gene_id] = (sample_id, seq_id, length, gene_name, gene_product, gene_index, strand)
            gene_index += 1
            # add to gene_position
            gene_position.setdefault(seq_id, []).append(gene_id)
    in_fh.close()
    return assembly_file, gene_annotation, gene_position, bed_records


def process_single_sample(sample, out_dir, table):
    sample_id = sample['id']
    sample_dir = os.path.join(out_dir, 'samples', sample_id)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)

    # parse gff file
    assembly_file, gene_annotation, gene_position, bed_records = parse_gff_file(
        ggf_file = sample['gff_file'],
        sample_dir = sample_dir,
        sample_id = sample_id
        )
    if sample['assembly'] != None:
        assembly_file = sample['assembly']

    fna_file = os.path.join(sample_dir, sample_id +'.fna.gz')
    faa_file = os.path.join(sample_dir, sample_id +'.faa.gz')
    seqs = {}
    for seq in SeqIO.parse(assembly_file, 'fasta'):
        seqs[seq.id.replace('-','_')] = seq
    os.remove(assembly_file)

    gene_seqs = []
    protein_seqs = []    
    for bed_record in bed_records:
        (seq_id, start, end, gene_id, strand) = bed_record
        if seq_id not in seqs:
            continue        
        seq = seqs[seq_id]  
        if end > len(seq):
            continue

        subseq = Seq(str(seq.seq)[start:end])
        if strand == '-':
            subseq = subseq.reverse_complement()        
        gene_seq = SeqIO.SeqRecord(subseq, gene_id)
        gene_seqs.append(gene_seq)
        protein_seq = gene_seq.translate(table=table, stop_symbol='') 
        protein_seq.id = gene_id   
        protein_seqs.append(protein_seq)
    with gzip.open(fna_file, 'wt') as fna_fh:
        SeqIO.write(gene_seqs, fna_fh, 'fasta')
    with gzip.open(faa_file, 'wt') as faa_fh:
        SeqIO.write(protein_seqs, faa_fh, 'fasta')
    
    annotation_fn = os.path.join(sample_dir, sample_id +'.annotation')    
    with open(annotation_fn, 'w') as ga_fp:
        #no header
        for gene_id in gene_annotation:
                (sample_id, seq_id, length, gene_name, gene_product, gene_index, strand) = gene_annotation[gene_id]
                ga_fp.write(f'{gene_id},{sample_id},{seq_id},{length},{gene_name},{gene_product},{gene_index},{strand}\n')
    
    gene_position_fn = os.path.join(sample_dir, sample_id +'.gene_position')
    with open(gene_position_fn, 'w') as gp_fp:
        for seq_id in gene_position:                
                gp_fp.write(f'{sample_id},{seq_id}')
                for gene_id in gene_position[seq_id]:
                    gp_fp.write(f',{gene_id}')
                gp_fp.write('\n')
    return annotation_fn, gene_position_fn 


def extract_proteins(samples, out_dir, gene_annotation, gene_position, table, threads):
    starttime = datetime.now()

    if threads == 0:
        threads = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(process_single_sample, out_dir=out_dir, table=table), samples)

    for sample, result in zip(samples, results):
        gene_annotation.update(result[0])
        gene_position[sample['id']] = result[1]

    elapsed = datetime.now() - starttime
    logging.info(f'Extract protein -- time taken {str(elapsed)}')


def extract_proteins_tofile(samples, out_dir, gene_annotation_fn, gene_position_fn, table, existing_gene_annotation_fn=None, existing_gene_position_fn=None, threads=0):
    """
    Extract annotations of all the samples, and store in the gene annotation and gene position files
    For now run in single thread, and will convert to asynchronous multi-threaded
    """
    starttime = datetime.now()   
    
    if threads == 0:
        threads = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(process_single_sample, out_dir=out_dir, table=table), samples)

    with gzip.open(gene_annotation_fn,'wt') as ga_fp, gzip.open(gene_position_fn,'wt') as gp_fp: 
        # If there are existing files then copy over
        if existing_gene_annotation_fn:
            with gzip.open(existing_gene_annotation_fn, 'rt') as ega_fp:
                for line in ega_fp.readlines():
                    ga_fp.write(line)
        else:
            ga_fp.write('gene_id,sample_id,seq_id,length,gene_name,gene_product,gene_index,strand\n')

        if existing_gene_position_fn:
            with gzip.open(existing_gene_position_fn, 'rt') as egp_fp:
                for line in egp_fp.readlines():
                    gp_fp.write(line)
    
        for result in results:
            s_annotation_fn, s_gene_position_fn = result
            with open(s_annotation_fn) as fn:
                for line in fn.readlines():
                    ga_fp.write(line)
            os.remove(s_annotation_fn)
            
            with open(s_gene_position_fn) as fn:
                for line in fn.readlines():
                    gp_fp.write(line) 
            os.remove(s_gene_position_fn)   
    elapsed = datetime.now() - starttime
    logging.info(f'Extract protein -- time taken {str(elapsed)}')


def combine_proteins(out_dir, samples):
    # starttime = datetime.now()

    combined_faa_file = os.path.join(out_dir, 'temp', 'combined.faa')
    #TODO: gzip this?
    protein_files = os.path.join(out_dir, 'temp', 'protein.txt')
    with open(combined_faa_file, 'w') as fh:
    #with open(protein_files, 'w') as fh:
        for sample in samples:
            sample_id = sample['id']
            faa_file = os.path.join(out_dir, 'samples', sample_id, sample_id + '.faa.gz')            
            if os.path.isfile(faa_file):
                with gzip.open(faa_file, 'rt') as in_fn:
                    seqs = list(SeqIO.parse(in_fn, 'fasta'))
                    SeqIO.write(seqs, fh, 'fasta')
                #fh.write(faa_file + '\n')
            else:
                raise Exception(f'{faa_file} does not exist')
    #cmd = f"cat {protein_files} | xargs cat > {combined_faa_file}"
    #os.system(cmd)

    # elapsed = datetime.now() - starttime
    # logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file
