import os
import re
import logging
import multiprocessing
from functools import partial
from datetime import datetime
import gzip
from Bio import SeqIO
from Bio.Seq import Seq

from collections import OrderedDict

logger = logging.getLogger(__name__)

def parse_fasta(fa_fh):
    name, desc, seq = None, '', []
    for line in fa_fh:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, desc,''.join(seq))
            toks = line[1:].partition(' ')
            name, desc = toks[0], toks[2]
            seq = []
        else:
            seq.append(line)
    if name:
        yield (name, desc, ''.join(seq))

def parse_gff(gff_fh, sample_id, min_protein_len=40):
    gene_annotation = OrderedDict()
    gene_position = OrderedDict()
    suffix = 1
    bed_records = []
    gene_index = 0
    seq_id = None
    min_cds_len = 3 * min_protein_len

    for line in gff_fh:
        if line.startswith('##FASTA'):
            #Done reading gff, move on to reading fasta
            break

        if line[0] == '#':
            continue
        line = line.strip()
        #print(line)
        cells = line.split('\t')
        if cells[2] != 'CDS':
            continue

        start = int(cells[3])
        end = int(cells[4])
        length = end - start + 1
        if length < min_cds_len:
            continue
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
            if tag.startswith('ID='):
                gene_id = tag[3:]
            elif tag.startswith('gene='):
                gene_name = tag[5:]
                gene_name = re.sub(r'\W', '_', gene_name)
            elif tag.startswith('product='):
                gene_product = tag[8:]
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

    return gene_annotation, gene_position, bed_records


def process_single_sample(sample, out_dir, table):
    sample_id = sample['id']
    sample_dir = os.path.join(out_dir, 'samples', sample_id)

    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    seqs = {}

    sample_assembly = sample['assembly']
    if sample_assembly:
        if sample_assembly.endswith('.gz'):
            fa_fh = gzip.open(sample_assembly,'rt')
        else:
            fa_fh = open(sample_assembly)
        for (seq_id, seq_desc, seq) in parse_fasta(fa_fh):
            seqs[seq_id.replace('-','_')] = seq
        fa_fh.close()

    gff_file = sample['gff_file']
    if gff_file.endswith('gff.gz'):
        in_fh = gzip.open(gff_file,'rt')
    else: #end with .gff
        in_fh = open(gff_file)

    gene_annotation, gene_position, bed_records = parse_gff(in_fh, sample_id)
    if len(seqs) == 0:
        for (seq_id, seq_desc, seq) in parse_fasta(in_fh):
            seqs[seq_id.replace('-','_')] = seq
    in_fh.close()

    fna_file = os.path.join(sample_dir, sample_id +'.fna')
    faa_file = os.path.join(sample_dir, sample_id +'.faa')

    gene_seqs = []
    protein_seqs = []
    for bed_record in bed_records:
        (seq_id, start, end, gene_id, strand) = bed_record
        if seq_id not in seqs:
            continue

        seq = seqs[seq_id]
        if end > len(seq):
            continue

        subseq = Seq(seq[start:end])
        if strand == '-':
            subseq = subseq.reverse_complement()
        gene_seq = SeqIO.SeqRecord(subseq, gene_id)
        gene_seqs.append(gene_seq)
        protein_seq = gene_seq.translate(table=table, stop_symbol='')
        protein_seq.id = gene_id
        protein_seqs.append(protein_seq)
    with open(fna_file, 'w') as fna_fh:
        SeqIO.write(gene_seqs, fna_fh, 'fasta')
    with open(faa_file, 'w') as faa_fh:
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

    elapsed = datetime.now() - starttime
    logging.info(f'Extract protein -- time taken {str(elapsed)}')

    with open(gene_annotation_fn,'w') as ga_fp, open(gene_position_fn,'w') as gp_fp:
        # If there are existing files then copy over
        if existing_gene_annotation_fn:
            with open(existing_gene_annotation_fn) as ega_fp:
                for line in ega_fp.readlines():
                    ga_fp.write(line)
        else:
            ga_fp.write('gene_id,sample_id,seq_id,length,gene_name,gene_product,gene_index,strand\n')

        if existing_gene_position_fn:
            with open(existing_gene_position_fn) as egp_fp:
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
    with open(combined_faa_file, 'w') as fh:
        for sample in samples:
            sample_id = sample['id']
            faa_file = os.path.join(out_dir, 'samples', sample_id, sample_id + '.faa')
            if os.path.isfile(faa_file):
                with open(faa_file) as in_fn:
                    seqs = list(SeqIO.parse(in_fn, 'fasta'))
                    SeqIO.write(seqs, fh, 'fasta')
            else:
                raise Exception(f'{faa_file} does not exist')
    # logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file

def combine_proteins_with_maps(out_dir, samples):
    starttime = datetime.now()
    combined_faa_file = os.path.join(out_dir, 'temp', 'combined.faa')
    combined_faa_map = os.path.join(out_dir, 'temp', 'combined.map')
    count = 0
    with open(combined_faa_file, 'w') as fh, open(combined_faa_map, 'w') as map_fh:
        for sample in samples:
            sample_id = sample['id']
            faa_file = os.path.join(out_dir, 'samples', sample_id, sample_id + '.faa')
            if os.path.isfile(faa_file):
                with open(faa_file) as in_fn:
                    for seq in SeqIO.parse(in_fn, 'fasta'):
                        fh.write(f'>{count}\n{seq.seq}\n')
                        map_fh.write(f'{seq.id}\n')
                        count += 1
            else:
                raise Exception(f'{faa_file} does not exist')

    elapsed = datetime.now() - starttime
    logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file, combined_faa_map
