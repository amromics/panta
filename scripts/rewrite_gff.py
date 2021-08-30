import argparse
import os
import re
import shutil
import csv
import logging
import json
from glob import glob

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)

def check_inputs(args):

    gff_samples = []
    for gff in args.gff_files:
        if not gff.endswith('gff'):
            raise Exception(f'{gff} should be a gff3 file')
        base_name = os.path.basename(gff)
        sample_id = base_name.rsplit('.', 1)[0]
        gff_samples.append(sample_id)
    
    for row in csv.reader(open(args.roary, 'r')):
        if row[0] != 'Gene':
            break
        result_samples = row[14:]
    
    if len(gff_samples) != len(result_samples):
        raise Exception(f'There are {len(gff_samples)} gff files while there are {len(result_samples)} samples from roary')

    logging.info(f'Check input')


def parse_roary_result(roary_file):
    cluster_id = 1
    gene_dict = {}
    for row in csv.reader(open(roary_file, 'r')):
        if row[0] == 'Gene':
            samples = row[14:]
            continue
        for cell, sample in zip(row[14:], samples):
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                gene_dict.setdefault(sample, {})[gene] = cluster_id
        cluster_id += 1

    logging.info(f'Parse roary result')
    return gene_dict


def rewrite_gff(infile, outfile, dictionary):
    with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
        for line in in_fh:
            cells = line.rstrip().split('\t')
            if len(cells) != 9:
                out_fh.write(line)
                continue
            if cells[2] != 'CDS':
                out_fh.write(line)
                continue

            tags = cells[8].split(';')
            found = False
            for tag in tags:
                ID = re.match(r"^ID=(.+)", tag)
                if ID != None:
                    gene_id = ID.group(1)
                    if gene_id in dictionary:
                        cluster_id = dictionary[gene_id]
                        found = True
                    break
            if found == True:
                tags.append('pangenome_id={}'.format(str(cluster_id)))
            
            cells[8] = ';'.join(tags)
            out_fh.write('\t'.join(cells) + '\n')


def rewrite_gffs(gff_files, gene_dict, outdir):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    for gff in gff_files:
        base_name = os.path.basename(gff)
        sample_id = base_name.rsplit('.', 1)[0]
        outfile = os.path.join(outdir, sample_id + '.gff')
        rewrite_gff(gff, outfile, gene_dict[sample_id])

    logging.info(f'Done')


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('gff_files', help='a.gff b.gff ... (*.gff)', type=str, nargs='+')
    parser.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    parser.add_argument('-r', '--roary', help='gene_presence_absence.csv from Roary', required=True, type=str)

    args = parser.parse_args()

    check_inputs(args)
    gene_dict = parse_roary_result(args.roary)
    rewrite_gffs(gff_files=args.gff_files, gene_dict=gene_dict, outdir=args.outdir)


if __name__ == "__main__":
    main()

    # python3 rewrite_gff.py -o a -r /home/ntanh1999/disk/results/roary/Sa110_split/gene_presence_absence.csv /home/ntanh1999/disk/datasets/data/Sa110/*.gff
