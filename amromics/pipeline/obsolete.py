# -*- coding: utf-8 -*-
import json
import os
import shutil
import gzip
import logging
import pandas as pd

from Bio import SeqIO
logger = logging.getLogger(__name__)

def run_alignment(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run phylogenetic analysis of gene clusters. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    gene_cluster_file = report['roary'] + '/gene_presence_absence.csv.gz'
    dict_cds = {}
    for sample in report['samples']:
        with gzip.open(os.path.join(sample['annotation'], sample['id'] + '.ffn.gz'), 'rt') as fn:
            for seq in SeqIO.parse(fn, 'fasta'):
                dict_cds[seq.id] = seq

    # make folder contains sequences for each gene
    alignment_dir = os.path.join(collection_dir, 'alignments')
    gene_df = pd.read_csv(gene_cluster_file, dtype=str, compression='gzip')
    gene_df.fillna('', inplace=True)

    sample_columns = list(gene_df.columns)[14:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                # roary can pool together genes from the same sample and tab-separate them
                for sample_gene in row[sample_column].split('\t'):
                    gene_list.append(sample_gene)
                    # TODO: make sure all samples in this gene have not updated

        gene_list = sorted(gene_list)
        # Only analyse if there are at least 3 genes
        if len(gene_list) < 3:
            logger.info('There are too few genes for {} skipping'.format(gene_id))
            continue

        gene_dir = os.path.join(alignment_dir, gene_id)
        # Check if done before
        gene_list_json = os.path.join(gene_dir, 'gene_list.json')
        # if os.path.isfile(os.path.join(gene_dir, 'parsnp.tree')) and (not overwrite):
        if not overwrite:
            if os.path.isfile(gene_list_json):
                with open(gene_list_json) as fn:
                    existing_gene_list = json.load(fn)
                    if gene_list == existing_gene_list:
                        logger.info('Phylogeny for gene {} done, skipping'.format(gene_id))
                        continue  # for _, row

        gene_file_dir = os.path.join(gene_dir, 'files')
        if not os.path.exists(gene_file_dir):
            os.makedirs(gene_file_dir)

        gene_files = []
        for sample_gene in gene_list:
            gene_file = os.path.join(gene_file_dir, sample_gene + '.fasta')
            SeqIO.write(dict_cds[sample_gene], gene_file, 'fasta')
            gene_files.append(gene_file)

        # Use the first gene as the reference
        cmd = 'parsnp -d {} -r {} -o {} -p {}'.format(
            ' '.join(gene_files[1:]), gene_files[0], gene_dir, threads)
        ret = run_command(cmd, timing_log)
        # if ret != 0:
        #     raise Exception('error')

        with open(gene_list_json, 'w') as fn:
            json.dump(gene_list, fn)
        run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.xmfa')))
        run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.ggr')))

        if os.path.exists(gene_file_dir):
            shutil.rmtree(gene_file_dir)
        #clean up
        run_command('rm -f ' + os.path.join(gene_dir, '*.ini ') + os.path.join(gene_dir, '*block* '))
        shutil.rmtree(os.path.join(gene_dir, 'blocks'), True)
        shutil.rmtree(os.path.join(gene_dir, 'tmp'), True)

    report['alignments'] = alignment_dir
    return report

def run_phylogeny(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run parsnp to create phylogeny tree. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    genome_dir = os.path.join(collection_dir, 'temp/fasta')
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)

    report['phylogeny'] = phylogeny_folder
    phylogeny_file = os.path.join(phylogeny_folder, 'parsnp.tree')
    if os.path.isfile(phylogeny_file) and (not overwrite):
        logger.info('phylogeny tree exists and input has not changed, skip phylogeny analysis')
        return report

    temp_folder = os.path.join(collection_dir, 'temp_phylo')
    if not os.path.isdir(temp_folder):
        os.makedirs(temp_folder)
    reference_genome = None
    sample_list = []
    for i, sample in enumerate(report['samples']):
        fasta_file = os.path.join(temp_folder, sample['id'] + '.fasta')
        cmd = 'gunzip -c {} > {}'.format(sample['assembly'], fasta_file)
        run_command(cmd)
        if i == 0:
            reference_genome = fasta_file
        else:
            sample_list.append(fasta_file)
    cmd = 'parsnp -r {} -d {} -o {} -p {}'.format(
        reference_genome,
        ' '.join(sample_list),
        phylogeny_folder, threads)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parsnp')
    run_command('gzip {}'.format(os.path.join(phylogeny_folder, 'parsnp.xmfa')))
    run_command('gzip {}'.format(os.path.join(phylogeny_folder, 'parsnp.ggr')))
    shutil.rmtree(temp_folder)
    return report