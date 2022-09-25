#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    The entry point
"""
import argparse
import os
import shutil
import logging
import sys
import gzip
import csv
from datetime import datetime
import multiprocessing

from pan_genome import wrapper
from pan_genome.utils import check_dir_exist, check_create_folder

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)

def collect_sample(args, previous_samples=None):
    """
    Collect sample from command-line input.

    There are three ways to input data: tsv file (-f/--tsv), gff files 
    (-g/--gff) and fasta files (-a/--fasta). Samples is collected from 
    only one of these, followed the mentioned order.

    Parameters
    ----------
    args : object
        Command-line input arguments.
    previous_samples : list of str, default 'None'
        Existing samples from the previous collection. 
        It is None if we run the init pipeline.
    
    Returns
    -------
    list 
        List of sample from the command line. 
        The list is sorted by sample ID. 
        Each sample is a dictionary {'id':, 'name':, gff_file':, 'assembly':} 
    """
    ## TODO Accept both gff and fasta file at the same time.
    samples = []
    if previous_samples == None:
        previous_samples = []
    sample_id = len(previous_samples) # continue previous id
    if args.tsv != None:
        with open(args.tsv,'r') as fh:
            csv_reader = csv.reader(fh, delimiter='\t')
            for row in csv_reader:
                sample_name = row[0]
                if sample_name in previous_samples:
                    logging.info(f'{sample_name} already exists -- skip')
                    continue
                else:
                    previous_samples.append(sample_name)
                gff = row[1]
                check_dir_exist(gff)
                assembly = row[2]
                if assembly == '':
                    assembly = None # then the GFF file must contain assembly data
                else:
                    check_dir_exist(assembly)
                samples.append({'id':str(sample_id), 'name':sample_name, 
                    'gff_file':gff, 'assembly':assembly})
                sample_id += 1
    elif args.gff != None:
        gff_list = args.gff
        for gff in gff_list:
            base_name = os.path.basename(gff)
            if base_name.endswith('.gz'):
                sample_name = base_name.rsplit('.', 2)[0]
            else:
                sample_name = base_name.rsplit('.', 1)[0]
            if sample_name in previous_samples:
                logging.info(f'{sample_name} already exists -- skip')
                continue
            else:
                previous_samples.append(sample_name)
            check_dir_exist(gff)
            samples.append({'id':str(sample_id), 'name':sample_name, 
                'gff_file':gff, 'assembly':None})
            sample_id += 1
    elif args.fasta != None:
        fasta_list = args.fasta
        for fasta in fasta_list:
            base_name = os.path.basename(fasta)
            if base_name.endswith('.gz'):
                sample_name = base_name.rsplit('.', 2)[0]
            else:
                sample_name = base_name.rsplit('.', 1)[0]
            if sample_name in previous_samples:
                logging.info(f'{sample_name} already exists -- skip')
                continue
            else:
                previous_samples.append(sample_name)
            check_dir_exist(fasta)
            samples.append({'id':str(sample_id), 'name':sample_name, 
                'gff_file':None, 'assembly':fasta})
            sample_id += 1   
    else:
        raise Exception(f'There is no input file')
    
    samples.sort(key= lambda x:x['name'])
    return samples


def get_previous_cluster(summary_file):
    """
    Re-create previous cluster.

    Parameters
    ----------
    summary_file : path
        summary file of previous pan-genome.
    
    Returns
    -------
    list
        a list of []
        each empty list correspond to a previous cluster.

    """
    with open(summary_file, 'r') as fh:
        for line in fh:
            cells = line.rstrip().split('\t')
            if cells[0] == "Total genes":
                num_of_clusters = int(cells[2])
    previous_clusters = [] 
    for i in range(0, num_of_clusters):
        previous_clusters.append([])
    return previous_clusters


def get_previous_samples(sample_file):
    """
    Get previous sample ID.

    Parameters
    ----------
    sample_file : path
        file contains samples of previous pan-genome.
    
    Returns
    -------
    list
        a list of previous sample.
    """
    previous_samples = []
    with open(sample_file, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            sample_name = line.split('\t')[0]
            previous_samples.append(sample_name)
    return previous_samples

def write_sample_file(samples, collection_dir):
    """
    Write sample id to a samples.tsv.

    Parameters
    ----------
    samples : list 
        List of sample from the command line. 
        The list is sorted by sample name. 
        Each sample is a dictionary {'id':, 'name':, 'gff_file':, 'assembly':}
    
    collection_dir : path
        directory of the collection.

    """
    sample_file = os.path.join(collection_dir, 'samples.tsv')
    with open(sample_file, 'a') as fh:
        for sample in samples:
            sample_id = sample['id']
            sample_name = sample['name']
            fh.write(sample_name + '\t' + sample_id + '\n')



def init_function(args):
    """
    Parse arguments and call Init pipeline.

    Parameters
    ----------
    args : object
        Command-line input arguments.
    """
    starttime = datetime.now()

    # parse arguments
    resume = [args.resume]

    collection_dir = args.outdir
    if os.path.exists(collection_dir):
        shutil.rmtree(collection_dir)
        os.makedirs(collection_dir)
    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        if resume[0] == False:
            shutil.rmtree(temp_dir)
            os.makedirs(temp_dir)
    else:
        os.makedirs(temp_dir)
    timing_log = os.path.join(collection_dir, 'time.log')
    baseDir = os.path.dirname(os.path.realpath(__file__))
    # collect samples
    samples = collect_sample(args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')
    
    # call Init pipeline
    wrapper.run_init_pipeline(
        samples, collection_dir, temp_dir, baseDir, args, timing_log, resume)

    # write sample file
    sample_file = os.path.join(collection_dir, 'samples.tsv')
    if os.path.exists(sample_file):
        os.remove(sample_file) # remove existing file
    write_sample_file(samples, collection_dir)

    shutil.rmtree(temp_dir)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

    
def add_function(args):
    """
    Parse arguments and call Add pipeline.

    Parameters
    ----------
    args : object
        Command-line input arguments
    """
    starttime = datetime.now()

    # parse arguments
    resume = [args.resume]

    collection_dir = args.outdir
    check_dir_exist(collection_dir)
    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        if resume[0] == False:
            shutil.rmtree(temp_dir)
            os.makedirs(temp_dir)
    else:
        os.makedirs(temp_dir)
    timing_log = os.path.join(collection_dir, 'time.log')
    baseDir = os.path.dirname(os.path.realpath(__file__))
    # Check required files
    old_representative_fasta = os.path.join(
        collection_dir, 'reference_pangenome.fasta')
    clusters_dir = os.path.join(collection_dir, 'clusters')
    old_cluster_info_file = os.path.join(collection_dir, 'cluster_info.csv')
    summary_file = os.path.join(collection_dir, 'summary_statistics.txt')
    sample_file = os.path.join(collection_dir, 'samples.tsv')
    check_dir_exist(old_representative_fasta)
    check_dir_exist(clusters_dir)
    check_dir_exist(old_cluster_info_file)
    check_dir_exist(summary_file)
    previous_clusters = get_previous_cluster(summary_file)
    # collect new samples
    previous_samples = get_previous_samples(sample_file)
    new_samples = collect_sample(args, previous_samples)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    
    # call Add pipeline
    wrapper.run_add_pipeline(
        new_samples, old_representative_fasta, previous_clusters, 
        collection_dir, temp_dir, baseDir, args, timing_log, resume)
    
    write_sample_file(new_samples, collection_dir)

    shutil.rmtree(temp_dir)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

def main():
    """Setup command-line interface"""
    parser = argparse.ArgumentParser()
    pipeline_help = """
        Analysis pipeline: 
        run initial pan-genome analysis (init), 
        add new samples to previous collection (add)
        """
    parser.add_argument(
        '-p', '--pipeline', help = pipeline_help,
        action='store', choices=['init', 'add'], required=True)
    parser.add_argument(
        '-g', '--gff', help='genome annotation input files (e.g. from Prokka)',
        default=None, nargs='*', type=str)
    parser.add_argument(
        '-a', '--fasta', help='genome assembly input files',
        default=None, nargs='*', type=str)
    tsv_help = """
        when GFF file and FASTA file are seperated (e.g. produced by tools 
        other than Prokka), they could be input through a tsv file. It should
        have 3 columns (without header): Sample ID, path to GFF, path to FASTA.  
        """
    parser.add_argument(
        '-f', '--tsv', help=tsv_help, default=None, type=str)
    parser.add_argument(
        '-o', '--outdir', 
        help='output directory or directory of previous collection ', 
        required=True, type=str) 
    parser.add_argument(
        '-d', '--diamond', 
        help='use Diamond for all-agaist-all alignment instead of Blastp', 
        default=False, action='store_true')
    parser.add_argument(
        '-i', '--identity', 
        help='minimum percentage identity (0..100)', default=95, type=float)
    parser.add_argument(
        '-c', '--coverage', 
        help='minimum percentage coverage (0..100)', default=0, type=float)
    parser.add_argument(
        '-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    parser.add_argument(
        '-t', '--threads', help='number of threads to use, 0 for all', 
        default=0, type=int)
    parser.add_argument(
        '--table', help='codon table', default=11, type=int)
    parser.add_argument(
        '-s', '--split', 
        help='Split paralogs from ortholog', 
        default=False, action='store_true')
    parser.add_argument(
        '-r', '--resume', 
        help='Resume the analysis from interuption', 
        default=False, action='store_true')
    parser.add_argument(
        '-as', '--addstrand2gene', 
        help='Add strand direction (+, -) to gene name', 
        default=False, action='store_true')
    # Execute parse_args()
    args = parser.parse_args()
    if args.threads <= 0:
        args.threads = multiprocessing.cpu_count()

    if args.pipeline == 'init':
        init_function(args)
    elif args.pipeline == 'add':
        add_function(args)

if __name__ == "__main__":
    main()

