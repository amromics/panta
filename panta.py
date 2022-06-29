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
from pan_genome import annotate
from pan_genome import alignment
from pan_genome import output
from pan_genome.utils import check_dir_exist, check_create_folder

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)

def collect_sample(args, previous_sample_id=None):
    """
    Collect sample from command-line input.

    There are three ways to input data: tsv file (-f/--tsv), gff files 
    (-g/--gff) and fasta files (-b/--fasta). Samples is extracted by 
    only one of these ways, followed the mentioned order.

    Parameters
    ----------
    args : object
        Command-line input arguments.
    previous_sample_id : list of str, default 'None'
        Existing sample ID from the previous collection. 
        It is None if we run the init pipeline.
    Returns
    -------
    list 
        List of sample from the command line. 
        The list is sorted by sample ID. 
        Each sample is a dictionary {'id':, 'gff_file':, 'assembly':} 
    """
    ## TODO Accept both gff and fasta file at the same time.
    samples = []
    if previous_sample_id == None:
        previous_sample_id = []
    if args.tsv != None:
        with open(args.tsv,'r') as fh:
            csv_reader = csv.reader(fh, delimiter='\t')
            for row in csv_reader:
                gff = row[1]
                sample_id = row[0]
                if sample_id in previous_sample_id:
                    logging.info(f'{sample_id} already exists -- skip')
                    continue
                else:
                    previous_sample_id.append(sample_id)
                assembly = row[2]
                if row[2] == '':
                    assembly = None
                samples.append(
                    {'id':sample_id, 'gff_file':gff, 'assembly':assembly}) 
    elif args.gff != None:
        gff_list = args.gff
        for gff in gff_list:
            base_name = os.path.basename(gff)
            if base_name.endswith('.gz'):
                sample_id = base_name.rsplit('.', 2)[0]
            else:
                sample_id = base_name.rsplit('.', 1)[0]
            if sample_id in previous_sample_id:
                logging.info(f'{sample_id} already exists -- skip')
                continue
            else:
                previous_sample_id.append(sample_id)
            samples.append({'id':sample_id, 'gff_file':gff, 'assembly':None})
    elif args.fasta != None:
            fasta_list = args.fasta
            for fasta in fasta_list:
                base_name = os.path.basename(fasta)
                if base_name.endswith('.gz'):
                    sample_id = base_name.rsplit('.', 2)[0]
                else:
                    sample_id = base_name.rsplit('.', 1)[0]
                if sample_id in previous_sample_id:
                    logging.info(f'{sample_id} already exists -- skip')
                    continue
                else:
                    previous_sample_id.append(sample_id)
                samples.append(
                    {'id':sample_id, 'gff_file':None, 'assembly':fasta})    
    else:
        raise Exception(f'There is no input file')
    
    samples.sort(key= lambda x:x['id'])
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


def get_previous_sample_id(presence_absence_file):
    """
    Get previous sample ID.

    Parameters
    ----------
    presence_absence_file : path
        Gene presence absence file of previous pan-genome.
    
    Returns
    -------
    list
        a list of previous sample ID.
    """
    with gzip.open(presence_absence_file, 'rt') as fh:
        csv.field_size_limit(sys.maxsize)
        reader = csv.reader(fh, delimiter=',')
        header = next(reader)
        previous_sample_id = header[1:] # exclude ID column
        return previous_sample_id


def init_function(args):
    """
    Parse arguments and run Init pipeline.

    Parameters
    ----------
    args : object
        Command-line input arguments.
    """
    starttime = datetime.now()

    # parse arguments
    collection_dir = args.outdir
    temp_dir = os.path.join(collection_dir, 'temp')
    check_create_folder(collection_dir)
    check_create_folder(temp_dir)
    timing_log = os.path.join(collection_dir, 'time.log')
    baseDir = os.path.dirname(os.path.realpath(__file__))
    # collect samples
    samples = collect_sample(args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')
    
    # call clustering pipeline
    clusters, gene_dictionary = wrapper.run_init_pipeline(
        samples, collection_dir, temp_dir, args, timing_log)

    # annotate clusters, create gene alignment and output
    if args.fasta == None:
        clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=clusters, 
            gene_dictionary=gene_dictionary)
        output.create_output(
            clusters, clusters_annotation, 
            gene_dictionary, samples, collection_dir)
        representative_fasta = alignment.main_create_msa(
            clusters, samples, collection_dir, baseDir, args.threads)
    else:        
        representative_fasta = alignment.main_create_msa(
            clusters, samples, collection_dir, baseDir, args.threads)
        clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=clusters,
            rep_fasta = representative_fasta,
            temp_dir=temp_dir,
            baseDir = baseDir,
            timing_log=timing_log,
            threads = args.threads)
        output.create_output(
            clusters, clusters_annotation, 
            gene_dictionary, samples, collection_dir)
    
    # shutil.rmtree(temp_dir)    
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

    
def add_function(args):
    """
    Parse arguments and run Add pipeline.

    Parameters
    ----------
    args : object
        Command-line input arguments
    """
    starttime = datetime.now()

    # parse arguments
    collection_dir = args.outdir
    temp_dir = os.path.join(collection_dir, 'temp')
    check_dir_exist(collection_dir)
    check_create_folder(temp_dir)
    timing_log = os.path.join(collection_dir, 'time.log')
    baseDir = os.path.dirname(os.path.realpath(__file__))
    # Check required files
    old_representative_fasta = os.path.join(
        collection_dir, 'reference_pangenome.fasta')
    clusters_dir = os.path.join(collection_dir, 'clusters')
    samples_dir = os.path.join(collection_dir, 'samples')
    old_presence_absence_file = os.path.join(
        collection_dir, 'gene_presence_absence.csv.gz')
    old_cluster_info_file = os.path.join(collection_dir, 'cluster_info.csv')
    summary_file = os.path.join(collection_dir, 'summary_statistics.txt')
    check_dir_exist(old_representative_fasta)
    check_dir_exist(clusters_dir)
    check_dir_exist(samples_dir)
    check_dir_exist(old_presence_absence_file)
    check_dir_exist(old_cluster_info_file)
    check_dir_exist(summary_file)
    previous_clusters = get_previous_cluster(summary_file)
    # collect new samples
    previous_sample_id = get_previous_sample_id(old_presence_absence_file)
    new_samples = collect_sample(args, previous_sample_id)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    
    # call clustering pipeline
    new_clusters, gene_dictionary = wrapper.run_add_pipeline(
        new_samples, old_representative_fasta, previous_clusters, 
        collection_dir, temp_dir, args, timing_log)
    
    # annotate clusters, create gene alignment and output
    if args.fasta == None:
        new_clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=new_clusters, 
            gene_dictionary=gene_dictionary)
        output.update_output(previous_clusters, new_clusters, 
            new_clusters_annotation, gene_dictionary, 
            new_samples, temp_dir, collection_dir)
        new_representative_fasta = alignment.add_create_msa(
            previous_clusters, new_clusters, new_samples, 
            collection_dir, baseDir, args.threads)
    else:
        new_representative_fasta = alignment.add_create_msa(
            previous_clusters, new_clusters, new_samples, 
            collection_dir, baseDir, args.threads)
        new_clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=new_clusters,
            rep_fasta = new_representative_fasta,
            temp_dir=temp_dir,
            baseDir = baseDir,
            timing_log = timing_log,
            threads = args.threads)
        output.update_output(
            previous_clusters, new_clusters, new_clusters_annotation, 
            gene_dictionary, new_samples, temp_dir, collection_dir)

    # shutil.rmtree(temp_dir)
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

def main():
    """Setup command-line interface"""
    parser = argparse.ArgumentParser()
    pipeline_help = """
        Select the pipeline: 
        run initial pan-genome analysis (init), 
        add new samples to previous collection (add)
        """
    parser.add_argument(
        '-p', '--pipeline', help = pipeline_help,
        action='store', choices=['init', 'add'], required=True)
    parser.add_argument(
        '-g', '--gff', help='gff input files',
        default=None, nargs='*', type=str)
    parser.add_argument(
        '-b', '--fasta', help='assembly input files',
        default=None, nargs='*', type=str)
    parser.add_argument(
        '-f', '--tsv', help='tsv input file', default=None, type=str)
    parser.add_argument(
        '-o', '--outdir', 
        help='output directory/previous collection directory', 
        required=True, type=str) 
    parser.add_argument(
        '-d', '--diamond', 
        help='use Diamond for all-agaist-all alignment instead of Blastp', 
        default=False, action='store_true')
    parser.add_argument(
        '-i', '--identity', 
        help='minimum percentage identity', default=0.95, type=float)
    parser.add_argument(
        '--LD', help='length difference cutoff between two sequences', 
        default=0, type=float)
    parser.add_argument(
        '--AL', help='alignment coverage for the longer sequence', 
        default=0, type=float)
    parser.add_argument(
        '--AS', help='alignment coverage for the shorter sequence', 
        default=0, type=float)
    parser.add_argument(
        '-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    parser.add_argument(
        '-t', '--threads', help='number of threads to use, 0 for all', 
        default=0, type=int)
    parser.add_argument(
        '--table', help='codon table', default=11, type=int)

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

