import argparse
import os
import shutil
import logging
import sys
import gzip
import json
import csv
from datetime import datetime
from pan_genome import *
from pan_genome import wrapper

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)

def collect_sample(sample_id_list, args):
    samples = []
    if args.tsv != None:
        with open(args.tsv,'r') as fh:
            csv_reader = csv.reader(fh, delimiter='\t')
            for row in csv_reader:
                gff = row[1]
                sample_id = row[0]
                if sample_id in sample_id_list:
                    logging.info(f'{sample_id} already exists -- skip')
                    continue
                else:
                    sample_id_list.append(sample_id)
                assembly = row[2]
                if row[2] == '':
                    assembly = None
                samples.append({'id':sample_id, 'gff_file':gff, 'assembly':assembly}) 
    elif args.gff != None:
        gff_list = args.gff
        for gff in gff_list:
            base_name = os.path.basename(gff)
            if base_name.endswith('.gz'):
                sample_id = base_name.rsplit('.', 2)[0]
            else:
                sample_id = base_name.rsplit('.', 1)[0]
            if sample_id in sample_id_list:
                logging.info(f'{sample_id} already exists -- skip')
                continue
            else:
                sample_id_list.append(sample_id)
            samples.append({'id':sample_id, 'gff_file':gff, 'assembly':None})
    elif args.fasta != None:
            fasta_list = args.fasta
            for fasta in fasta_list:
                base_name = os.path.basename(fasta)
                if base_name.endswith('.gz'):
                    sample_id = base_name.rsplit('.', 2)[0]
                else:
                    sample_id = base_name.rsplit('.', 1)[0]
                if sample_id in sample_id_list:
                    logging.info(f'{sample_id} already exists -- skip')
                    continue
                else:
                    sample_id_list.append(sample_id)
                samples.append({'id':sample_id, 'gff_file':None, 'assembly':fasta})    
    else:
        raise Exception(f'There is no input file')
    
    samples.sort(key= lambda x:x['id'])
    return samples

def init_function(args):
    starttime = datetime.now()

    # parse arguments
    collection_dir = args.outdir
    if not os.path.exists(collection_dir):
        os.makedirs(collection_dir)
    
    timing_log = os.path.join(collection_dir, 'time.log')

    temp_dir = os.path.join(collection_dir, 'temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)   
    
    baseDir = os.path.dirname(os.path.realpath(__file__))

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')
    
    # pipeline
    clusters, gene_dictionary = wrapper.run_main_pipeline(samples, collection_dir, temp_dir, baseDir, args, timing_log)

    if args.fasta == None:
        clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=clusters, 
            gene_dictionary=gene_dictionary)
        output.create_output(clusters, clusters_annotation, gene_dictionary, samples, collection_dir)
        representative_fasta = alignment.main_create_msa(clusters, samples, collection_dir, baseDir, args.threads)
    
    else:        
        representative_fasta = alignment.main_create_msa(clusters, samples, collection_dir, baseDir, args.threads)
        clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=clusters,
            rep_fasta = representative_fasta,
            temp_dir=temp_dir,
            baseDir = baseDir,
            timing_log=timing_log,
            threads = args.threads)
        output.create_output(clusters, clusters_annotation, gene_dictionary, samples, collection_dir)
    # shutil.rmtree(temp_dir)
        
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

    

def add_function(args):
    starttime = datetime.now()

    # parse arguments
    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    
    timing_log = os.path.join(collection_dir, 'time.log')
    
    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
    else:
        os.makedirs(temp_dir)    
    
    baseDir = os.path.dirname(os.path.realpath(__file__))

    # Check required files
    old_representative_fasta = os.path.join(collection_dir, 'reference_pangenome.fasta')
    if not os.path.isfile(old_representative_fasta):
        raise Exception(f'{old_representative_fasta} does not exist')
    
    clusters_dir = os.path.join(collection_dir, 'clusters')
    if not os.path.exists(clusters_dir):
        raise Exception(f'{clusters_dir} does not exist')

    samples_dir = os.path.join(collection_dir, 'samples')
    if not os.path.exists(samples_dir):
        raise Exception(f'{samples_dir} does not exist')

    old_presence_absence_file = os.path.join(collection_dir, 'gene_presence_absence.csv.gz')
    if not os.path.isfile(old_presence_absence_file):
        raise Exception(f'{old_presence_absence_file} does not exist')

    old_cluster_info_file = os.path.join(collection_dir, 'cluster_info.csv')
    if not os.path.isfile(old_cluster_info_file):
        raise Exception(f'{old_cluster_info_file} does not exist')
    
    summary_file = os.path.join(collection_dir, 'summary_statistics.txt')
    with open(summary_file, 'r') as fh:
        for line in fh:
            cells = line.rstrip().split('\t')
            if cells[0] == "Total genes":
                num_of_clusters = int(cells[2])
                old_clusters = [] # create a list of empty lists inside, used to grab new seqs
                for i in range(0, num_of_clusters):
                    old_clusters.append([])
    
    with gzip.open(old_presence_absence_file, 'rt') as fh:
        csv.field_size_limit(sys.maxsize)
        reader = csv.reader(fh, delimiter=',')
        header = next(reader)
        sample_id_list = header[1:] # exclude ID column
    
    # collect new samples
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    
    # pipeline
    new_clusters, gene_dictionary = wrapper.add_sample(new_samples, old_representative_fasta, old_clusters, collection_dir, temp_dir, args, timing_log)
    
    if args.fasta == None:
        new_clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=new_clusters, 
            gene_dictionary=gene_dictionary)
        output.update_output(old_clusters, new_clusters, new_clusters_annotation, gene_dictionary, new_samples, temp_dir, collection_dir)
        new_representative_fasta = alignment.add_create_msa(old_clusters, new_clusters, new_samples, collection_dir, baseDir, args.threads)
    else:
        new_representative_fasta = alignment.add_create_msa(old_clusters, new_clusters, new_samples, collection_dir, baseDir, args.threads)
        new_clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=new_clusters,
            rep_fasta = new_representative_fasta,
            temp_dir=temp_dir,
            baseDir = baseDir,
            timing_log = timing_log,
            threads = args.threads)
        output.update_output(old_clusters, new_clusters, new_clusters_annotation, gene_dictionary, new_samples, temp_dir, collection_dir)

    # shutil.rmtree(temp_dir)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    
    init_cmd = subparsers.add_parser(
        'init',
        description='Init pipeline: run initial pan-genome analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    init_cmd.set_defaults(func=init_function)
    init_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    init_cmd.add_argument('-b', '--fasta', help='assembly input files',default=None, nargs='*', type=str)
    init_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    init_cmd.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    init_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    init_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.95, type=float)
    init_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0, type=float)
    init_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    init_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0, type=float)
    init_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    init_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    init_cmd.add_argument('--table', help='codon table', default=11, type=int)

    add_cmd = subparsers.add_parser(
        'add',
        description='Add pipeline: add sample into previous collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=add_function)
    add_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    add_cmd.add_argument('-b', '--fasta', help='assembly input files',default=None, nargs='*', type=str)
    add_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    add_cmd.add_argument('-o', '--collection-dir', help='previous collection directory', required=True, type=str)
    add_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    add_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.95, type=float)
    add_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0, type=float)
    add_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    add_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0, type=float)
    add_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    add_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--table', help='codon table', default=11, type=int)

if __name__ == "__main__":
    main()

