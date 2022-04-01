import argparse
import os
import shutil
import logging
import sys
import gzip
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

def main_function(args):
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
    wrapper.run_main_pipeline(samples, collection_dir, temp_dir, baseDir, args, timing_log)

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
    old_represent_faa = os.path.join(collection_dir, 'representative.fasta')
    if not os.path.isfile(old_represent_faa):
        raise Exception(f'{old_represent_faa} does not exist')
    old_clusters, old_clusters_annotation, gene_to_old_cluster = new_pipeline.read_database(old_represent_faa)

    old_spreadsheet_file = os.path.join(collection_dir, 'gene_presence_absence.csv.gz')
    if not os.path.isfile(old_spreadsheet_file):
        raise Exception(f'{old_spreadsheet_file} does not exist')
    with gzip.open(old_spreadsheet_file, 'rt') as fh:
        csv.field_size_limit(sys.maxsize)
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        sample_id_list = header[8:]
    
    old_rtab_file = os.path.join(collection_dir, 'gene_presence_absence.Rtab.gz')
    if not os.path.isfile(old_rtab_file):
        raise Exception(f'{old_rtab_file} does not exist')


    # collect new samples
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    
    # pipeline

    new_clusters, gene_to_new_cluster, unmatched_represent_faa, gene_dictionary = wrapper.add_sample(new_samples, old_represent_faa, old_clusters, gene_to_old_cluster, collection_dir, temp_dir, args, timing_log)

    if args.fasta == None:
        new_clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=new_clusters, 
            gene_dictionary=gene_dictionary,
            start = 1)
    else:
        # create representative
        new_represent_fasta = os.path.join(temp_dir, 'representative.fasta')
        utils.create_fasta_include(
            fasta_file_list=[unmatched_represent_faa], 
            include_list=gene_to_new_cluster, 
            output_file=new_represent_fasta
            )
        
        new_clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=new_clusters,
            rep_fasta = new_represent_fasta,
            temp_dir=temp_dir,
            baseDir = baseDir,
            timing_log = timing_log,
            threads = args.threads,
            start = 1
            )

    output.create_representative_fasta(
        gene_to_cluster=gene_to_new_cluster, 
        clusters_annotation=new_clusters_annotation, 
        source_fasta=unmatched_represent_faa, 
        out_fasta=os.path.join(collection_dir, 'representative.fasta'), 
        mode='a')

    output.update_spreadsheet(old_spreadsheet_file, old_clusters, new_clusters, new_clusters_annotation, gene_dictionary, new_samples, temp_dir)
    rtab_file = output.update_rtab(old_rtab_file, old_clusters, new_clusters, new_clusters_annotation, gene_dictionary, new_samples, temp_dir)
    output.create_summary(rtab_file, collection_dir)

    # shutil.rmtree(temp_dir)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

def new_function(args):
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

    # collect new samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    # pipeline
    species_db = args.db
    if not os.path.isfile(species_db):
        raise Exception(f'{species_db} does not exist')
    old_clusters, old_clusters_annotation, gene_to_cluster = new_pipeline.read_database(species_db)

    new_clusters, gene_to_new_cluster, unmatched_represent_faa, gene_dictionary = wrapper.add_sample(samples, species_db, old_clusters, gene_to_cluster, collection_dir, temp_dir, args, timing_log)

    if args.fasta == None:
        clusters = new_pipeline.combine_clusters(old_clusters, new_clusters)
        clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=clusters, 
            gene_dictionary=gene_dictionary)
    else:
        # create new representative
        new_represent_fasta = os.path.join(temp_dir, 'representative.fasta')
        utils.create_fasta_include(
            fasta_file_list=[unmatched_represent_faa], 
            include_list=gene_to_new_cluster, 
            output_file=new_represent_fasta
            )  
        
        new_clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=new_clusters,
            rep_fasta = new_represent_fasta,
            temp_dir=temp_dir,
            baseDir = baseDir,
            timing_log=timing_log,
            threads = args.threads,
            start=len(old_clusters) + 1
            )

        clusters, clusters_annotation=new_pipeline.combine_result(
            old_clusters, old_clusters_annotation, new_clusters, new_clusters_annotation)

    # output
    output.create_spreadsheet(clusters, clusters_annotation, gene_dictionary, samples, collection_dir)
    rtab_file = output.create_rtab(clusters, clusters_annotation, gene_dictionary,samples,collection_dir)
    output.create_summary(rtab_file, collection_dir)

    # shutil.rmtree(temp_dir)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    
    main_cmd = subparsers.add_parser(
        'main',
        description='Main pipeline: run pan-genome analysis for the first time',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    main_cmd.set_defaults(func=main_function)
    main_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    main_cmd.add_argument('-b', '--fasta', help='assembly input files',default=None, nargs='*', type=str)
    main_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    main_cmd.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    main_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    main_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.95, type=float)
    main_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0, type=float)
    main_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    main_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0.95, type=float)
    main_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    main_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    main_cmd.add_argument('--table', help='codon table', default=11, type=int)

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
    add_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0.95, type=float)
    add_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    add_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--table', help='codon table', default=11, type=int)

    new_cmd = subparsers.add_parser(
        'new',
        description='new pipeline: run pan-genome analysis using species database',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    new_cmd.set_defaults(func=new_function)
    new_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    new_cmd.add_argument('-b', '--fasta', help='assembly input files',default=None, nargs='*', type=str)
    new_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    new_cmd.add_argument('--db', help='database directory', required=True, type=str)
    new_cmd.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    new_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    new_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.95, type=float)
    new_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0, type=float)
    new_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    new_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0.95, type=float)
    new_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    new_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    new_cmd.add_argument('--table', help='codon table', default=11, type=int)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

