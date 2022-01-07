import argparse
import os
import shutil
import logging
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
                if not gff.endswith('gff'):
                    raise Exception(f'{gff} should be a gff3 file')
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
            if not gff.endswith('gff'):
                raise Exception(f'{gff} should be a gff3 file')
            base_name = os.path.basename(gff)
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

    collection_dir = args.outdir
    threads = args.threads
    if args.fasta != None:
        annotate = True
    else:
        annotate = False
    
    temp_dir = os.path.join(collection_dir, 'temp')
    if not os.path.exists(collection_dir):
        os.makedirs(collection_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)    

    dir_path = os.path.dirname(os.path.realpath(__file__))
    db_dir = os.path.join(dir_path, 'db')


    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')
    
    gene_annotation = {}
    gene_position = {}
    
    # pipeline
    annotated_clusters = wrapper.run_main_pipeline(samples, gene_annotation, gene_position, collection_dir, temp_dir, db_dir, args, annotate, threads)
    wrapper.create_outputs(gene_annotation,annotated_clusters,samples,collection_dir)
    wrapper.run_gene_alignment(annotated_clusters, gene_annotation, samples, collection_dir, args.alignment, threads)

    # shutil.rmtree(temp_dir)
        
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

    

def add_function(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if args.fasta != None:
        annotate = True
    else:
        annotate = False

    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
    else:
        os.makedirs(temp_dir)

    dir_path = os.path.dirname(os.path.realpath(__file__))
    db_dir = os.path.join(dir_path, 'db')

    # Check required files
    if args.alignment != None:
        samples_dir = os.path.join(collection_dir, 'samples')
        if not os.path.exists(samples_dir):
            raise Exception(f'{samples_dir} does not exist')
    
    gene_annotation_file = os.path.join(collection_dir, 'gene_annotation.tsv')
    if not os.path.isfile(gene_annotation_file):
        raise Exception(f'{gene_annotation_file} does not exist')
    gene_annotation = output.import_gene_annotation(gene_annotation_file)

    gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))
    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))

    old_represent_faa = os.path.join(collection_dir, 'representative.fasta')
    if not os.path.isfile(old_represent_faa):
        raise Exception(f'{old_represent_faa} does not exist')

    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    
    # pipeline
    annotated_clusters = wrapper.run_add_pipeline(old_samples, new_samples, old_represent_faa, old_clusters, gene_annotation, gene_position, temp_dir, collection_dir, db_dir, annotate, threads, args)
    wrapper.create_outputs(gene_annotation,annotated_clusters,new_samples,collection_dir)
    wrapper.run_gene_alignment(annotated_clusters, gene_annotation, new_samples, collection_dir, args.alignment, threads)

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
    main_cmd.add_argument('-s', '--split', help='split paralog clusters', default=False, action='store_true')
    main_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    main_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.95, type=float)
    main_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0, type=float)
    main_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    main_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0, type=float)
    main_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    main_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    main_cmd.add_argument('--table', help='codon table', default=11, type=int)
    main_cmd.add_argument('-a', '--alignment', help='run alignment for each gene cluster', default=None, action='store', choices=['protein', 'nucleotide'])
    main_cmd.add_argument('-n', '--number', help='number of samples which are analysed first', default=0, type=int)
    main_cmd.add_argument('--genus', help='Genus',default=None, type=str)
    main_cmd.add_argument('--species', help='Species',default=None, type=str)

    add_cmd = subparsers.add_parser(
        'add',
        description='Add pipeline: add sample into previous collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=add_function)
    add_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    add_cmd.add_argument('-b', '--fasta', help='assembly input files',default=None, nargs='*', type=str)
    add_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    add_cmd.add_argument('-c', '--collection-dir', help='previous collection directory', required=True, type=str)
    add_cmd.add_argument('-s', '--split', help='split paralog clusters', default=False, action='store_true')
    add_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    add_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.95, type=float)
    add_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0, type=float)
    add_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    add_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0, type=float)
    add_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    add_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--table', help='codon table', default=11, type=int)
    add_cmd.add_argument('-a', '--alignment', help='run alignment for each gene cluster', default=None, action='store', choices=['protein', 'nucleotide'])
    add_cmd.add_argument('--genus', help='Genus',default=None, type=str)
    add_cmd.add_argument('--species', help='Species',default=None, type=str)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    # main()

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s : %(message)s',
        datefmt='%I:%M:%S')
    logger = logging.getLogger(__name__)


    # setup_db(db_dir="/home/ted/amromics/amromics/pan-genome/db", force=True)

    annotate.annotate_cluster(
        rep_fasta='/home/ted/test_prodigal/out/1/representative.fasta', 
        temp_dir='/home/ted/test_prodigal/out/1/temp', 
        threads=4, 
        genus="Staphylococcus"
    )