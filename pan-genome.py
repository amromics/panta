import argparse
import os
import subprocess
import shutil
import logging
import json
import csv
from datetime import datetime
from pan_genome import *

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

def add_samples(temp_dir, new_samples, old_represent_faa, old_clusters, gene_annotation, collection_dir, threads, args):
    new_combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=temp_dir,
        samples=new_samples)

    unmatched_faa, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1 = old_represent_faa,
        database_2 = new_combined_faa,
        out_dir = temp_dir,
        threads=threads)

    gene_to_cluster, old_clusters = add_sample_pipeline.add_gene_cd_hit_2d(old_clusters, cd_hit_2d_clusters)

    num_seq = subprocess.run(f'grep ">" {unmatched_faa} | wc -l', shell=True, capture_output=True, text=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return old_clusters, new_combined_faa

    unmatched_represent_faa, unmatched_clusters = main_pipeline.run_cd_hit(
        faa_file=unmatched_faa,
        out_dir=temp_dir,
        threads=threads)

    blast_1_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = old_represent_faa,
        query_fasta = unmatched_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast1'),
        evalue = args.evalue,
        threads=threads)

    remain_fasta, old_clusters = add_sample_pipeline.add_gene_blast(
        old_clusters=old_clusters,
        gene_to_cluster=gene_to_cluster,
        unmatched_clusters = unmatched_clusters,
        blast_result=blast_1_result, 
        fasta_file=unmatched_represent_faa, 
        out_dir=temp_dir,
        gene_annotation=gene_annotation, 
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    num_seq = subprocess.run(f'grep ">" {remain_fasta} | wc -l', shell=True, capture_output=True, text=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return old_clusters, new_combined_faa
    
    blast_2_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = remain_fasta,
        query_fasta = remain_fasta,
        out_dir = os.path.join(temp_dir, 'blast2'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_2_result, 
        gene_annotation=gene_annotation, 
        out_dir = temp_dir, 
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result)

    new_clusters = add_sample_pipeline.add_new_clusters(
        old_clusters = old_clusters,
        unmatched_clusters = unmatched_clusters,
        mcl_file=mcl_file
        )

    return new_clusters, new_combined_faa



def run_main_pipeline(args):
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

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')
    
    # data preparation
    gene_annotation = {}
    gene_position = {}
    data_preparation.extract_proteins(
        samples=samples,
        out_dir=collection_dir,
        gene_annotation = gene_annotation,
        gene_position = gene_position,
        table=args.table,
        annotate=annotate,
        threads=threads)

    number = args.number
    if number == 0:
        subset = samples[0:]
        remain = []
    else:
        subset = samples[0:number]
        remain = samples[number:]

    # run a subset of collection
    subset_dir = os.path.join(temp_dir, 'subset')
    if not os.path.exists(subset_dir):
        os.makedirs(subset_dir)
    
    subset_combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=subset_dir,
        samples=subset)

    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=subset_combined_faa,
        out_dir=subset_dir,
        threads=threads)

    blast_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(subset_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result, 
        gene_annotation=gene_annotation, 
        out_dir = subset_dir, 
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    subset_mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = subset_dir,
        blast_result = filtered_blast_result)

    subset_inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_clusters,
        mcl_file=subset_mcl_file)

    # run the remain of collection
    if len(remain) == 0:
        inflated_clusters = subset_inflated_clusters
        remain_combined_faa = ""
    else:
        remain_dir = os.path.join(temp_dir, 'remain')
        if not os.path.exists(remain_dir):
            os.makedirs(remain_dir)
        
        subset_representative_fasta = output.create_representative_fasta(
            clusters=subset_inflated_clusters, 
            gene_annotation=gene_annotation, 
            faa_fasta=cd_hit_represent_fasta, 
            out_dir=subset_dir)
        
        inflated_clusters, remain_combined_faa = add_samples(
            temp_dir=remain_dir, 
            new_samples=remain, 
            old_represent_faa=subset_representative_fasta,
            old_clusters=subset_inflated_clusters, 
            gene_annotation=gene_annotation, 
            collection_dir=collection_dir, 
            threads=threads, 
            args=args)


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_annotation=gene_annotation,
        gene_position=gene_position,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )

    if annotate == False:
        annotated_clusters = post_analysis.annotate_cluster_1(
            unlabeled_clusters=split_clusters, 
            gene_annotation=gene_annotation)
    else:
        annotated_clusters = post_analysis.annotate_cluster_2(
            unlabeled_clusters=split_clusters)
    
    output.create_outputs(gene_annotation,annotated_clusters,samples,collection_dir)

    if args.alignment != None:
        post_analysis.run_gene_alignment(annotated_clusters, gene_annotation, samples, collection_dir, args.alignment, threads)

    # output for next run
    rep_temp_file = os.path.join(temp_dir, 'representative_temp')
    os.system(f'cat {subset_combined_faa} {remain_combined_faa} > {rep_temp_file}')
    output.create_representative_fasta(
        clusters=split_clusters, 
        gene_annotation=gene_annotation, 
        faa_fasta=rep_temp_file, 
        out_dir=collection_dir)
    
    output.export_gene_annotation(gene_annotation, collection_dir)
    json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    json.dump(split_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)

    # shutil.rmtree(temp_dir)
        
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')

    

def run_add_sample_pipeline(args):
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
    
    # Check required files
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
    
    # data preparation
    data_preparation.extract_proteins(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation = gene_annotation,
        gene_position = gene_position,
        table=args.table,
        annotate=annotate,
        threads=threads
        )

    inflated_clusters, new_combined_faa = add_samples(
        temp_dir=temp_dir, 
        new_samples=new_samples, 
        old_represent_faa=old_represent_faa,
        old_clusters=old_clusters,
        gene_annotation=gene_annotation, 
        collection_dir=collection_dir, 
        threads=threads, 
        args=args)

    # post analysis
    new_samples.extend(old_samples)
    new_samples.sort(key= lambda x:x['id'])
    
    split_clusters = post_analysis.split_paralogs(
        gene_annotation=gene_annotation,
        gene_position=gene_position,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters, 
        gene_annotation=gene_annotation)

    output.create_outputs(gene_annotation,annotated_clusters,new_samples,collection_dir)

    if args.alignment != None:
        samples_dir = os.path.join(collection_dir, 'samples')
        if not os.path.exists(samples_dir):
            raise Exception(f'{samples_dir} does not exist')
        
        post_analysis.run_gene_alignment(annotated_clusters, gene_annotation, new_samples, collection_dir, args.alignment, threads)

    # output for next run
    output.export_gene_annotation(gene_annotation, collection_dir)
    json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(new_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    json.dump(split_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    
    rep_temp_file = os.path.join(temp_dir, 'representative_temp')
    os.system(f'cat {old_represent_faa} {new_combined_faa} > {rep_temp_file}')
    output.create_representative_fasta(
        clusters=split_clusters, 
        gene_annotation=gene_annotation, 
        faa_fasta=rep_temp_file, 
        out_dir=collection_dir)   
    
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
    main_cmd.set_defaults(func=run_main_pipeline)
    main_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    main_cmd.add_argument('-b', '--fasta', help='assembly input files',default=None, nargs='*', type=str)
    main_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    main_cmd.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    main_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
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

    add_cmd = subparsers.add_parser(
        'add',
        description='Add pipeline: add sample into previous collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=run_add_sample_pipeline)
    add_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    add_cmd.add_argument('-b', '--fasta', help='assembly input files',default=None, nargs='*', type=str)
    add_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    add_cmd.add_argument('-c', '--collection-dir', help='previous collection directory', required=True, type=str)
    add_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    add_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    add_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.95, type=float)
    add_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0, type=float)
    add_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    add_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0, type=float)
    add_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    add_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--table', help='codon table', default=11, type=int)
    add_cmd.add_argument('-a', '--alignment', help='run alignment for each gene cluster', default=None, action='store', choices=['protein', 'nucleotide'])

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()