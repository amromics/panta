import argparse
import os
import logging
import json
from pan_genome import *

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s : %(message)s',
                    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)

def run_main_pipeline(args):
    collection_dir = args.collection_dir
    timing_log = args.time_log
    threads = args.threads
    overwrite = args.over_write
    dontsplit = args.dont_split
    
    report = {}
    report['samples'] = []
    for path in args.inputs:
        base_name = os.path.basename(path)
        sample_id = base_name.split('.')[0]
        sample = {'id':sample_id, 'input_file':path}
        report['samples'].append(sample)
    
    pan_genome_folder = os.path.join(collection_dir, 'pan_genome')
    temp_dir = os.path.join(collection_dir, 'temp')
    report['pan_genome'] = pan_genome_folder
    report['temp_dir'] = temp_dir
    
    # Check if pan-genome has run
    pan_genome_output = os.path.join(pan_genome_folder,'summary_statistics.txt')
    if os.path.isfile(pan_genome_output) and (not overwrite):
        return

    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    report = data_preparation.extract_proteins(
        report, 
        gene_annotation = {}, 
        timing_log=timing_log
        )
    report = data_preparation.combine_proteins(report, timing_log=timing_log)

    report = main_pipeline.run_cd_hit_iterative(report, threads=threads, timing_log=timing_log)
    
    report['blast_result_file'] = main_pipeline.all_against_all_blast(
        out_dir = os.path.join(report['temp_dir'], 'blast'),
        database_fasta = report['cd_hit_cluster_fasta'],
        query_fasta = report['cd_hit_cluster_fasta'],
        threads=threads, 
        timing_log=timing_log
        )

    report['uninflated_mcl_clusters'] = main_pipeline.cluster_with_mcl(
        out_dir = report['temp_dir'],
        blast_results = report['blast_result_file'],
        threads=threads, 
        timing_log=timing_log)
    report = main_pipeline.reinflate_clusters(report)

    if dontsplit == False:
        report = post_analysis.split_paralogs(report)
        report['labeled_clusters'] = post_analysis.label_cluster(unlabeled_clusters=report['split_clusters'])
    else:
        report['labeled_clusters'] = post_analysis.label_cluster(unlabeled_clusters=report['inflated_unsplit_clusters'])
    
    report = post_analysis.annotate_cluster(report)

    report = output.create_spreadsheet(report)
    report = output.create_rtab(report)
    report = output.create_summary(report)
    report = output.create_representative_fasta(report)

    json.dump(report['gene_annotation'], open(os.path.join(report['pan_genome'], 'gene_annotation.json'), 'w'), indent=4, sort_keys=True)
    json.dump(report['inflated_unsplit_clusters'], open(os.path.join(report['pan_genome'], 'unsplit_clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(report['samples'], open(os.path.join(report['pan_genome'], 'samples.json'), 'w'), indent=4, sort_keys=True)

def run_add_sample_pipeline(args):
    collection_dir = args.collection_dir
    timing_log = args.time_log
    threads = args.threads
    dontsplit = args.dont_split

    report = {}
    report['samples'] = []
    for path in args.inputs:
        base_name = os.path.basename(path)
        sample_id = base_name.split('.')[0]
        sample = {'id':sample_id, 'input_file':path}
        report['samples'].append(sample)

    pan_genome_folder = os.path.join(collection_dir, 'pan_genome')
    temp_dir = os.path.join(collection_dir, 'new_temp')
    report['pan_genome'] = pan_genome_folder
    report['temp_dir'] = temp_dir

    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Check if last collection exist
    representative_fasta = os.path.join(pan_genome_folder, 'representative.fasta')
    if not os.path.isfile(representative_fasta):
        raise Exception(f'{representative_fasta} is not exist')
    report['representative_fasta'] = representative_fasta

    gene_annotation = json.load(open(os.path.join(pan_genome_folder, 'gene_annotation.json'), 'r'))
    
    old_clusters = json.load(open(os.path.join(pan_genome_folder, 'unsplit_clusters.json'), 'r'))
    report = data_preparation.extract_proteins(
        report, 
        gene_annotation = gene_annotation, 
        timing_log=timing_log
        )
    report = data_preparation.combine_proteins(report, timing_log=timing_log)

    report = add_sample_pipeline.run_cd_hit_2d(report, threads=threads, timing_log=timing_log)

    report['blast_1_result_file'] = main_pipeline.all_against_all_blast(
        out_dir = os.path.join(report['temp_dir'], 'blast1'),
        database_fasta = representative_fasta,
        query_fasta = report['not_match_fasta_file'],
        threads=threads, 
        timing_log=timing_log
        )

    report['blast_remain_fasta'] = add_sample_pipeline.filter_fasta(
        blast_result = report['blast_1_result_file'], 
        fasta_file = report['not_match_fasta_file'], 
        out_dir = temp_dir
        )

    report['blast_2_result_file'] = main_pipeline.all_against_all_blast(
        out_dir = os.path.join(report['temp_dir'], 'blast2'),
        database_fasta = report['blast_remain_fasta'],
        query_fasta = report['blast_remain_fasta'],
        threads=threads, 
        timing_log=timing_log
        )
    
    report['mcl_clusters'] = main_pipeline.cluster_with_mcl(
        out_dir = report['temp_dir'],
        blast_results = report['blast_2_result_file'],
        threads=threads, 
        timing_log=timing_log)

    report['inflated_unsplit_clusters'] = add_sample_pipeline.reinflate_clusters(
        old_clusters=old_clusters, 
        cd_hit_2d_cluster=report['cd_hit_2d_cluster'], 
        blast_1_result_file=report['blast_1_result_file'], 
        mcl_clusters=report['mcl_clusters']
        )
    
    old_samples = json.load(open(os.path.join(pan_genome_folder, 'samples.json'), 'r'))
    report['samples'].extend(old_samples)
    
    if dontsplit == False:
        report = post_analysis.split_paralogs(report)
        report['labeled_clusters'] = post_analysis.label_cluster(unlabeled_clusters=report['split_clusters'])
    else:
        report['labeled_clusters'] = post_analysis.label_cluster(unlabeled_clusters=report['inflated_unsplit_clusters'])
    
    report = post_analysis.annotate_cluster(report)

    report = output.create_spreadsheet(report)
    report = output.create_rtab(report)
    report = output.create_summary(report)
    report = output.create_representative_fasta(report)

    json.dump(report['gene_annotation'], open(os.path.join(report['pan_genome'], 'gene_annotation.json'), 'w'), indent=4, sort_keys=True)
    json.dump(report['inflated_unsplit_clusters'], open(os.path.join(report['pan_genome'], 'unsplit_clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(report['samples'], open(os.path.join(report['pan_genome'], 'samples.json'), 'w'), indent=4, sort_keys=True)

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    main_cmd = subparsers.add_parser(
        'main',
        description='main pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    main_cmd.set_defaults(func=run_main_pipeline)
    main_cmd.add_argument('inputs', help='Input files', type=str, nargs='+')
    main_cmd.add_argument('-c', '--collection-dir', help='Collection directory', required=True, type=str)
    main_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    main_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)
    main_cmd.add_argument('-w', '--over-write', help='over write', default=False, action='store_true')
    main_cmd.add_argument('-s', '--dont-split', help='dont-split', default=False, action='store_true')

    add_cmd = subparsers.add_parser(
        'add',
        description='add sample pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=run_add_sample_pipeline)
    add_cmd.add_argument('inputs', help='Input files', type=str, nargs='+')
    add_cmd.add_argument('-c', '--collection-dir', help='Collection directory', required=True, type=str)
    add_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)
    add_cmd.add_argument('-s', '--dont-split', help='dont-split', default=True, action='store_false')

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
