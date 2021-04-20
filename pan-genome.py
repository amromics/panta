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
        report['summary'] = os.path.join(pan_genome_folder,'summary_statistics.txt')
        return report

    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    report = data_preparation.extract_proteins(report, timing_log=timing_log)
    report = data_preparation.combine_proteins(report, timing_log=timing_log)

    report = main_pipeline.run_cd_hit_iterative(report, threads=threads, timing_log=timing_log)
    report = main_pipeline.all_against_all_blast(report, threads=threads, timing_log=timing_log)
    report = main_pipeline.cluster_with_mcl(report, threads=threads, timing_log=timing_log)
    report = main_pipeline.reinflate_clusters(report)

    report = post_analysis.split_paralogs(report)
    report = post_analysis.label_cluster(report)
    report = post_analysis.annotate_cluster(report)

    report = output.create_spreadsheet(report)
    report = output.create_rtab(report)
    report = output.create_summary(report)
    report = output.create_representative_fasta(report)

    json.dump(report, open(os.path.join(report['pan_genome'], 'report.json'), 'w'), indent=4, sort_keys=True)


def run_add_sample_pipeline(args):
    report = {}
    report['samples'] = []
    for path in args.inputs:
        base_name = os.path.basename(path)
        sample_id = base_name.split('.')[0]
        sample = {'id':sample_id, 'input_file':path}
        report['samples'].append(sample)


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

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
