#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    The entry point
"""
from __future__ import division, print_function, absolute_import

import argparse
import logging
import multiprocessing
import os
import sys
import pandas as pd

from amromics.pipeline.analysis import single_genome_analysis, pan_genome_analysis
from amromics.utils import valid_id, software_version
from amromics import __version__

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s [%(name)s] %(levelname)s : %(message)s')
logger = logging.getLogger(__name__)


def version_func(args):
    software_version([
        'python',
        'blast',
        'trimmomatic',
        'spades',
        'shovill',
        'prokka',
        'mlst',
        'abricate',
        'roary',
        'parsnp'
    ])


def is_true(t):
    return t.lower() in ['y', 'yes', '1', 'true', 't', 'ok']


def input_file_to_samples(input_file, sep='\t'):
    sample_report = []
    sample_df = pd.read_csv(input_file, sep=sep, dtype='str')
    sample_df.fillna('', inplace=True)

    if 'trim' not in sample_df.columns:
        sample_df['trim'] = False
    else:
        sample_df['trim'] = sample_df['trim'].apply(lambda x: is_true(x))

    if 'strain' not in sample_df.columns:
        sample_df['strain'] = ''

    if 'metadata' not in sample_df.columns:
        sample_df['metadata'] = ''

    # 1. validate input
    for i, row in sample_df.iterrows():
        sample_id = row['sample_id']
        if not valid_id(sample_id):
            raise Exception(
                '{} invalid: sample ID can only contain alpha-numerical or underscore charactors'.format(
                    sample_id))
        input_files = row['files'].split(';')
        input_files = [input_file.strip() for input_file in input_files]
        if len(input_files) <= 0:
            raise Exception(
                'No input file for sample {}'.format(sample_id))

        for input_file in input_files:
            if not os.path.isfile(input_file):
                raise Exception(
                    'Input file {} (sample {}) not found!'.format(input_file, sample_id))

        metadata = row['metadata'].split(';')
        mt = {}
        if len(metadata) > 0:
            for kv in metadata:
                if len(kv.split(':')) == 2:
                    k, v = kv.split(':')
                    mt[k] = v
        sample = {
            'id': sample_id,
            'name': row['sample_desc'].strip(),
            'input_type': row['input_type'].strip(),
            'files': ';'.join(input_files),  # Re-join to make sure no white characters slipped in
            'genus': row['genus'].strip(),
            'species': row['species'].strip(),
            'strain': row['strain'].strip(),
            'trim': row['trim'],
            'metadata': mt,
            'updated': False,
        }
        sample_report.append(sample)
    return sample_report


def single_genome_clean_func(args):
    """
    Remove a single genome from the file structure, report any collections that might be affected
    Parameters
    ----------
    args

    Returns
    -------

    """
    pass


def collection_clean_func(args):
    pass


def single_genome_analysis_func(args):
    """
    Parse command line and call single genome analysis
    Parameters
    ----------
    args

    Returns
    -------

    """
    work_dir = args.work_dir
    threads = args.threads
    memory = args.memory
    timing_log = args.time_log
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    # run single sample pipeline
    samples = input_file_to_samples(args.input)
    single_genome_analysis(samples, work_dir, threads, memory, timing_log)
    return samples


def pan_genome_analysis_func(args):
    """
    Parse command line and call single genome analysis and pan-genome analysis
    """
    collection_id = args.collection_id
    collection_name = args.collection_name
    if not collection_name:
        collection_name = collection_id

    work_dir = args.work_dir
    threads = args.threads
    memory = args.memory
    timing_log = args.time_log
    overwrite = False

    if not valid_id(collection_id):
        raise Exception('{} invalid: collection ID can only contain alpha-numerical or underscore charactors'.format(collection_id))

    if threads <= 0:
        threads = multiprocessing.cpu_count()

    samples = input_file_to_samples(args.input)

    # First run single analysis
    samples = single_genome_analysis(samples, work_dir, threads=threads, memory=memory, timing_log=timing_log)
    report = pan_genome_analysis(
        samples, work_dir,
        collection_id, collection_name, overwrite=overwrite,
        threads=threads, memory=memory, timing_log=timing_log)
    logger.info('Congratulations, collection {} is done!'.format(collection_id))


def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='amromics',
        description='Tool for managing and analyzing antibiotic resistant bacterial datasets')
    parser.add_argument('-V', '--version', action='version', version=__version__)

    subparsers = parser.add_subparsers(title='sub command', help='sub command help')
    version_cmd = subparsers.add_parser(
        'dep', description='Check dependency versions',
        help='Print versions dependencies if exist',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    version_cmd.set_defaults(func=version_func)

    pa_cmd = subparsers.add_parser(
        'pg',
        description='Pan-genome analysis of a collection',
        help='Pan-genome analysis of a collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa_cmd.set_defaults(func=pan_genome_analysis_func)
    pa_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    pa_cmd.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    pa_cmd.add_argument('-c', '--collection-id', help='Collection ID', required=True, type=str)
    pa_cmd.add_argument('-n', '--collection-name', help='Collection name', type=str, default='')
    pa_cmd.add_argument('-i', '--input', help='Input file', required=True, type=str)
    pa_cmd.add_argument('--work-dir', help='Working directory', default='data/work')
    pa_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)

    pa_cmd = subparsers.add_parser(
        'sg',
        description='Single-genome analysis',
        help='Single-genome analysis of a collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pa_cmd.set_defaults(func=single_genome_analysis_func)
    pa_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    pa_cmd.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    pa_cmd.add_argument('-i', '--input', help='Input file', required=True, type=str)
    pa_cmd.add_argument('--work-dir', help='Working directory', default='data/work')
    pa_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)

    args = parser.parse_args(arguments)
    return args.func(args)


if __name__ == "__main__":
    main()

