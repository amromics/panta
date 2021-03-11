#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    The entry point
"""
from __future__ import division, print_function, absolute_import

import argparse
import json
import logging
import multiprocessing
import os
import shutil
import socket
import sys
import pandas as pd

from amromics import extract_json
from amromics.pipeline import wrapper
from amromics.utils import valid_id, software_version

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s [%(name)s] %(levelname)s : %(message)s')
logger = logging.getLogger(__name__)

__version__ = '0.9.0'


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


def start_server_func(args):
    """
    Start the server at the specified  port
    """

    port = args.port
    webapp_dir = args.webapp_dir

    # Check if the port is open
    a_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    location = ("127.0.0.1", port)
    result_of_check = a_socket.connect_ex(location)
    a_socket.close()
    if result_of_check != 0:
        logger.info('Starting server on port {}'.format(port))
        cmd = 'cd {} && live-server --no-browser --port={}  --entry-file=index.html'.format(webapp_dir, port)
        ret = wrapper.run_command(cmd)
        if ret != 0:
            sys.exit(ret)
    else:
        logger.error('Port {} is not available, server cannot start!'.format(port))


def collection_pa_func(args):
    """

    """
    collection_id = args.collection_id
    collection_name = args.collection_name
    if not collection_name:
        collection_name = collection_id

    work_dir = args.work_dir
    webapp_dir = args.webapp_dir
    webapp_static_dir = os.path.join(webapp_dir, 'static')
    if not os.path.exists(webapp_static_dir):
        raise Exception('Webapp directory {} not available'.format(webapp_dir))
    webapp_data_dir = os.path.join(webapp_static_dir, 'data')
    if not os.path.exists(webapp_data_dir):
        os.makedirs(webapp_data_dir)

    threads = args.threads
    memory = args.memory
    timing_log = args.time_log
    overwrite = False

    if not valid_id(collection_id):
        raise Exception('{} invalid: collection ID can only contain alpha-numerical or underscore charactors'.format(collection_id))

    if threads <= 0:
        threads = multiprocessing.cpu_count()

    report = {'collection_id': collection_id,
              'samples': []}
    sample_df = pd.read_csv(args.input, sep='\t', dtype='str')
    sample_df.fillna('', inplace=True)

    if 'trim' not in sample_df.columns:
        sample_df['trim'] = False
    else:
        sample_df['trim'] = sample_df['trim'].apply(lambda x: x.lower() in ['y', 'yes', '1', 'true', 't', 'ok'])

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
                'No input file for sample {}'.format( sample_id))

        for input_file in input_files:
            if not os.path.isfile(input_file):
                raise Exception(
                    'Input file {} (sample {}) not found!'.format(input_file, sample_id))

        if i == 0 and row['input_type'].startswith('ass'):
            # TODO Check if this sequence can be the reference for parsnp (Quang)
            pass

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
        report['samples'].append(sample)

    # run single sample pipeline
    for sample in report['samples']:
        sample_id = sample['id']
        sample_dir = os.path.join(work_dir, 'samples', sample_id)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        wrapper.run_single_sample(
            sample, sample_dir=sample_dir, threads=threads,
            memory=memory, timing_log=timing_log)
    collection_dir = os.path.join(work_dir, 'collections', collection_id)

    dataset_sample_ids = []
    for sample in report['samples']:
        sample_id = sample['id']
        dataset_sample_ids.append(sample_id)
        overwrite = overwrite or sample['updated']

    if not os.path.exists(collection_dir):
        os.makedirs(collection_dir)

    # to check if the set of samples has not changed
    dataset_sample_ids = sorted(dataset_sample_ids)
    sample_set_file = os.path.join(collection_dir, 'sample_set.json')
    if os.path.isfile(sample_set_file):
        with open(sample_set_file) as fn:
            sample_set = json.load(fn)
        if sample_set != dataset_sample_ids:
            overwrite = True

    if overwrite:
        roary_folder = os.path.join(collection_dir, 'roary')
        if os.path.exists(roary_folder):
            shutil.rmtree(roary_folder)

        phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
        if os.path.exists(phylogeny_folder):
            shutil.rmtree(phylogeny_folder)
        # Note: Check for existing alignment is done within

    # Write the set of sample IDs
    with open(sample_set_file, 'w') as fn:
        json.dump(dataset_sample_ids, fn)

    report = wrapper.run_roary(report, collection_dir=collection_dir, threads=threads, overwrite=overwrite, timing_log=timing_log)
    report = wrapper.run_phylogeny_iqtree(report, collection_dir=collection_dir, threads=threads, overwrite=overwrite, timing_log=timing_log)
    report = wrapper.run_gene_phylogeny_parellel(report, collection_dir=collection_dir, threads=threads, overwrite=overwrite, timing_log=timing_log)
    with open(os.path.join(collection_dir, collection_id + '_dump.json'), 'w') as fn:
        json.dump(report, fn)

    # # clean up
    # if os.path.exists(collection_dir + "/temp"):
    #     shutil.rmtree(collection_dir + "/temp")

    extract_json.export_json(work_dir, webapp_data_dir,
                             collection_id, collection_name)
    logger.info('Congratulations, collection {} is imported to web-app!'.format(collection_id))


def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='amrviz',
        description='Tool for managing and analyzing antibiotic resistant bacterial datasets')
    parser.add_argument('-V', '--version', action='version', version=__version__)

    subparsers = parser.add_subparsers(title='sub command', help='sub command help')
    version_cmd = subparsers.add_parser(
        'dep', description='Check dependency versions',
        help='Print versions dependencies if exist',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    version_cmd.set_defaults(func=version_func)

    pa_cmd = subparsers.add_parser(
        'pa', description='NGS analysis pipeline', help='NGS analysis pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa_cmd.set_defaults(func=collection_pa_func)
    pa_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    pa_cmd.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    pa_cmd.add_argument('-c', '--collection-id', help='Collection ID', required=True, type=str)
    pa_cmd.add_argument('-n', '--collection-name', help='Collection name', type=str, default='')
    pa_cmd.add_argument('-i', '--input', help='Input file', required=True, type=str)
    pa_cmd.add_argument('--work-dir', help='Working directory', default='data/work')
    pa_cmd.add_argument('--webapp-dir', help='Webapp directory', default='web-app')
    pa_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)

    start_cmd = subparsers.add_parser(
        'start', description='Start amr_viz server', help='Start server',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    start_cmd.set_defaults(func=start_server_func)
    start_cmd.add_argument('-p', '--port', help='The port the server is running on', default=3000, type=int)
    start_cmd.add_argument('--webapp-dir', help='Webapp directory', default='web-app')

    args = parser.parse_args(arguments)
    return args.func(args)


if __name__ == "__main__":
    main()

