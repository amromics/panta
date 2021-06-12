import os,shutil
import logging
import multiprocessing
from Bio import SeqIO
import csv
import pandas as pd
import json
import gzip
from amromics.utils.command import run_command
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_roary(gff_folder,overwrite=False,threads=0, base_dir='.', timing_log=None):
    """
        Run roay make pangeome analysis (using prokka results in previous step)
        :param read_data: result holder
        :param base_dir: working directory
        :param threads: number of core CPU
        :return:
    """

    roary_folder=os.path.join(base_dir,'pangenome/roary')
    temp_folder = os.path.join(base_dir, 'pangenome/temp_roary')
    roary_output = os.path.join(roary_folder, 'summary_statistics.txt')
    if os.path.isfile(roary_output) and (not overwrite):
        logger.info('roary has run and the input has not changed, skip roarying')
        return roary_folder
    if not os.path.isdir(temp_folder):
        os.makedirs(temp_folder)
    gff_list = []
    for filename in os.listdir(gff_folder):
        if filename.endswith('.gz'):
            sample_id = filename.replace('.gff.gz','')
            #gffgz_file = os.path.join(sample['annotation'], sample_id + '.gff.gz')
            gff_file = os.path.join(temp_folder, sample_id + '.gff')
            if run_command('gunzip -c {} > {}'.format(os.path.join(gff_folder, filename), gff_file)) != 0:
                raise Exception('Cannot get {}'.format(os.path.join(gff_folder, filename)))
            gff_list.append(gff_file)
        else:
            gff_list.append(os.path.join(gff_folder, filename))

    # Make sure the directory is not there or roary will add timestamp
    if os.path.isfile(roary_folder):
        os.remove(roary_folder)
    if os.path.exists(roary_folder):
        shutil.rmtree(roary_folder)
    cmd = 'roary -p {} -f {} -v '.format(threads, roary_folder) + ' '.join(gff_list)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('roary fail to run!')

    cmd = 'gzip ' + os.path.join(roary_folder, 'gene_presence_absence.csv')
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running {}'.format(cmd))

    shutil.rmtree(temp_folder)

    return roary_folder
