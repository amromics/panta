import os
import shutil
import csv
import logging
import gzip
import glob
import multiprocessing
from amromics.utils.command import run_command
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def annotate_prokka(prefix_name,assembly,genus=None,species=None, strain=None,gram=None, base_dir='.',  overwrite=False,timing_log=None, threads=0):
    """
        Run annotation process using prokka
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, prefix_name+'_prokka' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    annotation_gbk=  os.path.join(path_out, prefix_name + '.gbk.gz')
    annotation_gff= path_out+'/'+str(prefix_name)+'.gff.gz'
    annotation_faa = path_out+'/'+str(prefix_name)+'.faa.gz'
    annotation_ffn = path_out+'/'+str(prefix_name)+'.ffn.gz'
    annotation_fna = path_out+'/'+str(prefix_name)+'.fna.gz'
    if os.path.isfile(annotation_gff) and os.path.isfile(annotation_gbk) and (not overwrite):
        # Dont run again if gff/gbk file exists
        logger.info('GFF and GBK files found, skip annotating')
        return annotation_gff,annotation_faa,annotation_ffn,annotation_fna,annotation_gbk
    gunzip_fasta=assembly
    if assembly.endswith('.gz'):
        gunzip_fasta = os.path.join(path_out, prefix_name + '.fin')
        cmd = 'gunzip -c {} > {}'.format(assembly, gunzip_fasta)
        run_command(cmd)
    cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
    cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=prefix_name, path=path_out)
    if not genus ==None:
        cmd += ' --genus ' +genus
    if not species  ==None:
        cmd += ' --species ' + species
    if not strain  ==None:
        cmd += ' --strain ' + strain
    if not gram  ==None:
        cmd += ' --gram ' + gram
    cmd += ' ' + gunzip_fasta
    cmd = "bash -c '{}'".format(cmd)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Command {} returns non-zero ()!'.format(cmd, ret))

    for file_name in glob.glob(os.path.join(path_out, '*')):
        ext = file_name[-3:]
        if ext in ['gff', 'gbk', 'ffn','faa','fna']: # fna?
            run_command('gzip {}'.format(file_name))
        else:
            os.remove(file_name)

    return annotation_gff,annotation_faa,annotation_ffn,annotation_fna,annotation_gbk
