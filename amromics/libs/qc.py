import os
import shutil
import csv
import logging
import multiprocessing
import gzip
from amromics.utils.command import run_command
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def qc_reads(prefix_name, reads,base_dir = '.', threads=0, timing_log=None, **kargs):
    """
        Run QC process for pair-end input using fastqc and multiqc
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to qc output file
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    out_fastqc = os.path.join(base_dir, prefix_name + '_fastqc')
    if not os.path.exists(out_fastqc):
        os.makedirs(out_fastqc)
    cmd = 'fastqc -t {threads} -o {outdir} {pe1} {pe2}'.format(threads=threads, outdir=out_fastqc, pe1=reads['pe1'], pe2=reads['pe2'])
    fastqc_ret = run_command(cmd, timing_log)

    if fastqc_ret != 0:
        return None
    out_multiqc = os.path.join(base_dir, prefix_name + '_multiqc')
    if not os.path.exists(out_multiqc):
        os.makedirs(out_multiqc)
    cmd = 'multiqc -o {outdir} {indir}'.format(outdir=out_multiqc, indir=out_fastqc)
    multiqc_ret = run_command(cmd, timing_log)
    if multiqc_ret != 0:
        return None
    #read_data['fastqc']=out_fastqc

    return os.path.join(out_multiqc,'multiqc_data','multiqc_fastqc.txt')

def assembly_eval(prefix_name,assembly, base_dir = '.', threads=0, timing_log=None, **kargs):
    """
        Run QC process for assembly output
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to report file
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT


    out_quast = os.path.join(base_dir, prefix_name + '_quast')
    if not os.path.exists(out_quast):
        os.makedirs(out_quast)

    cmd = 'quast.py -t {threads} -o {outdir} {input}'.format(threads=threads, outdir=out_quast, input=assembly)
    quast_ret = run_command(cmd, timing_log)

    if quast_ret != 0:
        return None


    return os.path.join(out_quast,'report.tsv')

def map_reads_to_assembly_bwamem(prefix_name,assembly,reads, base_dir = '.', threads=0, memory=50, timing_log=None, **kargs):
    if not os.path.isfile(assembly + '.sa'):
        cmd = 'bwa index ' +assembly
        ret = run_command(cmd, timing_log)
        if ret != 0:
            return None
    path_out = os.path.dirname(assembly)

    cmd = 'bwa mem -t {threads} {index}'.format(threads=threads, index=assembly)
    if 'pe1' in reads and 'pe2' in reads:
        pe_sam = os.path.join(path_out, prefix_name + '_pe.sam')
        pe_bam = os.path.join(path_out, prefix_name + '_pe.bam')
        cmd_bwa_pe = cmd + ' ' + reads['pe1'] + ' ' + reads['pe2'] + ' > ' + pe_sam
        run_command(cmd_bwa_pe, timing_log)
        cmd_st_pe = 'samtools view -u {sam} | samtools sort -@{threads} -o {bam} - ;samtools index {bam}'.format(
            sam=pe_sam,
            threads=threads,
            bam=pe_bam
        )
        run_command(cmd_st_pe, timing_log)
        return pe_bam

    if 'se' in reads:
        se_sam = os.path.join(path_out, prefix_name + '_se.sam')
        se_bam = os.path.join(path_out, prefix_name + '_se.bam')
        cmd_bwa_se = cmd + ' ' + reads['se'] + ' > ' + se_sam
        run_command(cmd_bwa_se, timing_log)
        cmd_st_se = 'samtools view -u {sam} | samtools sort -@{threads} -o {bam} - ;samtools index {bam}'.format(
            sam=se_sam,
            threads=threads,
            bam=se_bam
        )
        run_command(cmd_st_se, timing_log)
        return se_bam
