# -*- coding: utf-8 -*-
import json
import os
import shutil
from Bio import SeqIO
import datetime
import pandas as pd
import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s [%(name)s] %(levelname)s : %(message)s')
logger = logging.getLogger(__name__)

# TODOs:
# - Can make it faster with using fasttree (parsnps need to have this option specifically set
# - By default, parsnp use bootstrap of 1000. See if we can change the value and get the boottrap values
# - Can provide the genbank of the reference (using prokka annotation)


def run_command(cmd, timing_log=None):
    """
    Run a command line, return the returning code of the command
    :param cmd:
    :param timing_log:
    :return:
    """
    if timing_log is not None:
        cmd = '/usr/bin/time --append -v -o {} bash -c "{}"'.format(timing_log, cmd)
    logger.info('Run "' + cmd + '"')
    ret = os.system(cmd)
    return ret


def assemble_shovill(sample, sample_dir, threads=4, memory=50, overwrite=False, timing_log=None):
    """
    Run assembly using shovill

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    memory: float
        maximum memory used
    timing_log: str
        log file
    Returns str
    -------
        Path to the assembly
    """
    sample_id = sample['id']
    path_out = os.path.join(sample_dir, 'assembly')
    assembly_file = os.path.join(path_out, sample_id + '_contigs.fasta')
   
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    elif os.path.isfile(assembly_file) and (not overwrite):
        return assembly_file

    cmd = 'shovill  --ram {memory} --cpus {threads} --outdir {path_out}'.format(
        memory=int(memory), threads=threads, path_out=path_out)

    if 'trim' in sample and sample['trim']:
        cmd += ' --trim'

    pe_files = sample['files'].split(';')
    if len(pe_files) > 0:
        pe1 = pe_files[0]
        pe2 = pe_files[1]
        cmd += ' --force  --R1 {pe1} --R2 {pe2}'.format(pe1=pe1, pe2=pe2)
    else:
       raise Exception('Only support pair-end reads!')

    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running shovill!')

    # Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fa'), "fasta"))
    contigs = sorted(contigs, key=len, reverse=True)
    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id = sample['id']+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")

    sample['updated'] = True
    return assembly_file


def get_assembly(sample, sample_dir, overwrite=False):
    """
    Get the assembly from user input

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    Returns
    -------
        path to assembly file
    """
    path_out = os.path.join(sample_dir, 'assembly')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    contigs = list(SeqIO.parse(sample['files'], "fasta"))
    assembly_file = os.path.join(path_out, sample['id'] + '_contigs.fasta')
    if os.path.isfile(assembly_file) and (not overwrite):
        return assembly_file

    contigs = sorted(contigs, key=len, reverse=True)

    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id = sample['id'] + '_C' + str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")

    sample['updated'] = True
    return assembly_file


def annotate_prokka(sample, sample_dir,  threads=8, overwrite=False, timing_log=None):
    """
    Annotate the sample using prokka

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing annotation
    timing_log: str
        log file
    Returns
    -------
        path to prokka folder
    """

    path_out = os.path.join(sample_dir, 'prokka')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    gff_file_out = os.path.join(path_out, sample['id'] + '.gff')
    gbk_file_out = os.path.join(path_out, sample['id'] + '.gbk')

    if os.path.isfile(gff_file_out) and os.path.isfile(gbk_file_out) and (not overwrite):
        # Dont run again if gff/gbk file exists
        logger.info('GFF and GBK files found, skip annotating')
        return path_out

    cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
    cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=sample['id'], path=path_out)
    if sample['genus']:
        cmd += ' --usegenus --genus ' + sample['genus']
    if sample['species']:
        cmd += ' --species ' + sample['species']
    if sample['strain']:
        cmd += ' --strain ' + sample['strain']

    # Disable this for now so that we dont have to install signalp
    # if sample['gram']:
    #    cmd += ' --gram ' + sample['gram']
    cmd += ' ' + sample['assembly']
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Command {} returns non-zero ()!'.format(cmd, ret))

    for ext in ['err', 'faa', 'fsa', 'log', 'sqn', 'tbl', 'tsv', 'txt']:
        file_name = os.path.join(path_out, sample['id'] + '.' + ext)
        if os.path.isfile(file_name):
            os.remove(file_name)

    sample['updated'] = True
    return path_out


def mlst(sample, sample_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run mlst

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to mlst result file
    """
    path_out = os.path.join(sample_dir, 'mlst')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    mlst_out = os.path.join(path_out, sample['id'] + '_mlst.tsv')
    if os.path.isfile(mlst_out) and (not overwrite):
        logger.info('MLST for {} exists, skip mlsting'.format(sample['id']))
        return mlst_out

    cmd = 'mlst --quiet --threads {threads} --nopath {infile} > {outfile}'.format(
        threads=threads,
        infile=sample['assembly'],
        outfile=mlst_out)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mlst')
    sample['updated'] = True
    return mlst_out


def detect_amr(sample, sample_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run abricate to identify resistant genes

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to resistant gene file
    """
    path_out = os.path.join(sample_dir, 'abricate')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    # TODO: replace by consensus db later
    amr_out = os.path.join(path_out, sample['id'] + '_resistome.tsv')
    if os.path.isfile(amr_out) and (not overwrite):
        logger.info('Resistome for {} exists, skip analysis'.format(sample['id']))
        return amr_out

    cmd = 'abricate --quiet --threads {threads} --nopath --db card {infile} > {outfile}'.format(
        threads=threads,
        infile=sample['assembly'],
        outfile=amr_out)
    if run_command(cmd, timing_log) != 0:
        raise Exception('Error running amr')
    sample['updated'] = True
    return amr_out


def detect_virulome(sample, sample_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run abricate to identify virulent genes using VFDB

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to virulent gene file
    """
    path_out = os.path.join(sample_dir, 'abricate')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    vir_out = os.path.join(path_out, sample['id'] + '_virulome.tsv')
    if os.path.isfile(vir_out) and (not overwrite):
        logger.info('Virulome for {} exists, skip analysis'.format(sample['id']))
        return vir_out

    cmd = 'abricate --quiet --threads {threads} --nopath --db vfdb {infile} > {outfile}'.format(
        threads=threads,
        infile=sample['assembly'],
        outfile=vir_out)
    if run_command(cmd, timing_log) != 0:
        raise Exception('Error running virulent detection')
    sample['updated'] = True
    return vir_out


def detect_plasmid(sample, sample_dir,  threads=8, overwrite=False, timing_log=None):
    """
    Detect plasmids

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to origin of replication file
    """

    path_out = os.path.join(sample_dir, 'abricate')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    # Plasmid finder
    oriREP_out = os.path.join(path_out, sample['id'] + '_plasmid.tsv')
    if os.path.isfile(oriREP_out) and (not overwrite):
        logger.info('ORI for {} exists, skip analysis'.format(sample['id']))
        return oriREP_out

    cmd = 'abricate --quiet --threads {threads} --nopath --db plasmidfinder {infile} > {outfile}'.format(
        threads=threads, infile=sample['assembly'], outfile=oriREP_out)
    if run_command(cmd, timing_log) != 0:
        return None

    sample['updated'] = True
    return oriREP_out


def run_single_sample(sample, sample_dir, threads=8, memory=50, timing_log=None):
    """
    Run the pipeline to analyse one single sample. Steps:
        0. Assembly if the input are read data
        1. Annotate the assembly using prokka
        2. Identify resistant, virulent genes, origins of replication, and plasmids

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    timing_log: str
        log file
    Returns
    -------
        the sample object
    """

    sample['execution_start'] = str(datetime.datetime.now())

    if sample['input_type'] not in ['asm', 'assembly']:
        sample['assembly'] = assemble_shovill(sample, sample_dir=sample_dir, threads=threads, memory=memory, timing_log=timing_log)
    else:
        sample['assembly'] = get_assembly(sample, sample_dir=sample_dir)

    sample['annotation'] = annotate_prokka(
        sample, sample_dir=sample_dir,  threads=threads, timing_log=timing_log)
    sample['mlst'] = mlst(sample, sample_dir=sample_dir, threads=threads)
    sample['resistome'] = detect_amr(
        sample, sample_dir=sample_dir, threads=threads, timing_log=timing_log)
    sample['virulome'] = detect_virulome(
        sample, sample_dir=sample_dir, threads=threads, timing_log=timing_log)
    sample['plasmid'] = detect_plasmid(
        sample, sample_dir=sample_dir, threads=threads, timing_log=timing_log)
    sample['execution_end'] = str(datetime.datetime.now())
    return sample


def run_roary(report, collection_dir='.', threads=8, overwrite=False, timing_log=None):
    """
    Run roary for pangenome analysis. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """

    gff_list = []
    for sample in report['samples']:
        sample_id = sample['id']
        gff_file = os.path.join(sample['annotation'], sample_id + '.gff')
        assert os.path.isfile(gff_file)
        gff_list.append(gff_file)
        # If at least one sample has changed, we need to re-run roary
        overwrite = overwrite or sample['updated']

    roary_folder = os.path.join(collection_dir, 'roary')
    roary_output = os.path.join(roary_folder, 'core_alignment_header.embl')

    # Check if roary has run for the same dataset ID and the same set of samples
    report['roary'] = roary_folder
    if os.path.isfile(roary_output) and (not overwrite):
        logger.info('roary has run and the input has not changed, skip roarying')
        return report

    # Make sure the directory is not there or roary will add timestamp
    if os.path.isfile(roary_folder):
        os.remove(roary_folder)
    if os.path.exists(roary_folder):
        shutil.rmtree(roary_folder)

    cmd = 'roary -p {} -f {} -e -n -v '.format(threads, roary_folder) + ' '.join(gff_list)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('roary fail to run!')
    return report


def run_phylogeny(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run parsnp to create phylogeny tree. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    genome_dir = os.path.join(collection_dir, 'temp/fasta')
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)

    report['phylogeny'] = phylogeny_folder
    phylogeny_file = os.path.join(phylogeny_folder, 'parsnp.tree')
    if os.path.isfile(phylogeny_file) and (not overwrite):
        logger.info('phylogeny tree exists and input has not changed, skip phylogeny analysis')
        return report

    reference_genome = None
    sample_list = []
    for i, sample in enumerate(report['samples']):
        if i == 0:
            reference_genome = sample['assembly']
        else:
            sample_list.append(sample['assembly'])

    cmd = 'parsnp -r {} -d {} -o {} -p {}'.format(
        reference_genome,
        ' '.join(sample_list),
        phylogeny_folder, threads)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parsnp')
    return report


def run_alignment(report, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run phylogenetic analysis of gene clusters. If the list of samples has not changed, and
    none of the samples has changed, the existing tree will be kept unless overwrite is
    set to True

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    gene_cluster_file = report['roary'] + '/gene_presence_absence.csv'
    dict_cds = {}
    for sample in report['samples']:
        for seq in SeqIO.parse(os.path.join(sample['annotation'], sample['id'] + '.ffn'), 'fasta'):
            dict_cds[seq.id] = seq

    # make folder contains sequences for each gene
    alignment_dir = os.path.join(collection_dir, 'alignments')
    gene_df = pd.read_csv(gene_cluster_file, dtype=str)
    gene_df.fillna('', inplace=True)

    sample_columns = list(gene_df.columns)[14:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_dir = os.path.join(alignment_dir, gene_id)
        gene_file_dir = os.path.join(gene_dir, 'files')

        if not os.path.exists(gene_file_dir):
            os.makedirs(gene_file_dir)
        gene_files = []
        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                gene_file = os.path.join(gene_file_dir, row[sample_column] + '.fasta')
                SeqIO.write(dict_cds[row[sample_column]], gene_file, 'fasta')
                gene_files.append(gene_file)
                gene_list.append(row[sample_column])
                # TODO: make sure all samples in this gene have not updated

        gene_list = sorted(gene_list)
        # Only analyse if there are more than 3 genes
        if len(gene_list) < 3:
            logger.info('There are too few genes for {} skipping'.format(gene_id))
            continue

        # Check if done before
        gene_list_json = os.path.join(gene_dir, 'gene_list.json')
        #if os.path.isfile(os.path.join(gene_dir, 'parsnp.tree')) and (not overwrite):
        if not overwrite:
            if os.path.isfile(gene_list_json):
                with open(gene_list_json) as fn:
                    existing_gene_list = json.load(fn)
                    if gene_list == existing_gene_list:
                        logger.info('Phylogeny for gene {} done, skipping'.format(gene_id))
                        continue  # for _, row

        # Use the first gene as the reference
        cmd = 'parsnp -d {} -r {} -o {} -p {}'.format(
            ' '.join(gene_files[1:]), gene_files[0], gene_dir, threads)
        ret = run_command(cmd, timing_log)
        #if ret != 0:
        #    raise Exception('error')
        with open(gene_list_json, 'w') as fn:
            json.dump(gene_list, fn)

    report['alignments'] = alignment_dir
    return report
