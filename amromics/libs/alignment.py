import os,shutil
import logging
import multiprocessing
import gzip
from Bio import SeqIO
import csv
import pandas as pd
import json
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import amromics.libs.bioseq as bioseq
from amromics.utils.command import run_command
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_alignment_by_parsnp(roary_folder,ffn_dir,base_dir, overwrite=False,  timing_log=None,threads=0):
    """
        Run aligment process to create both multi-alignment  and phylogeny tree for each gene in gene clusters
        :param read_data: result holder
        :param ffn_dir: path to folder of .ffn (output of prokka)
        :param base_dir: working directory
        :return:
    """
    gene_cluster_file=roary_folder+'/gene_presence_absence.csv.gz'
    dict_cds={}
    for root, dirs, files in os.walk(ffn_dir):
        for _file in files:
            if _file.endswith('.ffn.gz'):
                with gzip.open(os.path.join(root, _file), 'rt') as fn:
                    for seq in SeqIO.parse(fn, 'fasta'):
                        dict_cds[seq.id] = seq


    #make folder contains sequences for each gene
    alignment_dir=os.path.join(base_dir,'alignments')
    if (not overwrite) and os.path.exists(alignment_dir):
        return alignment_dir
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)

    gene_df = pd.read_csv(gene_cluster_file, dtype=str)
    gene_df.fillna('', inplace=True)

    sample_columns = list(gene_df.columns)[14:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                # roary can pool together genes from the same sample and tab-separate them
                for sample_gene in row[sample_column].split('\t'):
                    gene_list.append(sample_gene)
                    # TODO: make sure all samples in this gene have not updated

        gene_list = sorted(gene_list)
        # Only analyse if there are more than 3 genes
        if len(gene_list) < 3:
            logger.info('There are too few genes for {} skipping'.format(gene_id))
            continue

        gene_dir = os.path.join(alignment_dir, gene_id)
        # Check if done before
        gene_list_json = os.path.join(gene_dir, 'gene_list.json')
        # if os.path.isfile(os.path.join(gene_dir, 'parsnp.tree')) and (not overwrite):
        if os.path.isfile(gene_list_json):
            with open(gene_list_json) as fn:
                existing_gene_list = json.load(fn)
                if gene_list == existing_gene_list:
                    logger.info('Phylogeny for gene {} done, skipping'.format(gene_id))
                    continue  # for _, row

        gene_file_dir = os.path.join(gene_dir, 'files')
        if not os.path.exists(gene_file_dir):
            os.makedirs(gene_file_dir)

        gene_files = []
        for sample_gene in gene_list:
            gene_file = os.path.join(gene_file_dir, sample_gene + '.fasta')
            SeqIO.write(dict_cds[sample_gene], gene_file, 'fasta')
            gene_files.append(gene_file)

        # Use the first gene as the reference
        cmd = 'parsnp -d {} -r {} -o {} -p {}'.format(
            ' '.join(gene_files[1:]), gene_files[0], gene_dir, threads)
        ret = run_command(cmd)
        # if ret != 0:
        #     raise Exception('error')

        with open(gene_list_json, 'w') as fn:
            json.dump(gene_list, fn)
        #run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.xmfa')))
        #run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.ggr')))

        if os.path.exists(gene_file_dir):
            shutil.rmtree(gene_file_dir)
        #clean up
        run_command('rm -f ' + os.path.join(gene_dir, '*.ini ') + os.path.join(gene_dir, '*block* '))
        shutil.rmtree(os.path.join(gene_dir, 'blocks'), True)
        shutil.rmtree(os.path.join(gene_dir, 'tmp'), True)


    return alignment_dir
def run_protein_alignment(roary_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Align protein sequence by mafft

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
    alignment_dir = os.path.join(collection_dir, 'alignments')

    gene_cluster_file =roary_folder + '/gene_presence_absence.Rtab'
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    cmds_file = os.path.join(alignment_dir,"align_cmds")
    with open(cmds_file,'w') as cmds:
        for gene_id, row in gene_df.iterrows():
            # Only align if there are at least 2 sequences
            if row.sum() < 2:
                continue

            gene_id = re.sub(r'\W+', '', gene_id)
            gene_dir = os.path.join(alignment_dir, gene_id)

            # check if done before
            gene_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln.gz')
            if (not overwrite) and os.path.isfile(gene_aln_file):
                continue

            gene_seq_file = os.path.join(gene_dir, gene_id + '.faa')
            if not os.path.isfile(gene_seq_file):
                logger.info('{} does not exist'.format(gene_aln_file))
                continue

            cmd = f"mafft --auto --quiet --thread 1 {gene_seq_file} | gzip > {gene_aln_file}"
            cmds.write(cmd + '\n')

    cmd = f"parallel --bar -j {threads} -a {cmds_file}"
    ret = run_command(cmd, timing_log)
    #report['alignments'] = alignment_dir
    return alignment_dir


def create_nucleotide_alignment(roary_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Create nucleotide alignment according to protein alignment

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
    logger.info('Creating nucleotide alignment')
    alignment_dir = os.path.join(collection_dir, 'alignments')

    gene_cluster_file = roary_folder + '/gene_presence_absence.Rtab'
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    for gene_id, row in gene_df.iterrows():
        # Only run if there are at least 2 sequences
        if row.sum() < 2:
            continue
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)

        # check if done before
        nucleotide_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln.gz')
        if (not overwrite) and os.path.isfile(nucleotide_aln_file):
            continue

        protein_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln.gz')
        if not os.path.isfile(protein_aln_file):
            logger.info('{} does not exist'.format(protein_aln_file))
            continue
        protein_dict = {}
        with gzip.open(protein_aln_file, 'rt') as fh:
            for seq_record in SeqIO.parse(fh, 'fasta'):
                protein_dict[seq_record.id] = str(seq_record.seq)

        nucleotide_seq_file = os.path.join(gene_dir, gene_id + '.fna')
        nucleotide_dict = {}
        for seq_record in SeqIO.parse(nucleotide_seq_file, 'fasta'):
            nucleotide_dict[seq_record.id] = str(seq_record.seq)

        with gzip.open(nucleotide_aln_file, 'wt') as fh:
            for seq_id in protein_dict.keys():
                protein = protein_dict[seq_id]
                nucleotide = nucleotide_dict[seq_id]
                result = ''
                codon_pos = 0
                for c in protein:
                    if c == '-':
                        result += '---'
                    else:
                        result += nucleotide[codon_pos * 3: codon_pos * 3 + 3]
                        codon_pos += 1
                new_record = SeqRecord(Seq(result), id = seq_id, description = '')
                SeqIO.write(new_record, fh, 'fasta')
        
        # remove input files
        os.remove(nucleotide_seq_file)
        protein_seq_file = os.path.join(gene_dir, gene_id + '.faa')
        os.remove(protein_seq_file)

    return alignment_dir


def create_core_gene_alignment(roary_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Concatenate all the nucleotide alignment of core genes to create core gene alignment

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
    logger.info('Creating core gene alignment')
    alignment_dir = os.path.join(collection_dir, 'alignments')
    phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
    #report['phylogeny'] = phylogeny_folder
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    # check if done before
    core_gene_aln_file = os.path.join(phylogeny_folder, 'core_gene_alignment.aln.gz')
    if os.path.isfile(core_gene_aln_file) and (not overwrite):
        logger.info('Core gene alignment exists and input has not changed, skipping')
        return phylogeny_folder

    gene_cluster_file = roary_folder + '/gene_presence_absence.Rtab'
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    seq_dict = {}
    samples=gene_df.columns.tolist()
    for sample in samples:
        seq_dict[sample]= ''
    sample_list = seq_dict.keys()
    print(sample_list)
    for gene_id, row in gene_df.iterrows():
        # Only run if it is core gene
        if len(row[row == 0]) != 0:
            continue
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)

        nucleotide_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln.gz')
        if not os.path.isfile(nucleotide_aln_file):
            logger.info('{} does not exist'.format(nucleotide_aln_file))
            continue
        cluster_dict = {}
        with gzip.open(nucleotide_aln_file, 'rt') as fh:
            for seq_record in SeqIO.parse(fh, 'fasta'):
                sample_name = re.findall(r'^(.+)_', seq_record.id)
                sample_name = sample_name[0]
                if sample_name not in sample_list:
                    raise Exception(f'Error concatenating gene alignment: {sample_name} is not a sample id')
                cluster_dict[sample_name] = str(seq_record.seq)
        
        for sample_name in cluster_dict:
            seq_dict[sample_name] += cluster_dict[sample_name]

    with gzip.open(core_gene_aln_file, 'wt') as fh:
        for sample in sample_list:
            new_record = SeqRecord(Seq(seq_dict[sample]), id = sample, description = '')
            SeqIO.write(new_record, fh, 'fasta')

    return phylogeny_folder
def get_gene_sequences(roary_folder,ffn_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Create protein sequences and nucleotide sequences for each gene cluster

    Parameters
    ----------
    roary_folder: string
        roary_folder
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
    logger.info('Getting sequences of gene clusters')
    gene_cluster_file = roary_folder + '/gene_presence_absence.csv.gz'
    dict_nucleotide = {}
    for ffn_file in os.listdir(ffn_folder):
        ffn_file_path=os.path.join(ffn_folder,ffn_file)
        if ffn_file.endswith('.gz'):
            with gzip.open(ffn_file_path, 'rt') as fn:
                for seq_record in SeqIO.parse(fn, 'fasta'):
                    seq_record.seq = seq_record.seq[:-3]
                    seq_record = SeqRecord(seq_record.seq, id=seq_record.id, description = '')
                    dict_nucleotide[seq_record.id] = seq_record
        else:
            for seq_record in SeqIO.parse(ffn_file_path, 'fasta'):
                seq_record.seq = seq_record.seq[:-3]
                seq_record = SeqRecord(seq_record.seq, id=seq_record.id, description = '')
                dict_nucleotide[seq_record.id] = seq_record

    # make folder contains sequences for each gene
    alignment_dir = os.path.join(collection_dir, 'alignments')
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    #report['alignments'] = alignment_dir

    gene_df = pd.read_csv(gene_cluster_file, dtype=str, compression='gzip')
    gene_df.fillna('', inplace=True)
    sample_columns = list(gene_df.columns)[14:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)

        # check if done before
        protein_seq_file = os.path.join(gene_dir, gene_id + '.faa')
        nucleotide_seq_file = os.path.join(gene_dir, gene_id + '.fna')
        if (not overwrite) and os.path.isfile(protein_seq_file) and os.path.isfile(nucleotide_seq_file):
            continue

        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                # roary can pool together genes from the same sample and tab-separate them
                for sample_gene in row[sample_column].split('\t'):
                    gene_list.append(sample_gene)
        gene_list = sorted(gene_list)

        with open(protein_seq_file, 'w') as prot_fh, open(nucleotide_seq_file, 'w') as nucl_fh:
            for sample_gene in gene_list:
                nu_seq_record = dict_nucleotide[sample_gene]
                SeqIO.write(nu_seq_record, nucl_fh, 'fasta')
                pro_seq = nu_seq_record.seq.translate(table=11)
                pro_seq_record = SeqRecord(pro_seq, id = nu_seq_record.id, description = '')
                SeqIO.write(pro_seq_record, prot_fh, 'fasta')

    return alignment_dir
def runGeneAlignment(roary_folder,ffn_dir, collection_dir, threads=8, overwrite=False, timing_log=None):
    alignment_dir=get_gene_sequences(roary_folder, ffn_dir,overwrite=overwrite,collection_dir=collection_dir, threads=threads,timing_log=timing_log)
    alignment_dir=run_protein_alignment(roary_folder, collection_dir=collection_dir, overwrite=overwrite,threads=threads,timing_log=timing_log)
    alignment_dir=create_nucleotide_alignment(roary_folder, collection_dir=collection_dir, overwrite=overwrite,threads=threads,timing_log=timing_log)
    return alignment_dir
