import os,shutil
import logging
import multiprocessing
from Bio import SeqIO
import csv
import pandas as pd
import json
import gzip
import re
import amromics.libs.bioseq as bioseq
from amromics.utils.command import run_command
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_phylogeny_parsnp(base_dir,ref_genome, genome_dir='.',threads=0):
    """
        Run parsnp to create phylogeny tree
        :param read_data: result holder
        :param ref_genome: path to reference genome, if equal None, one of genome in genome directory will be choosed to be reference.
        :param base_dir: working directory
        :param threads: number of core CPU
        :return:
    """
    phylogeny_folder=os.path.join(base_dir,'pangenome/phylogeny')
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    else:

        return phylogeny_folder
    temp_folder = os.path.join(phylogeny_folder, 'temp_phylo')
    if not os.path.isdir(temp_folder):
        os.makedirs(temp_folder)
    #take first genome to get reference genome
    sample_list = []
    files=os.listdir(genome_dir)
    if ref_genome==None:
        #pick a file in genome_dir to make ref
        for i,f in enumerate(files):
            fasta_file = os.path.join(temp_folder, os.path.basename(f))
            cmd = 'gunzip -c {} > {}'.format(os.path.join(genome_dir,f), fasta_file)
            run_command(cmd)
            cmd = 'cat {} > {}'.format( os.path.join(genome_dir,f), fasta_file)
            run_command(cmd)
            if i == 0:
                ref_genome = fasta_file
            else:
                sample_list.append(fasta_file)
    else:
        for i,f in enumerate(files):
            fasta_file = os.path.join(temp_folder, os.path.basename(f))
            cmd = 'gunzip -c {} > {}'.format(os.path.join(genome_dir,f), fasta_file)
            run_command(cmd)
            cmd = 'zcat {} > {}'.format(f, fasta_file)
            run_command(cmd)
            sample_list.append(fasta_file)
    myCmd = 'parsnp -r {} -d {} -o {} -p {}'.format(
        ref_genome,
        ' '.join(sample_list),
        phylogeny_folder, threads)
    print(myCmd)
    et = run_command(myCmd, timing_log)
    if ret != 0:
        raise Exception('Error running parsnp')
    run_command('gzip {}'.format(os.path.join(phylogeny_folder, 'parsnp.xmfa')))
    run_command('gzip {}'.format(os.path.join(phylogeny_folder, 'parsnp.ggr')))
    shutil.rmtree(temp_folder)

    return phylogeny_folder
def run_species_phylogeny_iqtree(roary_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Run iqtree to create phylogeny tree from core gene alignment. If the list of samples has
    not changed, and none of the samples has changed, the existing tree will be kept unless
    overwrite is set to True
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
    #report['phylogeny'] = phylogeny_folder

    phylogeny_file = os.path.join(phylogeny_folder, 'core_gene_alignment.treefile')
    if os.path.isfile(phylogeny_file) and (not overwrite):
        logger.info('phylogeny tree exists and input has not changed, skip phylogeny analysis')
        return phylogeny_folder

    aln_file = os.path.join(phylogeny_folder, 'core_gene_alignment.aln.gz')
    if not os.path.isfile(aln_file):
        aln_file = os.path.join(report['roary'], 'core_gene_alignment.aln.gz')
    cmd = 'iqtree -s {alignment} --prefix {prefix} -B 1000 -T {threads} -czb -keep-ident'.format(
        alignment=aln_file, prefix=phylogeny_folder+'/core_gene_alignment', threads=threads)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('iqtree fail to create phylogeny tree from core gene alignment!')

    return phylogeny_folder


def run_gene_phylogeny_iqtree(roary_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
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
    alignment_dir = os.path.join(collection_dir, 'alignments')
    gene_cluster_file = roary_folder + '/gene_presence_absence.Rtab'
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    cmds_file = os.path.join(alignment_dir,"phylo_cmds")
    with open(cmds_file,'w') as cmds:
        for gene_id, row in gene_df.iterrows():
            # Only analyse if there are at least 3 genes
            if row.sum() < 3:
                continue

            gene_id = re.sub(r'\W+', '', gene_id)
            gene_dir = os.path.join(alignment_dir, gene_id)
            if not os.path.exists(gene_dir):
                os.makedirs(gene_dir)
            # check if done before
            iqtree_output = os.path.join(gene_dir, gene_id + '.treefile')
            if (not overwrite) and os.path.isfile(iqtree_output):
                continue

            gene_aln_file_roary = os.path.join(roary_folder,'pan_genome_sequences', gene_id + '.fa.aln')
            gene_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln')
            if os.path.isfile(gene_aln_file_roary):
                shutil.move(gene_aln_file_roary,gene_aln_file)
            if not os.path.isfile(gene_aln_file):
                logger.info('{} does not exist'.format(gene_aln_file))
                continue

            cmd = f"iqtree -s {gene_aln_file} --prefix {gene_dir+'/'+gene_id} -m GTR -quiet -T 1 -B 1000 2> /dev/null"
            cmd += f" || iqtree -s {gene_aln_file} --prefix {gene_dir+'/'+gene_id} -m GTR -quiet -T 1"
            # translate to protein alignment
            #protein_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln')
            #with open(protein_aln_file, 'w') as fh:
            #    for record in SeqIO.parse(gene_aln_file, 'fasta'):
            #        trans = translate_dna(str(record.seq))
            #        new_record = SeqRecord(Seq(trans), id=record.id,)
            #        SeqIO.write(new_record, fh, 'fasta')
            #cmd = f"iqtree -s {protein_aln_file} --prefix {gene_dir+'/'+gene_id} -m LG -quiet -T 1"
            #cmd = f"fasttree -nt -gtr -quiet {gene_aln_file} > {gene_dir+'/'+gene_id+'.treefile'} && echo '{gen_list_string}' > {gene_list_json}"
            cmds.write(cmd + '\n')

    cmd = f"parallel --bar -j {threads} -a {cmds_file}"
    ret = run_command(cmd, timing_log)



    return alignment_dir
