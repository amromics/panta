import os
import logging
import shutil
import subprocess
from datetime import datetime

import pan_genome.utils as utils

logger = logging.getLogger(__name__)


def run_cd_hit(faa_file, out_dir, threads, timing_log, resume):
    """
    Cluster by CD-HIT and parse the resulting cluster file.

    Parameters
    ----------
    faa_file : path
        fasta file of sequences
    out_dir : path
        directory of the output files (temp_dir)
    threads : int
        number of threads
    timing_log : path
        path of time.log
    resume : list
        A boolean inside a list
        If True, resume previous analysis

    Returns
    -------
    cd_hit_represent_fasta : path
        fasta file containing representative sequences
    cd_hit_clusters : dict
        dictionary of CD-HIT clusters
        {representative seq : [other sequences ids]}
    """
    starttime = datetime.now()
    # run CD-HIT
    cd_hit_represent_fasta = os.path.join(out_dir, 'cd-hit.fasta')
    cd_hit_cluster_file = cd_hit_represent_fasta + '.clstr'
    cmd = (f'cd-hit -i {faa_file} -o {cd_hit_represent_fasta} -s 0.98 -c 0.98 '
           f'-T {threads} -M 0 -g 1 -d 256 > /dev/null')
    
    statement = (os.path.isfile(cd_hit_represent_fasta) and 
            os.path.isfile(cd_hit_cluster_file) and 
            resume[0] == True)
    if not statement:
        utils.run_command(cmd, timing_log)
        resume[0] == False 

    # Parse cluster result
    cd_hit_clusters = utils.parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity -- time taken {str(elapsed)}')
    return cd_hit_represent_fasta, cd_hit_clusters


def run_blast(database_fasta, query_fasta, out_dir, timing_log, evalue, 
              max_target_seqs, identity, query_coverage, threads, resume):
    """
    Make blast database and then run BLASTP.

    Parameters
    ----------
    database_fasta : path
        fasta file of database
    query_fasta : path
        fasta file of query
    out_dir : path
        directory of blast output files (temp dir)
    timing_log : path
        path of time.log
    evalue : float
        evalue threshold
    max_target_seq : int
        number of max target sequences
    identity : float (0..100)
        the minimum percentage of identity
    query_coverage : float (0..100)
        the minimum percentage of the query protein that has to form 
        an alignment against the reference
    threads : int
        number of threads
    resume : list
        A boolean inside a list
        If True, resume previous analysis

    Returns
    -------
    path
        path of blast result file
    """
    starttime = datetime.now()
    blast_result = os.path.join(out_dir, 'blast_results')
    if os.path.isfile(blast_result) and resume[0] == True:
        logging.info(f'Resume - Run BLASTP')
        return blast_result
    else:
        resume[0] = False

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # make blast database
    blast_db = os.path.join(out_dir, 'blast_db')
    cmd = (f"makeblastdb -in {database_fasta} -dbtype "
           f"prot -out {blast_db} -logfile /dev/null")
    utils.run_command(cmd, timing_log)
            
    # run blast
    blast_result_temp = os.path.join(out_dir, 'blast_results.tmp')
    cmd = (f'blastp -query {query_fasta} -db {blast_db} -evalue {evalue} '
           f'-num_threads {threads} -mt_mode 1 '
           f'-outfmt 6 '
           f'-max_target_seqs {max_target_seqs} '
           f'-qcov_hsp_perc {query_coverage} '
           "| awk '{ if ($3 > " + str(identity) + ") print $0;}' "
           f'2> /dev/null 1> {blast_result_temp}')
    utils.run_command(cmd, timing_log)
    shutil.move(blast_result_temp, blast_result) # confirm finish
            
    elapsed = datetime.now() - starttime
    logging.info(f'Run BLASTP -- time taken {str(elapsed)}')
    return blast_result


def run_diamond(database_fasta, query_fasta, out_dir, timing_log, evalue, 
                max_target_seqs, identity, query_coverage, threads, resume):
    """
    Make DIAMOND database and then run DIAMOND.

    Parameters
    ----------
    database_fasta : path
        fasta file of database
    query_fasta : path
        fasta file of query
    out_dir : path
        directory of blast output files (temp dir)
    timing_log : path
        path of time.log
    evalue : float
        evalue threshold
    max_target_seq : int
        number of max target sequences
    identity : float (0..100)
        the minimum percentage of identity
    query_coverage : float (0..100)
        the minimum percentage of the query protein that has to form 
        an alignment against the reference
    threads : int
        number of threads
    resume : list
        A boolean inside a list
        If True, resume previous analysis

    Returns
    -------
    path
        path of DIAMOND result file
    """    
    starttime = datetime.now()
    diamond_result = os.path.join(out_dir, 'diamond_results')
    if os.path.isfile(diamond_result) and resume[0] == True:
        logging.info(f'Resume - Run BLASTP')
        return diamond_result
    else:
        resume[0] = False

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # make diamond database
    diamond_db = os.path.join(out_dir, 'diamond_db')
    cmd = (f'diamond makedb --in {database_fasta} -d {diamond_db} '
          f'-p {threads} --quiet')
    utils.run_command(cmd, timing_log)
            
    # run diamond blastp
    diamond_result_tmp = os.path.join(out_dir, 'diamond_results.tmp')
    cmd = (f"diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} "
           f"--evalue {evalue} --max-target-seqs {max_target_seqs} "
           f'-qcov_hsp_perc {query_coverage} '
           f'--outfmt 6 '
           "| awk '{ if ($3 > " + str(identity) + ") print $0;}' "
           f" 2> /dev/null 1> {diamond_result_tmp}")
    utils.run_command(cmd, timing_log)
    shutil.move(diamond_result_tmp, diamond_result) # confirm finish
            
    elapsed = datetime.now() - starttime
    logging.info(f'Run Diamond -- time taken {str(elapsed)}')
    return diamond_result


def pairwise_alignment(
        diamond, database_fasta, query_fasta, out_dir, timing_log, evalue, 
        max_target_seqs, identity, query_coverage, threads, resume):
    """
    Make search tool, BLAST or DIAMOND.

    Parameters
    ----------
    diamond : bool
        True - Run by DIAMOND, False - Run by BLASTP
    database_fasta : path
        fasta file of database
    query_fasta : path
        fasta file of query
    out_dir : path
        directory of blast output files (temp dir)
    timing_log : path
        path of time.log
    evalue : float
        evalue threshold
    max_target_seq : int
        number of max target sequences
    identity : float (0..100)
        the minimum percentage of identity
    query_coverage : float (0..100)
        the minimum percentage of the query protein that has to form 
        an alignment against the reference
    threads : int
        number of threads
    resume : list
        A boolean inside a list
        If True, resume previous analysis

    Returns
    -------
    path
        path of result file
    """    
    # print out the number of sequences
    database_result = subprocess.run(
        f'grep ">" {database_fasta} | wc -l', 
        capture_output=True, text=True, shell=True)
    query_result = subprocess.run(
        f'grep ">" {query_fasta} | wc -l', 
        capture_output=True, text=True, shell=True)
    logging.info(
        f'Comparing {query_result.stdout.rstrip()} sequences '
        f'with {database_result.stdout.rstrip()} sequences')

    if diamond == False:
        result = run_blast(
            database_fasta, query_fasta, out_dir, timing_log, evalue, 
            max_target_seqs, identity, query_coverage, threads, resume)
    else:
        result = run_diamond(
            database_fasta, query_fasta, out_dir, timing_log, evalue, 
            max_target_seqs, identity, query_coverage, threads, resume)
    return result

            
def cluster_with_mcl(blast_result, out_dir, timing_log, resume):
    """
    Run MCL to cluster sequences

    Parameters
    ----------
    blast_result : path
        path of filtered blast result file
    out_dir
        directory for MCL output file
    timing_log : path
        path of time.log
    resume : list
        A boolean inside a list
        If True, resume previous analysis

    Returns
    -------
    path
        path of MCL result file
    """
    starttime = datetime.now()
    mcl_file = os.path.join(out_dir, 'mcl_clusters')
    if os.path.isfile(mcl_file) and resume[0] == True:
        logging.info(f'Resume - Run MCL')
        return mcl_file
    else:
        resume[0] = False
    
    cmd = (f"mcxdeblast --m9 --score r --line-mode=abc {blast_result} "
           f"2> /dev/null | mcl - --abc -I 1.5 -o {mcl_file} > /dev/null 2>&1")
    utils.run_command(cmd, timing_log)
            
    elapsed = datetime.now() - starttime
    logging.info(f'Cluster with MCL -- time taken {str(elapsed)}')
    return mcl_file


def reinflate_clusters(cd_hit_clusters, mcl_file):
    """
    Bring back sequences that are excluded in CD-HIT step.

    Parameters
    ----------
    cd_hit_clusters : dict
        dictionary of CD-HIT clusters
        {representative seq : [other sequences ids]}
    mcl_file : path
        path of MCL result file

    Returns
    -------
    list 
        list of complete clusters

    """
    starttime = datetime.now()

    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for i,line in enumerate(fh):
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                if gene in cd_hit_clusters:
                    inflated_genes.extend(cd_hit_clusters[gene])
                    del cd_hit_clusters[gene]
            inflated_clusters.append(inflated_genes)

    # Inflate any cd-hit clusters that were not sent to mcl
    for gene in cd_hit_clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(cd_hit_clusters[gene])
        inflated_clusters.append(inflated_genes)

    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate clusters -- time taken {str(elapsed)}')
    return inflated_clusters
