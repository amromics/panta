import os
import logging
import subprocess
from datetime import datetime

import pan_genome.utils as utils

logger = logging.getLogger(__name__)


def run_cd_hit(faa_file, out_dir, threads, timing_log):
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
    utils.run_command(cmd, timing_log)

    # Parse cluster result
    cd_hit_clusters = utils.parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity -- time taken {str(elapsed)}')
    return cd_hit_represent_fasta, cd_hit_clusters


def run_blast(database_fasta, query_fasta, out_dir, timing_log, 
              evalue=1E-6, max_target_seqs=2000, threads=4):
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
    threads : int
        number of threads
    
    Returns
    -------
    path
        path of blast result file
    """
    starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    blast_db = os.path.join(out_dir, 'blast_db')
    cmd = (f"makeblastdb -in {database_fasta} -dbtype "
           f"prot -out {blast_db} -logfile /dev/null")
    utils.run_command(cmd, timing_log)
            
    
    # run blast
    blast_result = os.path.join(out_dir, 'blast_results')
    cmd = (f'blastp -query {query_fasta} -db {blast_db} -evalue {evalue} '
           f'-num_threads {threads} -mt_mode 1 '
           f'-max_target_seqs {max_target_seqs} '
           f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart '
            'qend sstart send evalue bitscore qlen slen" '
           f'2> /dev/null 1> {blast_result}')
    utils.run_command(cmd, timing_log)
            
    elapsed = datetime.now() - starttime
    logging.info(f'Run BLASTP -- time taken {str(elapsed)}')
    return blast_result


def run_diamond(database_fasta, query_fasta, out_dir, timing_log, 
                evalue=1E-6, max_target_seqs=2000, threads=4):
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
    threads : int
        number of threads
    
    Returns
    -------
    path
        path of DIAMOND result file
    """    
    starttime = datetime.now()
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # make diamond database
    diamond_db = os.path.join(out_dir, 'diamond_db')
    cmd = (f'diamond makedb --in {database_fasta} -d {diamond_db} '
          f'-p {threads} --quiet')
    utils.run_command(cmd, timing_log)
            
    
    # run diamond blastp
    diamond_result = os.path.join(out_dir, 'diamond_results')
    cmd = (f"diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} "
           f"--evalue {evalue} --max-target-seqs {max_target_seqs} "
            '--outfmt "6 qseqid sseqid pident length mismatch gapopen qstart '
            'qend sstart send evalue bitscore qlen slen"'
           f" 2> /dev/null 1> {diamond_result}")
    utils.run_command(cmd, timing_log)
            

    elapsed = datetime.now() - starttime
    logging.info(f'Run Diamond -- time taken {str(elapsed)}')
    return diamond_result


def pairwise_alignment(
        diamond, database_fasta, query_fasta, out_dir, timing_log, 
        evalue=1E-6, max_target_seqs=2000, threads=4):
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
    threads : int
        number of threads
    
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
        blast_result = run_blast(
            database_fasta = database_fasta,
            query_fasta = query_fasta,
            out_dir = out_dir,
            timing_log=timing_log,
            evalue = evalue,
            max_target_seqs= max_target_seqs,
            threads=threads
            )
    else:
        blast_result = run_diamond(
            database_fasta = database_fasta,
            query_fasta = query_fasta,
            out_dir = out_dir,
            timing_log=timing_log,
            evalue = evalue,
            max_target_seqs= max_target_seqs,
            threads=threads
            )
    return blast_result


def filter_blast_result(blast_result, out_dir, args):
    """
    Filter BLAST result to find matched sequences. There are 4 
    criteria: identity, LD, AL, AS. Output the filtered blast file.

    Parameters
    ----------
    blast_result : path
        blast 1 result file
    out_dir : path
        directory of filtered fasta output
    args : object
        Command-line input arguments
    
    Returns
    -------
    filtered_blast_result : path
        BLAST result file after filtering
    """    
    filtered_blast_result = os.path.join(out_dir, 'filtered_blast_results')

    with open(filtered_blast_result, 'w') as fh:
        for line in open(blast_result, 'r'):
            cells = line.rstrip().split('\t')
            qlen = int(cells[12])
            slen = int(cells[13])
            pident = float(cells[2]) / 100
            alignment_length = int(cells[3])

            short_seq = min(qlen, slen)
            long_seq = max(qlen, slen)
            len_diff = short_seq / long_seq
            align_short = alignment_length / short_seq
            align_long = alignment_length / long_seq
            
            if (pident <= args.identity 
                    or len_diff <= args.LD 
                    or align_short <= args.AS 
                    or align_long <= args.AL):
                continue

            fh.write(line)

    return filtered_blast_result

            
def cluster_with_mcl(blast_result, out_dir, timing_log):
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
    
    Returns
    path
        path of MCL result file
    -------
    """
    starttime = datetime.now()
    mcl_file = os.path.join(out_dir, 'mcl_clusters')
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
