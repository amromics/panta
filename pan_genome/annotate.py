import os
import logging
from Bio import SearchIO
from datetime import datetime
from glob import glob
import multiprocessing

import pan_genome.utils as utils

logger = logging.getLogger(__name__)

def setup_db(baseDir, timing_log, force=False):
    """
    Prepare databases if they have not been done (makeblastdb and hmmpress).

    Parameters
    ----------
    baseDir : path
        path of panta, which contains the database
    timing_log : path
        path of time.log
    force : bool
        True - force to prepare the database again, even if it has been done.
        False - do not prepare the database again.
    
    Returns
    -------
    bacteria_db : list of dict
        information of database and its parameters
    hmm_db : path
        path to HAPMAP hmm database
    """
    db_dir = os.path.join(baseDir, 'db')
    
    # Prepare blast database
    bacteria_db = [
      {'name':"AMR",'dir':None,'tool':'blastp','MINCOV':90,'EVALUE':1E-300},
      {'name':"IS",'dir':None,'tool': 'blastp','MINCOV':90,'EVALUE':1E-30},
      {'name':"sprot",'dir':None,'tool':'blastp','MINCOV':None,'EVALUE':None}]
    for db in bacteria_db:
        name = db['name']
        fasta = os.path.join(db_dir, 'bacteria', name)
        db['dir'] = fasta
        database_file = fasta + '.pin'
        if os.path.isfile(database_file) and force==False:
            # logging.info(f'Making BLASTP database {fasta} - skipping')
            continue
        if not os.path.isfile(fasta):
            raise Exception(f'Database {fasta} does not exist')
        cmd = (f'makeblastdb -hash_index -dbtype prot -in {fasta} '
               '-logfile /dev/null')
        utils.run_command(cmd, timing_log)

    # Prepare hmm database                
    hmm_db = os.path.join(db_dir, 'hmm', 'HAMAP.hmm')
    database_file = hmm_db + '.h3i'
    if not os.path.isfile(database_file) or force==True:
        # logging.info(f'Pressing HMM database {hmm} - skipping')
        if not os.path.isfile(hmm_db):
            raise Exception(f'Database {hmm_db} does not exist')
        cmd = f'hmmpress {hmm_db} > /dev/null'
        utils.run_command(cmd, timing_log)

    return bacteria_db, hmm_db


def search_sequence(faa_file, threads, database, out_dir, timing_log):
    """
    Search sequence in database. Use BLASTP or HMMER3.

    Parameters
    ----------
    faa_file : path
        fasta file contains search sequences
    threads : int
        number of threads
    database : dict
        a data structure contain information of the database, parameters
    out_dir : path
        output directory
    timing_log : path
        path to time.log

    Returns
    -------
    out_file : path
        path of search result file.
    """
    name = database['name']
    out_file = os.path.join(out_dir, name + '.out')
    faa_bytes = os.path.getsize(faa_file)
    bsize = int(faa_bytes / threads / 2)
    if database['tool'] == 'blastp':
        db_dir = database['dir']
        evalue = database['EVALUE']
        mincov = database['MINCOV']
        cmd = (f'blastp -query {faa_file} -db {db_dir} -num_threads {threads} '
               f'-mt_mode 1 -evalue {evalue} -qcov_hsp_perc {mincov} '
               '-num_descriptions 1 -num_alignments 1 -seg no '
               f'> {out_file} 2> /dev/null')
    elif database['tool'] == 'hmmer3':
        db_dir = database['dir']
        evalue = database['EVALUE']
        cmd = (f"cat {faa_file} | parallel --gnu --plain -j {threads} "
               f"--block {bsize} --recstart '>' --pipe hmmscan --noali "
               f"--notextw --acc -E {evalue} --cpu 1 {db_dir} /dev/stdin "
               f"> {out_file} 2> /dev/null")
    utils.run_command(cmd, timing_log)
    return out_file

def parse_search_result(result_file, tool, search_result):
    """
    Parse the search result.

    Parameters
    ----------
    result_file : path
        path to search result file
    tool : {blastp, hmmer3}
        to choose the correct format of result search file
    search_result : dict
        a dictionary contain search result
        {query_id: {sseqid, gene_name, gene_product}}
    """
    if tool == 'blastp':
        fmt = 'blast-text'
    elif tool == 'hmmer3':
        fmt = 'hmmer3-text'
    
    for qresult in SearchIO.parse(result_file, fmt):
        qseqid = qresult.id
        for hit in qresult:
            sseqid = hit.id
            desc = hit.description.split('~~~')
            if desc[1] != '':
                gene_name = desc[1]
            else:
                gene_name = None
            if desc[2] != '':
                product = desc[2]
            else:
                product = None
            search_result[qseqid] = {
                'id':sseqid, 'gene':gene_name, 'product':product}

    
def annotate_cluster_fasta(unlabeled_clusters, rep_fasta, temp_dir, 
                           baseDir, timing_log, threads):
    """
    Annotate gene clusters if the input data is genome assembly.

    A representative sequence is chosen from each cluster. It will be
    searched against some known database to find annotation information 
    for gene clusters. The workflow and databases is from Prokka.
    + prepare the databases
    + search iteratively through the databases
    + parse search result and extract annotation information

    Paramters
    ---------
    unlabeled_clusters : list of list
        list of sequences IDs of each cluster
    rep_fasta : path
        a fasta file containing representative sequences of clusters
    temp_dir : path
        temporary directory
    baseDir : path
        path of panta, which contains the database 
    timing_log : path
        path of time.log
    threads : int
        number of threads

    Returns
    -------
    clusters_annotation : list
        list of annotation information of each cluster
        list of [cluster_name, cluster_product]   
    """
    starttime = datetime.now()
    
    out_dir = os.path.join(temp_dir, 'annotate')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    evalue = 1E-9
    mincov = 80

    # setup database
    bacteria_db, hmm = setup_db(baseDir, timing_log)
    bacteria_db[2]['MINCOV'] = mincov
    bacteria_db[2]['EVALUE'] = evalue
    hmm_db = {'name':"hmm",'dir':hmm,'tool': 'hmmer3','EVALUE': evalue}

    # order databases
    ordered_database = []    
    ordered_database.extend(bacteria_db)
    ordered_database.append(hmm_db)

    # search iteratively
    faa_file = rep_fasta
    search_result = {}
    logging.info(f'Total number of genes: {len(unlabeled_clusters)}')
    for database in ordered_database:
        # search
        out_file = search_sequence(
            faa_file, threads, database, out_dir, timing_log)
        # parse result
        parse_search_result(out_file, database['tool'], search_result)
        # filter out found sequences
        filter_faa = os.path.join(out_dir, database['name'] + '.filter.faa')
        utils.create_fasta_exclude(
            [faa_file], list(search_result.keys()), filter_faa)
        faa_file = filter_faa
        logging.info(f'Number of genes found: {len(search_result)}')

    # Extract cluster annotation
    clusters_annotation = []
    for cluster in unlabeled_clusters:
        cluster_name = ''
        product = 'hypothetical protein'
        for gene in cluster:
            if gene in search_result:
                if search_result[gene]['gene'] != None:
                    cluster_name = search_result[gene]['gene']
                elif search_result[gene]['id'] != None:
                    cluster_name = search_result[gene]['id']

                if search_result[gene]['product'] != None:
                    product = search_result[gene]['product']
        clusters_annotation.append([cluster_name, product])
        

    elapsed = datetime.now() - starttime
    logging.info(f'Annotate clusters -- time taken {str(elapsed)}')
    return clusters_annotation


def annotate_cluster_gff(unlabeled_clusters, gene_dictionary):
    """
    Annotate gene clusters if input data is genome annotation (GFF).
    + cluster's name is the most popular gene name in this cluster.
    + cluster's product includes all gene products in this cluster.

    Parameters
    ----------
    unlabeled_clusters : list of list
        list of sequences IDs of each cluster
    gene_dictionary : dict
        contain information of each gene of all samples
        {gene_id: (sample_id, contig, length, gene_name, gene_product)}

    Returns
    -------
    clusters_annotation : list
        list of annotation information of each cluster
        list of [cluster_name, cluster_product]        
    """
    starttime = datetime.now()

    clusters_annotation = []
    for cluster in unlabeled_clusters:
        gene_name_count = {}
        max_number = 0
        cluster_name = ''
        cluster_product = set()
        for gene_id in cluster:
            this_gene_list = gene_dictionary[gene_id]
            if this_gene_list[3] != '':
                gene_name = this_gene_list[3]
                gene_name_count[gene_name]=gene_name_count.get(gene_name, 0)+1
                if gene_name_count[gene_name] > max_number:
                    cluster_name = gene_name
                    max_number = gene_name_count[gene_name]
            if this_gene_list[4] != '':
                gene_product = this_gene_list[4]
                cluster_product.add(gene_product)
        
        cluster_product = ', '.join(cluster_product)
        clusters_annotation.append([cluster_name, cluster_product])
    
    elapsed = datetime.now() - starttime
    logging.info(f'Annotate clusters -- time taken {str(elapsed)}')
    return clusters_annotation