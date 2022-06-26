import os
import logging
from Bio import SearchIO
from datetime import datetime
from glob import glob
import multiprocessing

import pan_genome.utils as utils

logger = logging.getLogger(__name__)

def setup_db(baseDir, timing_log, force=False):
    db_dir = os.path.join(baseDir, 'db')
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
                    


    hmm_db =glob(os.path.join(db_dir, 'hmm', '*.hmm'))
    for hmm in hmm_db:
        database_file = hmm + '.h3i'
        if os.path.isfile(database_file) and force==False:
            # logging.info(f'Pressing HMM database {hmm} - skipping')
            continue
        if not os.path.isfile(hmm):
            raise Exception(f'Database {hmm} does not exist')
        cmd = f'hmmpress {hmm} > /dev/null'
        utils.run_command(cmd, timing_log)
                    

    return bacteria_db, hmm_db[0]


def run_parallel(faa_file, threads, database, out_dir, timing_log):
    name = database['name']
    out_file = os.path.join(out_dir, name + '.out')
    faa_bytes = os.path.getsize(faa_file)
    if threads > 0:
        bsize = int(faa_bytes / threads / 2)
    else:
        ncpu = multiprocessing.cpu_count()
        bsize = int(faa_bytes / ncpu / 2)
    
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

def parse_search_result(result_file, tool, dictionary):
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
            dictionary[qseqid] = {
                'id':sseqid, 'gene':gene_name, 'product':product}

    
def annotate_cluster_fasta(unlabeled_clusters, rep_fasta, temp_dir, 
                           baseDir, timing_log, threads):
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

    # order database
    ordered_database = []    
    ordered_database.extend(bacteria_db)
    ordered_database.append(hmm_db)

    # search iteratively
    faa_file = rep_fasta
    search_result = {}
    logging.info(f'Total number of genes: {len(unlabeled_clusters)}')
    for database in ordered_database:
        out_file = run_parallel(
            faa_file, threads, database, out_dir, timing_log)
        parse_search_result(out_file, database['tool'], search_result)
        filter_faa = os.path.join(out_dir, database['name'] + '.filter.faa')
        utils.create_fasta_exclude(
            [faa_file], list(search_result.keys()), filter_faa)
        faa_file = filter_faa
        logging.info(f'Number of genes found: {len(search_result)}')

    
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