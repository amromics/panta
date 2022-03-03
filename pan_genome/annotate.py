import os
import re
import logging
from Bio import SearchIO
from datetime import datetime
from glob import glob
import multiprocessing
import pan_genome.utils as utils

logger = logging.getLogger(__name__)

def setup_db(db_dir, force=False):

    bacteria_db = [
        {'name':"AMR",'dir':None,'tool': 'blastp','MINCOV': 90,'EVALUE': 1E-300},
        {'name':"IS",'dir':None,'tool': 'blastp','MINCOV': 90,'EVALUE': 1E-30},
        {'name':"sprot",'dir':None,'tool': 'blastp','MINCOV': None,'EVALUE': None}
    ]
    
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
        cmd = f'makeblastdb -hash_index -dbtype prot -in {fasta} -logfile /dev/null'
        utils.run_command(cmd)
        logging.info(f'Making BLASTP database {fasta}')


    hmm_db =glob(os.path.join(db_dir, 'hmm', '*.hmm'))
    for hmm in hmm_db:
        database_file = hmm + '.h3i'
        if os.path.isfile(database_file) and force==False:
            # logging.info(f'Pressing HMM database {hmm} - skipping')
            continue
        if not os.path.isfile(hmm):
            raise Exception(f'Database {hmm} does not exist')
        cmd = f'hmmpress {hmm} > /dev/null'
        utils.run_command(cmd)
        logging.info(f'Pressing HMM database {hmm}')

    return bacteria_db, hmm_db[0]


def run_parallel(faa_file, threads, database, out_dir):
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
        cmd = f'blastp -query {faa_file} -db {db_dir} -num_threads {threads} -mt_mode 1 -evalue {evalue} -qcov_hsp_perc {mincov} -num_descriptions 1 -num_alignments 1 -seg no > {out_file} 2> /dev/null'
    elif database['tool'] == 'hmmer3':
        db_dir = database['dir']
        evalue = database['EVALUE']
        cmd = f"cat {faa_file} | parallel --gnu --plain -j {threads} --block {bsize} --recstart '>' --pipe " 
        cmd += f'hmmscan --noali --notextw --acc -E {evalue} --cpu 1 {db_dir} /dev/stdin '
        cmd += f"> {out_file} 2> /dev/null"

    logging.info(f"Running: {cmd}")
    utils.run_command(cmd)

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
            dictionary[qseqid] = {'id':sseqid, 'gene':gene_name, 'product':product}

    
def annotate_cluster_fasta(unlabeled_clusters, rep_fasta, temp_dir, db_dir, threads, start=1):
    starttime = datetime.now()
    
    out_dir = os.path.join(temp_dir, 'annotate')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    evalue = 1E-9
    mincov = 80

    # setup database
    bacteria_db, hmm = setup_db(db_dir)
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
        out_file = run_parallel(faa_file, threads, database, out_dir)
        parse_search_result(out_file, database['tool'], search_result)
        filter_faa = os.path.join(out_dir, database['name'] + '.filter.faa')
        utils.create_fasta_exclude([faa_file], list(search_result.keys()), filter_faa)
        faa_file = filter_faa
        logging.info(f'Number of genes found: {len(search_result)}')

    
    clusters_annotation = []
    clusters_name_count = []

    for i,cluster in enumerate(unlabeled_clusters, start):
        cluster_name = 'groups_' + str(i)
        product = 'hypothetical protein'
        for gene in cluster:
            if gene in search_result:
                if search_result[gene]['gene'] != None:
                    cluster_name = search_result[gene]['gene']
                elif search_result[gene]['id'] != None:
                    cluster_name = search_result[gene]['id']

                if search_result[gene]['product'] != None:
                    product = search_result[gene]['product']

        # check if cluster_name already exists
        if cluster_name in clusters_name_count:
            cluster_name += '_{}'.format(str(i))
        else:
            clusters_name_count.append(cluster_name)
        
        clusters_annotation.append([cluster_name, product])
        

    elapsed = datetime.now() - starttime
    logging.info(f'Annotate clusters -- time taken {str(elapsed)}')
    return clusters_annotation

def annotate_cluster_gff(unlabeled_clusters, gene_dictionary,start=1):
    starttime = datetime.now()

    clusters_annotation = []
    clusters_name_count = []

    for i, gene_id_list in enumerate(unlabeled_clusters,start):
        cluster_name = 'groups_' + str(i)
        cluster_product = None
        gene_name_count = {}
        max_number = 0
        for gene_id in gene_id_list:
            this_gene = gene_dictionary[gene_id]
            if this_gene[3] != '':
                gene_name = this_gene[3]
                gene_name_count[gene_name] = gene_name_count.get(gene_name, 0) + 1
                if gene_name_count[gene_name] > max_number:
                    cluster_name = gene_name
                    max_number = gene_name_count[gene_name]
        cluster_product = []
        for gene_id in gene_id_list:
            this_gene = gene_dictionary[gene_id]
            if this_gene[4] != '':
                gene_product = this_gene[4]
                if gene_product not in cluster_product:
                    cluster_product.append(gene_product)
        if len(cluster_product) > 0:
            cluster_product = ', '.join(cluster_product)
        else:
            cluster_product = ''
        
        # # check if cluster_name already exists
        # if cluster_name in clusters_name_count:
        #     cluster_name += '_{}'.format(str(i))
        # else:
        #     clusters_name_count.append(cluster_name)
        
        clusters_annotation.append([cluster_name, cluster_product])
    
    elapsed = datetime.now() - starttime
    logging.info(f'Annotate clusters -- time taken {str(elapsed)}')
    return clusters_annotation


# if __name__ == "__main__":

#     logging.basicConfig(
#         level=logging.DEBUG,
#         format='%(asctime)s %(levelname)s : %(message)s',
#         datefmt='%I:%M:%S')
#     logger = logging.getLogger(__name__)


    # setup_db(db_dir="/home/ted/amromics/amromics/pan-genome/db", force=True)

    # old_clusters = json.load(open('/home/ted/test_prodigal/out/1/clusters.json', 'r'))

    # annotated_clusters=annotate_cluster(
    #     unlabeled_clusters= old_clusters,
    #     rep_fasta='/home/ted/test_prodigal/out/1/representative.fasta', 
    #     temp_dir='/home/ted/test_prodigal/out/1/temp', 
    #     db_dir = "/home/ted/amromics/amromics/pan-genome/db",
    #     threads=4
    # )

    # # dictionary = parse_search_result(
    # #     '/home/ted/test_prodigal/out/1/temp/annotate/sprot.out',
    # #     'blastp'
    # # )

    # json.dump(annotated_clusters, open('/home/ted/test_prodigal/out/1/temp/annotate/annotated_clusters.json', 'w'), indent=4, sort_keys=True)