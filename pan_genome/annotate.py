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
        {'name':"IS",'dir':None,'tool': 'blastp','MINCOV': 90,'EVALUE': 1E-30},
        {'name':"AMR",'dir':None,'tool': 'blastp','MINCOV': 90,'EVALUE': 1E-300},
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

    genus_db = []
    for genus in os.listdir(os.path.join(db_dir,'genus')):
        if re.search(r'\.', genus) != None:
            continue
        fasta = os.path.join(db_dir,'genus',genus)
        genus_db.append(fasta)
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

    return bacteria_db, genus_db, hmm_db[0]


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
        cmd = f'blastp -query - -db {db_dir} -evalue {evalue} -qcov_hsp_perc {mincov} -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no '
    elif database['tool'] == 'hmmer3':
        db_dir = database['dir']
        evalue = database['EVALUE']
        cmd = f'hmmscan --noali --notextw --acc -E {evalue} --cpu 1 {db_dir} /dev/stdin '

    parallel_cmd = f"cat {faa_file} | parallel --gnu --plain -j {threads} --block {bsize} --recstart '>' --pipe " 
    parallel_cmd += cmd
    parallel_cmd += f"> {out_file} 2> /dev/null"

    logging.info(f"Running: {parallel_cmd}")
    utils.run_command(parallel_cmd)

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

    
def annotate_cluster(unlabeled_clusters, rep_fasta, temp_dir, db_dir, threads, genus=None):
    starttime = datetime.now()
    
    out_dir = os.path.join(temp_dir, 'annotate')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    evalue = 1E-9
    mincov = 80

    # setup database
    bacteria_db, genus_db, hmm = setup_db(db_dir)
    bacteria_db[2]['MINCOV'] = mincov
    bacteria_db[2]['EVALUE'] = evalue
    hmm_db = {'name':"hmm",'dir':hmm,'tool': 'hmmer3','EVALUE': evalue}

    # order database
    ordered_database = []
    if genus != None:
        genus_fasta = os.path.join(db_dir, 'genus', genus)
        if genus_fasta in genus_db:
            db = {'name':"genus",'dir':genus_fasta,'tool': 'blastp','MINCOV': mincov,'EVALUE': evalue}
            ordered_database.append(db)
        else:
            logging.info(f'Skipping genus-specific proteins as cant see {genus_fasta}')
    
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
    clusters_name_count = {}
    suffix = 1
    for cluster in unlabeled_clusters:
        cluster_name = 'groups_' + str(suffix)
        suffix += 1
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
            clusters_name_count[cluster_name] += 1
            cluster_name += '_{}'.format(str(clusters_name_count[cluster_name]))
        else:
            clusters_name_count[cluster_name] = 0
        
        clusters_annotation.append([cluster_name, product])
        

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
    #     threads=4, 
    #     genus="Staphylococcus"
    # )

    # # dictionary = parse_search_result(
    # #     '/home/ted/test_prodigal/out/1/temp/annotate/sprot.out',
    # #     'blastp'
    # # )

    # json.dump(annotated_clusters, open('/home/ted/test_prodigal/out/1/temp/annotate/annotated_clusters.json', 'w'), indent=4, sort_keys=True)