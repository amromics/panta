import os
import re
import logging
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime
from glob import glob
import multiprocessing

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
        os.system(cmd)
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
        os.system(cmd)
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
        os.system(cmd)
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
    parallel_cmd += f"> {out_file} " # 2> /dev/null"

    logging.info(f"Running: {parallel_cmd}")
    os.system(parallel_cmd)

    return out_file

def annotate_cluster(rep_fasta, temp_dir, db_dir, threads, genus=None):
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
    for database in ordered_database:
        out_file = run_parallel(faa_file, threads, database, out_dir)
        










    
    # annotated_clusters = {}
    # suffix = 1
    # for cluster in unlabeled_clusters:
    #     cluster_name = 'groups_{:05d}'.format(suffix)
    #     annotated_clusters[cluster_name] = {'gene_id':cluster, 'product':'unknown'}
    #     suffix += 1

    # elapsed = datetime.now() - starttime
    # logging.info(f'Annotate clusters -- time taken {str(elapsed)}')
    # return annotated_clusters



# def create_nucleotide_representative_fasta():
#     representative_list = set()
#     for cluster in unlabeled_clusters:
#         length_max = 0
#         representative = None
#         for gene_id in cluster:
#             length = gene_annotation[gene_id][2]
#             if length > length_max:
#                 representative = gene_id
#                 length_max = length
#         representative_list.add(representative)
    
#     rep_fna = os.path.join(temp_dir, 'rep.fna')
#     with open(rep_fna, 'w') as out_fh:
#         for sample in samples:
#             sample_id = sample['id']
#             fna_file = os.path.join(collection_dir, 'samples', sample_id, sample_id + '.fna')
#             with open(fna_file, 'r') as in_fh:
#                 for line in in_fh:
#                     result = re.match(r"^>(\S+)", line)
#                     if result != None:
#                         skip = False
#                         seq_id = result.group(1)
#                         if seq_id not in representative_list:
#                             skip = True
#                             continue
#                         out_fh.write(line)
#                     else:
#                         if skip == True:
#                             continue
#                         else:
#                             out_fh.write(line)

#     # prokka
#     cmd = f"prokka --force --cpus {threads} --prefix prokka --outdir {temp_dir} --quiet --norrna --notrna --usegenus"
#     if genus != None:
#         cmd += ' --genus ' + genus
#     if species  != None:
#         cmd += ' --species ' + species
#     cmd += ' ' + rep_fna
#     ret = os.system(cmd)
#     if ret != 0:
#         raise Exception('Error running prokka')



if __name__ == "__main__":

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s : %(message)s',
        datefmt='%I:%M:%S')
    logger = logging.getLogger(__name__)


    # setup_db(db_dir="/home/ted/amromics/amromics/pan-genome/db", force=True)

    annotate_cluster(
        rep_fasta='/home/ted/test_prodigal/out/1/representative.fasta', 
        temp_dir='/home/ted/test_prodigal/out/1/temp', 
        db_dir = "/home/ted/amromics/amromics/pan-genome/db",
        threads=4, 
        genus="Staphylococcus"
    )