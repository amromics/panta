import os
import logging
import multiprocessing
from datetime import datetime
from panta.utils import run_command, parse_cluster_file, chunk_fasta_file

logger = logging.getLogger(__name__)


def run_cd_hit(faa_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    cd_hit_represent_fasta = os.path.join(out_dir, 'cd-hit.fasta')
    cd_hit_cluster_file = cd_hit_represent_fasta + '.clstr'
    cmd = f'cd-hit -i {faa_file} -o {cd_hit_represent_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running cd-hit')
    cd_hit_clusters = parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity -- time taken {str(elapsed)}')
    return cd_hit_represent_fasta, cd_hit_clusters


def run_cd_hit_with_map(faa_file, map_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    cd_hit_represent_fasta = os.path.join(out_dir, 'cd-hit_tmp.fasta')
    cd_hit_cluster_file = cd_hit_represent_fasta + '.clstr'
    cmd = f'cd-hit -i {faa_file} -o {cd_hit_represent_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running cd-hit')        

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity part 1 -- time taken {elapsed}')

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    cd_hit_represent_corrected_fasta = os.path.join(out_dir, 'cd-hit.fasta')
    with open(cd_hit_represent_corrected_fasta,'w') as ofh, open(cd_hit_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                ofh.write(f'>{gene_map[line[1:].strip()]}\n')
            else:
                ofh.write(line)      

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity part 1 -- time taken {elapsed}')
    with open(cd_hit_cluster_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line[0].startswith('>'):
                cluster_name = line[1:]
                clusters[cluster_name] = {'gene_names':[]}                
            else:
                _,_, line = line.partition(', >')            
                gene_name,_,identity = line.partition('... ')
                gene_name, identity                    
                if identity == '*':
                    clusters[cluster_name]['representative'] = gene_map[gene_name]
                elif identity: # make sure it is a valid string
                    clusters[cluster_name]['gene_names'].append(gene_map[gene_name])                
    
    del gene_map    

    # convert to a simple dictionary
    clusters_new = {}
    for cluster_name in clusters:
        clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity -- time taken {str(elapsed)}')
    return cd_hit_represent_corrected_fasta, clusters_new



def run_blast(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4):
    starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    blast_db = os.path.join(out_dir, 'output_contigs')
    cmd = f"makeblastdb -in {database_fasta} -dbtype prot -out {blast_db} -logfile /dev/null"
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running makeblastdb')
    
    # chunk fasta file
    chunk_dir = os.path.join(out_dir, 'chunk_files')
    chunked_file_list = chunk_fasta_file(query_fasta, chunk_dir)

    # run parallel all-against-all blast
    #blast_cmds_file = os.path.join(out_dir,"blast_cmds.txt")    
    blast_output_file_list = []
    pool = multiprocessing.Pool(processes=threads)
    results = []      

    #with open(blast_cmds_file,'w') as fh:
    for chunked_file in chunked_file_list:
        blast_output_file = os.path.splitext(chunked_file)[0] + '.out'
        blast_output_file_list.append(blast_output_file)
        cmd = f'blastp -query {chunked_file} -db {blast_db} -evalue {evalue} -num_threads 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -max_target_seqs 2000 2> /dev/null 1> {blast_output_file}'
        results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    
    for result in results:
        if result.get() != 0:
            raise Exception('Error running all-against-all blast')        
    
    # combining blast results
    blast_result = os.path.join(out_dir, 'blast_results')
    if os.path.isfile(blast_result):
        os.remove(blast_result)
    for blast_output_file in blast_output_file_list:
        os.system(f'cat {blast_output_file} >> {blast_result}')
        os.remove(blast_output_file)

    elapsed = datetime.now() - starttime
    logging.info(f'All-against-all BLASTP -- time taken {str(elapsed)}')
    return blast_result


def run_diamond(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4):
    starttime = datetime.now()
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # make diamond database
    diamond_db = os.path.join(out_dir, 'diamond_db')
    cmd = f'diamond makedb --in {database_fasta} -d {diamond_db} -p {threads} --quiet'
    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond makedb')
    
    # run diamond blastp
    diamond_result = os.path.join(out_dir, 'diamond.tsv')
    cmd = f'diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --max-target-seqs 2000 2> /dev/null 1> {diamond_result}'
    #subprocess.call(cmd, shell=True)
    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond makedb')


    elapsed = datetime.now() - starttime
    logging.info(f'Protein alignment with Diamond -- time taken {str(elapsed)}')
    return diamond_result


def pairwise_alignment(diamond, database_fasta, query_fasta, out_dir, evalue=1E-6 ,threads=4):
    if diamond == False:
        blast_result = run_blast(
            database_fasta = database_fasta,
            query_fasta = query_fasta,
            out_dir = out_dir,
            evalue = evalue,
            threads=threads
            )
    else:
        blast_result = run_diamond(
            database_fasta = database_fasta,
            query_fasta = query_fasta,
            out_dir = out_dir,
            evalue = evalue,
            threads=threads
            )
    return blast_result


def filter_blast_result(blast_result, 
                        # gene_annotation, 
                        out_dir, identity, length_difference, alignment_coverage_short, alignment_coverage_long):
    filtered_blast_result = os.path.join(out_dir, 'filtered_blast_results')

    with open(filtered_blast_result, 'w') as fh:
        for line in open(blast_result, 'r'):
            cells = line.rstrip().split('\t')            

            qlen = int(cells[12])# * 3 + 3
            slen = int(cells[13])# * 3 + 3

            pident = float(cells[2]) / 100
            alignment_length = int(cells[3]) # * 3

            short_seq = min(qlen, slen)
            long_seq = max(qlen, slen)
            len_diff = short_seq / long_seq
            align_short = alignment_length / short_seq
            align_long = alignment_length / long_seq
            
            if pident <= identity or len_diff <= length_difference or align_short <= alignment_coverage_short or align_long <= alignment_coverage_long:
                continue

            fh.write(line)

    return filtered_blast_result

            
def cluster_with_mcl(blast_result, out_dir, threads=4):
    starttime = datetime.now()
    if threads > 1:
        threads = threads - 1
    
    mcl_file = os.path.join(out_dir, 'mcl_clusters')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_result} 2> /dev/null | mcl - --abc -I 1.5 -te {threads} -o {mcl_file} > /dev/null 2>&1"
    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mcl')
    elapsed = datetime.now() - starttime
    logging.info(f'Cluster with MCL -- time taken {str(elapsed)}')
    return mcl_file


def reinflate_clusters(cd_hit_clusters, mcl_file):
    """
    Return
    ------
        - inflated_clusters: list of list of genes
        -clusters: dict(cluster_id->[gene_ids])
    """    
    starttime = datetime.now()
    clusters = {}
    clusters.update(cd_hit_clusters)

    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                if gene in cd_hit_clusters:
                    inflated_genes.extend(cd_hit_clusters[gene])
                    del cd_hit_clusters[gene]
            inflated_clusters.append(inflated_genes)
    
    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in cd_hit_clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(cd_hit_clusters[gene])
        inflated_clusters.append(inflated_genes)
    
    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate clusters -- time taken {str(elapsed)}')
    return inflated_clusters, clusters
