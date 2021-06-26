import os
import shutil
import logging
import subprocess
from datetime import datetime
from pan_genome.utils import *

logger = logging.getLogger(__name__)

# def exclude_full_clusters(cd_hit_cluster_file, faa_file, number_of_samples, excluded_cluster, greater_than):
    
#     clusters = parse_cluster_file(cd_hit_cluster_file)

#     full_cluster_gene_names =[]
#     for cluster_represent in clusters:
#         other_genes = clusters[cluster_represent]
#         this_cluster = [cluster_represent] + other_genes
#         if greater_than==True and len(this_cluster) >= number_of_samples:
#             full_cluster_gene_names.extend(this_cluster)
#             excluded_cluster[cluster_represent] = clusters[cluster_represent]
#         if greater_than==False and len(this_cluster) == number_of_samples:
#             full_cluster_gene_names.extend(this_cluster)
#             excluded_cluster[cluster_represent] = clusters[cluster_represent]
#     cluster_filtered_faa_file = faa_file + '.filtered'
#     full_cluster_gene_names = set(full_cluster_gene_names)
#     create_fasta_exclude(
#         fasta_file=faa_file, 
#         exclude_list=full_cluster_gene_names, 
#         output_file=cluster_filtered_faa_file
#         )
#     shutil.move(cluster_filtered_faa_file, faa_file)


def run_cd_hit(faa_file, out_dir, threads=4, timing_log=None):
    starttime = datetime.now()
    
    cd_hit_represent_fasta = os.path.join(out_dir, 'cd-hit.fasta')
    cd_hit_cluster_file = cd_hit_represent_fasta + '.clstr'
    cmd = f'cd-hit -i {faa_file} -o {cd_hit_represent_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit')
    cd_hit_clusters = parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity -- time taken {str(elapsed)}')
    return cd_hit_represent_fasta, cd_hit_clusters


def run_blast(database_fasta, query_fasta, out_dir, identity=95, threads=4, timing_log=None):
    starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    blast_db = os.path.join(out_dir, 'output_contigs')
    cmd = f"makeblastdb -in {database_fasta} -dbtype prot -out {blast_db} -logfile /dev/null"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running makeblastdb')
    
    # chunk fasta file
    chunk_dir = os.path.join(out_dir, 'chunk_files')
    chunked_file_list = chunk_fasta_file(query_fasta, chunk_dir)

    # run parallel all-against-all blast
    blast_cmds_file = os.path.join(out_dir,"blast_cmds.txt")
    blast_output_file_list = []
    with open(blast_cmds_file,'w') as fh:
        for chunked_file in chunked_file_list:
            blast_output_file = os.path.splitext(chunked_file)[0] + '.out'
            blast_output_file_list.append(blast_output_file)
            cmd = f"blastp -query {chunked_file} -db {blast_db} -evalue 1E-6 -num_threads 1 -outfmt 6 -max_target_seqs 2000 " 
            cmd += "| awk '{ if ($3 > " + str(identity) + ") print $0;}' 2> /dev/null 1> " + blast_output_file
            fh.write(cmd + '\n')
    cmd = f"parallel -j {threads} -a {blast_cmds_file}"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parallel all-against-all blast')

    # combining blast results
    blast_result = os.path.join(out_dir, 'blast_results')
    if os.path.isfile(blast_result):
        os.remove(blast_result)
    for blast_output_file in blast_output_file_list:
        os.system(f'cat {blast_output_file} >> {blast_result}')

    elapsed = datetime.now() - starttime
    logging.info(f'All-against-all BLASTP -- time taken {str(elapsed)}')
    return blast_result


def run_diamond(database_fasta, query_fasta, out_dir, identity=95, threads=4, timing_log=None):
    starttime = datetime.now()
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # make diamond database
    diamond_db = os.path.join(out_dir, 'diamond_db')
    cmd = f'diamond makedb --in {database_fasta} -d {diamond_db} -p {threads} --quiet'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running diamond makedb')
    
    # run diamond blastp
    diamond_result = os.path.join(out_dir, 'diamond.tsv')
    cmd = f"diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} --ultra-sensitive --outfmt 6 "
    cmd += "| awk '{ if ($3 > " + str(identity) + ") print $0;}' 2> /dev/null 1> " + diamond_result
    subprocess.call(cmd, shell=True)

    elapsed = datetime.now() - starttime
    logging.info(f'Protein alignment with Diamond -- time taken {str(elapsed)}')
    return diamond_result


def pairwise_alignment(diamond, database_fasta, query_fasta, out_dir, identity=95, threads=4, timing_log=None):
    if diamond == False:
        blast_result = run_blast(
            database_fasta = database_fasta,
            query_fasta = query_fasta,
            out_dir = out_dir,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )
    else:
        blast_result = run_diamond(
            database_fasta = database_fasta,
            query_fasta = query_fasta,
            out_dir = out_dir,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )
    return blast_result

def cluster_with_mcl(blast_result, out_dir, threads=4, timing_log=None):
    starttime = datetime.now()
    mcl_file = os.path.join(out_dir, 'mcl_clusters')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_result} 2> /dev/null | mcl - --abc -I 1.5 -o {mcl_file} > /dev/null 2>&1"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    elapsed = datetime.now() - starttime
    logging.info(f'Cluster with MCL -- time taken {str(elapsed)}')
    return mcl_file


def reinflate_clusters(cd_hit_clusters, mcl_file):
    starttime = datetime.now()
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
    return inflated_clusters


# def create_representative_fasta(cd_hit_clusters, excluded_cluster, in_fasta, outdir):
#     starttime = datetime.now()

#     representative_set = set()
#     for gene in cd_hit_clusters:
#         representative_set.add(gene)

#     for gene in excluded_cluster:
#         representative_set.add(gene)

#     representative_fasta = os.path.join(outdir, 'representative.fasta')
#     create_fasta_include(
#         fasta_file=in_fasta, 
#         include_list=representative_set, 
#         output_file=representative_fasta
#         )
#     logging.info(f'Number of representative sequences: {len(representative_set)}')
#     elapsed = datetime.now() - starttime
#     logging.info(f'Create representative fasta -- time taken {str(elapsed)}')
#     return representative_fasta