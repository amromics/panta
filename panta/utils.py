import os
import psutil
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
import shutil
import pickle
import json
import numpy as np
logger = logging.getLogger(__name__)


def run_command(cmd, timing_log=None):
    """
    Run a command line, return the returning code of the command
    :param cmd:
    :param timing_log:
    :return:
    """
    if timing_log is not None:
        cmd = '/usr/bin/time --append -v -o {} bash -c "{}"'.format(timing_log, cmd)
    logger.info('Running "{}'.format(cmd))
    ret = os.system(cmd)
    return ret

def mem_report(value, point='POINT'):
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    mem_usage = mem_info.rss/1000000
    logger.info(f'MEM at {point}: {mem_usage} {value} inc = {mem_usage-value}')
    return mem_usage

def get_seq_ids(gene_id):
    toks = gene_id.split('-',2)
    #sample_id, contig_id, gene_id
    if len(toks) < 2:
        logger.error(f'See {gene_id}')
    return toks[0], toks[1]

def parse_cluster_file(cd_hit_cluster_file): 
    """
    Parse cdhit cluster file
    Return:
        a dictionary of clusters: dict (cluster_name > [gene_id])
    """    
    clusters = {}
    with open(cd_hit_cluster_file, 'r') as fh:
        for line in fh:
            result = re.match(r"^>(.+)$", line)
            if result != None:
                cluster_name = result.group(1)
                clusters[cluster_name] = {}
                clusters[cluster_name]['gene_names'] = []
            else:
                result = re.match(r'[\d]+\t[\w]+, >(.+)\.\.\. (.+)$', line)
                if result != None:
                    gene_name = result.group(1)
                    identity = result.group(2)
                    if identity == '*':
                        clusters[cluster_name]['representative'] = gene_name
                    else:
                        percent = re.findall(r'([0-9\.]+)', identity)
                        percent = float(percent[0])
                        clusters[cluster_name]['gene_names'].append(gene_name)    
    # convert to a simple dictionary
    clusters_new = {}
    for cluster_name in clusters:
        clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']
    return clusters_new


def parse_cluster_file_with_map(cd_hit_cluster_file, map_file=None): 
    """
    Parse cdhit cluster file
    Return:
        a dictionary of clusters: dict (cluster_name > [gene_id])
    """

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        with line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1

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
                    #percent = float(identity[3:-1])
                    #percent = re.findall(r'([0-9\.]+)', identity)
                    #percent = float(percent[0])
                    clusters[cluster_name]['gene_names'].append(gene_map[gene_name])
    
    del gene_map    
    # convert to a simple dictionary
    clusters_new = {}
    for cluster_name in clusters:
        clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']
    return clusters_new



def chunk_fasta_file(fasta_file, out_dir):
    # starttime = datetime.now()
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)
    else:
        os.makedirs(out_dir)
    
    chunked_file_list = []
    chunk_number = 0
    current_chunk_length = 0
    chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
    chunked_fh = open(chunked_file, 'w')
    chunked_file_list.append(chunked_file)
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if current_chunk_length > 200000:
            chunked_fh.close()
            chunk_number += 1
            current_chunk_length = 0
            chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
            chunked_file_list.append(chunked_file)
            chunked_fh = open(chunked_file, 'w')
            SeqIO.write(seq_record, chunked_fh, 'fasta')
        else:
            chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
            SeqIO.write(seq_record, chunked_fh, 'fasta')
            current_chunk_length += len(seq_record.seq)
    
    chunked_fh.close()
    # elapsed = datetime.now() - starttime
    # logging.info(f'Chunk fasta -- time taken {str(elapsed)}')
    return chunked_file_list

def create_fasta_exclude(fasta_file, exclude_list, output_file):
    with open(output_file,'w') as fh_out:
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            if seq.id not in exclude_list:
                fh_out.write(SeqIO.FastaIO.as_fasta(seq))

    # with open(fasta_file, 'r') as fh_in, open(output_file,'w') as fh_out:
    #     for line in fh_in:
    #         result = re.match(r"^>(\S+)", line)
    #         if result != None:
    #             skip = False
    #             seq_id = result.group(1)
    #             if seq_id in exclude_list:
    #                 skip = True
    #                 continue
    #             fh_out.write(line)
    #         else:
    #             if skip == True:
    #                 continue
    #             else:
    #                 fh_out.write(line)


def create_fasta_include(fasta_file, include_list, output_file):
    with open(output_file,'w') as fh_out:
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            if seq.id in include_list:
                fh_out.write(SeqIO.FastaIO.as_fasta(seq))

    # with open(fasta_file, 'r') as fh_in, open(output_file,'w') as fh_out:
    #     for line in fh_in:
    #         result = re.match(r"^>(\S+)", line)
    #         if result != None:
    #             skip = False
    #             seq_id = result.group(1)
    #             if seq_id not in include_list:
    #                 skip = True
    #                 continue
    #             fh_out.write(line)
    #         else:
    #             if skip == True:
    #                 continue
    #             else:
    #                 fh_out.write(line)

# def translate_protein(nu_fasta, pro_fasta, table):
#     with open(nu_fasta, 'r') as fh_in, open(pro_fasta,'w') as fh_out:
#         for line in fh_in:
#             line = line.rstrip()
#             if re.match(r"^>", line) != None:  
#                 line = re.sub(r'\([-+]\)', '', line)
#                 result = re.match(r"^(>[^:]+)", line)
#                 seq_id = result.group(1)
#             else:
#                 dna = Seq(line)
#                 pro = dna.translate(table=table, stop_symbol='')
#                 pro = str(pro)
                
#                 ls = [pro[i:i+60] for i in range(0,len(pro), 60)]
#                 fh_out.write(seq_id + '\n')
#                 fh_out.write('\n'.join(ls) + '\n')
def getIdentAlignFromCell(cells):
    pident = float(cells[2]) / 100
   
    alignment_length = int(cells[3]) # * 3
    qlen = int(cells[12])# * 3 + 3
    slen = int(cells[13])# * 3 + 3
    short_seq = min(qlen, slen)
    long_seq = max(qlen, slen)
    len_diff = short_seq / long_seq
    align_short = alignment_length / short_seq
    align_long = alignment_length / long_seq
    return pident,len_diff, align_short,align_long,qlen
def concat2fasta(fasta1,fasta2,combined_fasta):
    #new_merged_seq_fasta=os.path.join(out_dir,"unique_concat_seq.fasta")
    cmd=f'cat {fasta1} {fasta2} > {combined_fasta} '
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error concat unique sequences')
    return combined_fasta
def appendTextfile(old_file,new_file):
    #new_merged_seq_fasta=os.path.join(out_dir,"unique_concat_seq.fasta")
    cmd=f'cat {new_file} >> {old_file} '
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error concat unique sequences')
    return old_file
def write_array(filename,data):
    json.dump(data, open(filename, 'w'), indent=4, sort_keys=True)
    #with open(filename, 'wb') as file:
    #    pickle.dump(data, file, protocol=pickle.HIGHEST_PROTOCOL)
    return filename
def read_array(filename):
    #loaded_data=None
    loaded_data=json.load(open(filename, 'r'))
    #with open(filename, 'rb') as file:
    #    loaded_data = pickle.load(file)
    return loaded_data
def save_clusters(clusters,out_dir):
    cluster_dir=os.path.join(out_dir,'clusters')
    if not os.path.exists(cluster_dir):
        os.mkdir(cluster_dir)
    for c in clusters.keys():
        write_array(os.path.join(cluster_dir,c+".seq.json"),clusters[c]['gene_id'])
        clusters[c]['gene_id']=os.path.join(cluster_dir,c+".seq.json")
        write_array(os.path.join(cluster_dir,c+".useq.json"),clusters[c]['unique_seq'])
        clusters[c]['unique_seq']=os.path.join(cluster_dir,c+".useq.json")
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    return os.path.join(out_dir, 'clusters.json')
def save_update_clusters(clusters_file,out_dir):
    cluster_dir=os.path.join(out_dir,'clusters')
    if not os.path.exists(cluster_dir):
        os.mkdir(cluster_dir)
    clusters= json.load(open(clusters_file, 'r'))
    for c in clusters.keys():
        if type(clusters[c]['gene_id']) is list: 
            write_array(os.path.join(cluster_dir,c+".seq.json"),clusters[c]['gene_id'])
            clusters[c]['gene_id']=os.path.join(cluster_dir,c+".seq.json")
        if type(clusters[c]['unique_seq']) is list: 
            write_array(os.path.join(cluster_dir,c+".useq.json"),clusters[c]['unique_seq'])
            clusters[c]['unique_seq']=os.path.join(cluster_dir,c+".useq.json")
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    return os.path.join(out_dir, 'clusters.json')
def save_unique_seqs(unique_groups,out_dir):
    group_dir=os.path.join(out_dir,'groups')
    if os.path.exists(group_dir):
        shutil.rmtree(group_dir)
    if not os.path.exists(group_dir):
        os.mkdir(group_dir)
    for g in unique_groups.keys():
        
        write_array(os.path.join(group_dir,g+".seq.json"),unique_groups[g])
        unique_groups[g]=os.path.join(group_dir,g+".seq.json")
        
    json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    return os.path.join(out_dir, 'unique_groups.json')
def add_unique_groups(groupname,value,out_dir):
    group_dir=os.path.join(out_dir,'groups')
    write_array(os.path.join(group_dir,groupname+".seq.json"),value)
    return os.path.join(group_dir,groupname+".seq.json")

def check_clusters(annotated_clusters_file,out_dir):
    annotated_clusters= json.load(open(annotated_clusters_file, 'r'))
    cluster_dir=os.path.join(out_dir,'clusters')
    for c in annotated_clusters:
        annotated_clusters[c]['updated']=0
        if type(annotated_clusters[c]['gene_id']) is list:
            write_array(os.path.join(cluster_dir,c+".seq.json"),annotated_clusters[c]['gene_id'])
            annotated_clusters[c]['gene_id']=os.path.join(cluster_dir,c+".seq.json")
        if type(annotated_clusters[c]['unique_seq']) is list:
            write_array(os.path.join(cluster_dir,c+".useq.json"),annotated_clusters[c]['unique_seq'])
            annotated_clusters[c]['unique_seq']=os.path.join(cluster_dir,c+".useq.json")
    json.dump(annotated_clusters, open(annotated_clusters_file, 'w'), indent=4, sort_keys=True)
def l2_distance(vector1, vector2):
   
    if len(vector1) != len(vector2):
        raise ValueError("2 vectors not the same lenght")
    
    # Tính khoảng cách L2
    distance = np.sqrt(np.sum((np.array(vector1) - np.array(vector2)) ** 2))
    return distance
def calculate_cosine_distances_np(vectors):
    vectors = np.array(vectors)
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    normalized_vectors = vectors / norms
    cosine_similarity_matrix = np.dot(normalized_vectors, normalized_vectors.T)
    cosine_distance_matrix = 1 - cosine_similarity_matrix
    return cosine_distance_matrix
def check_path(value):
    if isinstance(value,str):
        return True
    else:
        return False 