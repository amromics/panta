import os
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
import shutil

logger = logging.getLogger(__name__)


def parse_cluster_file(cd_hit_cluster_file): 
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
    with open(fasta_file, 'r') as fh_in, open(output_file,'w') as fh_out:
        for line in fh_in:
            result = re.match(r"^>(\S+)", line)
            if result != None:
                skip = False
                seq_id = result.group(1)
                if seq_id in exclude_list:
                    skip = True
                    continue
                fh_out.write(line)
            else:
                if skip == True:
                    continue
                else:
                    fh_out.write(line)


def create_fasta_include(fasta_file, include_list, output_file):
    with open(fasta_file, 'r') as fh_in, open(output_file,'w') as fh_out:
        for line in fh_in:
            result = re.match(r"^>(\S+)", line)
            if result != None:
                skip = False
                seq_id = result.group(1)
                if seq_id not in include_list:
                    skip = True
                    continue
                fh_out.write(line)
            else:
                if skip == True:
                    continue
                else:
                    fh_out.write(line)


def translate_protein(nu_fasta, pro_fasta, table):
    with open(nu_fasta, 'r') as fh_in, open(pro_fasta,'w') as fh_out:
        for line in fh_in:
            line = line.rstrip()
            if re.match(r"^>", line) != None:  
                line = re.sub(r'\([-+]\)', '', line)
                result = re.match(r"^(>[^:]+)", line)
                seq_id = result.group(1)
            else:
                dna = Seq(line)
                pro = dna.translate(table=table, stop_symbol='')
                pro = str(pro)
                
                ls = [pro[i:i+60] for i in range(0,len(pro), 60)]
                fh_out.write(seq_id + '\n')
                fh_out.write('\n'.join(ls) + '\n')