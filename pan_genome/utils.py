import os
import logging
import subprocess
import re
from Bio import SeqIO
from Bio.Seq import Seq
import shutil

logger = logging.getLogger(__name__)


def run_command(cmd, timing_log=None):
    if timing_log == None:
        ret = os.system(cmd)
    else:
        logger.info(f'Run: {cmd}')
        cmd = f'/usr/bin/time --append -v -o {timing_log} {cmd}'
        ret = os.system(cmd)

    if ret != 0:
        raise Exception(f'Error running {cmd}')

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

def create_fasta_exclude(fasta_file_list, exclude_list, output_file):
    with open(output_file,'w') as fh_out:
        for fasta_file in fasta_file_list:
            if fasta_file == None:
                continue
            with open(fasta_file, 'r') as fh_in:
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


def create_fasta_include(fasta_file_list, include_list, output_file):
    with open(output_file,'w') as fh_out:
        for fasta_file in fasta_file_list:
            if fasta_file == None:
                continue
            with open(fasta_file, 'r') as fh_in:
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
    premature = []
    startstop = []
    unknown = []
    with open(nu_fasta, 'r') as fh_in, open(pro_fasta,'w') as fh_out:
        for line in fh_in:
            line = line.rstrip()
            if re.match(r"^>", line) != None:  
                line = re.sub(r'\([-+]\)', '', line)
                result = re.match(r"^(>[^:]+)", line)
                seq_id = result.group(1)
            else:
                dna = Seq(line)
                pro = dna.translate(table=table)
                pro = str(pro)
                
                # filter seq with premature codon
                results = re.findall(r'\*', pro)
                if len(results) > 1:
                    premature.append(seq_id)
                    continue
                
                # filter seq lacking start and stop codon
                if pro[0] != 'M' and pro[-1] != '*':
                    startstop.append(seq_id)
                    continue

                # filter seq which has more than 5% of unknown
                results = re.findall(r'X', pro)
                if len(results) / len (pro) > 0.05:
                    unknown.append(seq_id)
                    continue
                
                ls = [pro[i:i+60] for i in range(0,len(pro), 60)]
                fh_out.write(seq_id + '\n')
                fh_out.write('\n'.join(ls) + '\n')
    # if len(premature)!= 0:
    #     logger.info('Have premature codon - exclude {}'.format(', '.join(premature)))
    # if len(startstop)!= 0:
    #     logger.info('Lack both start and stop codon - exclude {}'.format(', '.join(startstop)))
    # if len(unknown)!= 0:
    #     logger.info('Too many unknowns - exclude {}'.format(', '.join(unknown)))