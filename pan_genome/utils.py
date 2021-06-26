import os
import logging
import re
from Bio import SeqIO
import shutil
from datetime import datetime

logger = logging.getLogger(__name__)

def run_command(cmd, timing_log=None):
    """
    Run a command line, return the returning code of the command
    :param cmd:
    :param timing_log:
    :return:
    """
    #logger.info('Run "' + cmd + '"')
    if timing_log is not None:
        cmd = '/usr/bin/time --append -v -o {} bash -c "{}"'.format(timing_log, cmd)
    ret = os.system(cmd)
    return ret


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
            result = re.match(r"^>(\w+)", line)
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
            result = re.match(r"^>(\w+)", line)
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


def translate_dna(sequence):
    """
    :param sequence: (str) a DNA sequence string
    :return: (str) a protein string from the forward reading frame 1
    """

    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
                }
    seq = sequence.upper()
    prot = []
    num_unknown = 0
    for n in range(0, len(seq), 3):
        codon = seq[n:n + 3]
        if codon in codontable:
            residue = codontable[codon]
        else:
            residue = "X"
            num_unknown += 1
        prot.append(residue)
    
    # filter seq which has more than 5% of unknown
    if num_unknown / len (prot) > 0.05:
        return None
    
    return "".join(prot)

def translate_protein(nu_fasta, pro_fasta):
    with open(nu_fasta, 'r') as fh_in, open(pro_fasta,'w') as fh_out:
        for line in fh_in:
            line = line.rstrip()
            if re.match(r"^>", line) != None:  
                line = re.sub(r'\([-+]\)', '', line)
                result = re.match(r"^(>[^:]+)", line)
                seq_id = result.group(1)
            else:
                line = translate_dna(line)
                if line == None:
                    # logger.info(f'Exclude {seq_id} - too many unknowns')
                    continue
                ls = [line[i:i+60] for i in range(0,len(line), 60)]
                fh_out.write(seq_id + '\n')
                fh_out.write('\n'.join(ls) + '\n')
