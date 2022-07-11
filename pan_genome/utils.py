import os
import logging
import re


from Bio.Seq import Seq


logger = logging.getLogger(__name__)

def check_dir_exist(path):
    """
    Check and raise exception if a directory does not exit.

    Parameters
    ----------
    path : str
        path
    """
    if not os.path.exists(path):
        raise Exception(f'{path} does not exist')

def check_create_folder(path):
    """
    Check and create folder if it does not exit.

    Parameters
    ----------
    path : str
        path
    """
    if not os.path.exists(path):
        os.makedirs(path)


def run_command(cmd, timing_log=None):
    """
    Run external command. Log error if it is failed.

    Parameters
    ----------
    cmd : str
        command
    time_log : path
        path of time.log
    """
    if timing_log == None:
        ret = os.system(cmd)
    else:
        # logger.info(f'Run: {cmd}')
        cmd = f'/usr/bin/time --append -v -o {timing_log} {cmd}'
        ret = os.system(cmd)

    if ret != 0:
        logger.error(f'Error running {cmd}')
        raise Exception(f'Error running {cmd}')

def parse_cluster_file(cd_hit_cluster_file):
    """
    Parse CD-HIT cluster file.

    Parameters
    ----------
    cd_hit_cluster_file : path
        CD-HIT or CD-HIT-2D cluster file.
    
    Returns
    -------
    dictionary
        dictionary of sequence clusters.
        {representative sequence id: [other sequence ids]}
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
        rep = clusters[cluster_name]['representative']
        clusters_new[rep] = clusters[cluster_name]['gene_names']
    return clusters_new


def create_fasta_exclude(fasta_file_list, exclude_list, output_file):
    """
    Create a new fasta from some input fasta, but exclude some sequences.

    Parameters:
    fasta_file_list : list
        list of input fasta file path
    exclude_list : list
        list of excluding sequences
    output_file : path
        path of output fasta
    """
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
    """
    Create a new fasta from some input fasta, only include specific sequences.

    Parameters:
    fasta_file_list : list
        list of input fasta file path
    include_list : list
        list of including sequences
    output_file : path
        path of output fasta
    """
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
    """
    Translate nucleotide fasta file to protein fasta file.
    Filter out sequences:
        - with premature codon
        - lack start and stop codon
        - have more than 5% of unknown

    Parameters
    ----------
    nu_fasta : path
        input nucleotide fasta file
    pro_fasta : path
        output protein fasta file
    table : int
        codon table
    """
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
                    continue
                
                # filter seq lacking start and stop codon
                if pro[0] != 'M' and pro[-1] != '*':
                    continue

                # filter seq which has more than 5% of unknown
                results = re.findall(r'X', pro)
                if len(results) / len (pro) > 0.05:
                    continue
                
                ls = [pro[i:i+60] for i in range(0,len(pro), 60)]
                fh_out.write(seq_id + '\n')
                fh_out.write('\n'.join(ls) + '\n')