

import re
import subprocess
import gzip


def get_compress_type(filepath):
    with open(filepath, 'rb') as f:
        first_two_bytes = f.read(2)
    if first_two_bytes == b'\x1f\x8b':
        return 'gzip'
    # TODO: fill in for bzip

    # Dont know type
    return None


def get_open_func(filepath):
    """
    Determine compression type (by looking the first 2 bytes) and return the
    right open function to open the file
    Parameters
    ----------
    filepath

    Returns
    -------

    """
    compress_type = get_compress_type(filepath)
    if compress_type == 'gzip':
        return gzip.open
    return open


def valid_id(sid):
    """
    Check if string can be a valid id, that is, it can only the following characters:
        - alphanumerical
        - Underscores
        - Dot
    Parameters
    ----------
    sid

    Returns
    -------

    """
    return re.match(r'^[.\w]+$', sid)


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

    for n in range(0, len(seq), 3):
        if seq[n:n + 3] in codontable:
            residue = codontable[seq[n:n + 3]]
        else:
            residue = "-"

        prot.append(residue)

    return "".join(prot)


def software_version(software_list=None):
    # List of known software and their version
    cmd_versions = {
        # Basic
        'java': 'java -version 2>&1 | head -n 1',
        'python': 'python --version 2>&1',

        # workflow language
        'cromwell': 'cromwell.sh --version 2>&1',

        # Fundamental tools
        'samtools': 'samtools --version  2>&1| head -n 2 | tr "\n" "  "',
        'blast': 'blastn -version 2>&1 | head -n 1',

        # Assemblers
        'spades': 'spades.py -v 2>&1',
        'skesa': 'skesa --version 2>&1 | tail -1',
        'shovill': 'shovill --version 2>&1',

        # Annotations
        'prokka': 'prokka -version 2>&1',
        'mlst': 'mlst --version  2>&1',
        'abricate': 'abricate --version  2>&1|tr "\n" " " && abricate --list|awk \'BEGIN{printf("| Database: ");}NR>1{printf("%s ",$1)}\'',
        'snippy': 'snippy --version 2>&1',

        # Pangenome tools
        'roary': 'roary --version 2>&1 | tail -n 1',
        'parsnp': 'parsnp --version 2>&1 | tail -1',

        # misc
        'trimmomatic': 'trimmomatic -version 2>&1',

        # Others:
        'nodejs':'node --version 2>&1 | tail -1'
        }
    if software_list is None:
        software_list = cmd_versions.keys()

    for sw in software_list:
        if sw in cmd_versions:
            cmd = cmd_versions[sw]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            ret = p.wait()
            if ret == 0:
                print('{:12}= {}'.format(sw, output.decode().strip()))
            else:
                print('{:12}= {}'.format(sw, 'NOT FOUND'))

    for sw in software_list:
        if sw not in cmd_versions:
            print('Cannot check software version for {}'.format(sw))

