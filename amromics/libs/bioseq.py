# -*- coding: utf-8 -*-
"""
    Define sample class


Revision history:
----------------
2019-08-17: Amromics created

"""
from __future__ import division, print_function, absolute_import

import gzip
import bz2
from collections import defaultdict

complement_dict = defaultdict(lambda: 'N')
complement_dict['A'] = 'T'
complement_dict['C'] = 'G'
complement_dict['G'] = 'C'
complement_dict['T'] = 'A'

STOP_AA = '_'
codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
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
        'TAC': 'Y', 'TAT': 'Y', 'TAA': STOP_AA, 'TAG': STOP_AA,
        'TGC': 'C', 'TGT': 'C', 'TGA': STOP_AA, 'TGG': 'W',
    }


def complement(seq):
    bases = list(seq)
    comp = [complement_dict[b] for b in bases]
    return ''.join(comp)


def reverse_complement(seq):
    return complement(seq[::-1])


def identical_protein(s1, s2):
    if len(s1) != len(s2):
        return False
    for i, aa in enumerate(s1):
        if not identical_aa(aa, s2[i]):
            return False
    return True


def identical_aa(x,y):
    """
    Compare two amino acids
    :param x:
    :param y:
    :return:
    """
    if x == y:
        return True
    if x == 'X' or y == 'X':
        return True
    return False


def codon_to_aa(codon):
    if codon in codon_table.keys():
        return codon_table[codon]
    else:
        return 'X'


def _translate_orf(cds):
    """
    Translate the CDS to protein, starting at position 0 (assuming ORF=0)
    :param cds:
    :return:
    """
    protein = ''
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i + 3]
        aa = codon_to_aa(codon)
        protein += aa
        if aa == STOP_AA:
            break
    return protein


def translate(cds, reverse=False):
    """
    Translate the CDS to protein, try 3 ORFs and pick the longest one.
    If reverse is set to True, also try the reverse complement
    :param cds:
    :return:
    """
    ORF=0
    trans = _translate_orf(cds)

    _trans = _translate_orf(cds[1:])
    if len(_trans) > len(trans):
        ORF = 1
        trans = _trans

    _trans = _translate_orf(cds[2:])
    if len(_trans) > len(trans):
        ORF = 2
        trans = _trans

    if reverse:
        cds = reverse_complement(cds)
        _trans = _translate_orf(cds)
        if len(_trans) > len(trans):
           ORF = -1
           trans = _trans

        _trans = _translate_orf(cds[1:])
        if len(_trans) > len(trans):
           ORF = -2
           trans = _trans

        _trans = _translate_orf(cds[2:])
        if len(_trans) > len(trans):
           ORF = -3
           trans = _trans
    return trans, ORF


class Sample:
    """
    Definition of a sample
    """
    def __init__(self, sample_id='AB0000001', genus='Escherichia',  species='coli', strain='', gram=''):
        self.sample_id = sample_id
        self.genus = genus
        self.species = species
        self.strain = strain
        #self.path = self.genus + '/' + self.genus[0] + '.' + self.species + '/' + self.sample_id
        self.path = self.sample_id
        self.gram = gram

    def get_path(self):
        return self.path

    def get_id(self):
        return self.sample_id


class Sequence(object):
    def __init__(self, name, sequence='', desc='', quality=None):
        self.quality = quality
        self.name = name
        self.sequence = sequence
        self.desc = desc

    def __len__(self):
        return len(self.sequence)

    def length(self):
        return len(self.sequence)

    def set_name(self, name):
        self.name = name

    def set_desc(self, desc):
        self.desc = desc

    def get_name(self):
        return self.name

    def get_desc(self):
        return self.desc

    def same_sequence(self, seq):
        return self.sequence == seq.sequence

    def format_fasta(self, max_line=100):
        """
        Format the sequence in fasta format
        :param max_line:
        :return:
        """
        ret = '>' + self.name + (' ' + self.desc if self.desc else '') + '\n'
        for i in range(0, self.length(), max_line):
            ret += self.sequence[i:i+max_line] + '\n'
        return ret

    def format_fastq(self):
        return '@' + self.name + (' ' + self.desc if self.desc else '') + '\n' + self.sequence \
               + '\n+\n' + (self.quality if self.quality else '') + '\n'


class AMRGene():
    def __init__(self, amr_id, seq):
        self.gene_sequence = seq
        self.cds_start = 0
        self.cds_end = 0
        self.aka = {}

    def add_alias(self, database, id):
        if database in self.aka and self.aka[database] != id:
            raise Exception('Existing id {} found for {}'.format(self.aka[database], id))
        self.aka[database] = id


#Get the compression type of a file
#https://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type-and-uncompress
magic_dict = {
    "\x1f\x8b\x08": "gz",
    "\x42\x5a\x68": "bz2",
    "\x50\x4b\x03\x04": "zip"
    }
_magic_max_len = max(len(x) for x in magic_dict)


def get_compression_type(filename):
    with open(filename) as f:
        file_start = f.read(_magic_max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return None


def readfq(fp):
    """
    Fast fastq/fasta reader written by Heng Li
    See: https://github.com/lh3/readfq/blob/master/readfq.py
    :param fp:
    :return:
    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        partitions = last[1:].partition(" ")
        name, desc, seqs, last = partitions[0], partitions[2], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, desc, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, desc, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, desc, seq, None  # yield a fasta record instead
                break


def read_sequence_file(file_name):
    """
    Read a sequence file in either fastq of fasta format. The file can be gzip'ed/b2zip'ed
    :param file_name:
    :return:
    """
    compression_type = get_compression_type(file_name)

    if compression_type is None:
        open_method = open
    elif compression_type == 'gz':
        open_method = gzip.open
    elif compression_type == 'bz2':
        open_method = bz2.BZ2File
    else:
        raise Exception('Unknown compression type {}'.format(compression_type))

    with open_method(file_name, 'r') as f:
        for name, desc, seq, qual in readfq(f):
            yield Sequence(name, desc=desc, sequence=seq, quality=qual)


def write_fasta(file_name, seqs, max_line=100):
    """
    Write
    :param file_name:
    :param seqs:
    :param max_line:
    :return:
    """
    with open(file_name, 'w') as fc:
        for seq in seqs:
            fc.write(seq.format_fasta(max_line=max_line))


def write_fastq(file_name, seqs):
    with open(file_name, 'w') as fc:
        for seq in seqs:
            fc.write(seq.format_fastq())

