# -*- coding: utf-8 -*-
"""
    Statistics of sequence, files

Revision history:
----------------
2019-11-05: Amromics created

"""
from __future__ import division, print_function, absolute_import

from amromics.utils.bioseq import read_sequence_file


def fastq_stats(fastq_file):
    """
    Get the statistics of a fastq file
    TODO: quality
    :param fastq_file:
    :return:
    """

    num_read = 0
    num_bases = 0
    for seq in read_sequence_file(fastq_file):
        num_read += 1
        num_bases += len(seq)
    return {'num_read': num_read, 'num_bases': num_bases}


def assembly_stats(fasta_file, min_contig_len=200):
    """
    Get the statistis of an assembly in a fasta file
    :param fasta_file:
    :return:
    """
    ll = []
    contig_count = 0
    for seq in read_sequence_file(fasta_file):
        contig_count += 1
        if len(seq) >= min_contig_len:
            ll.append(len(seq))




