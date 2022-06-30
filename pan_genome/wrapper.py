# -*- coding: utf-8 -*-
"""
    Wrapper for pipelines
"""
import os
import subprocess
import logging

from pan_genome import data_preparation
from pan_genome import main_pipeline
from pan_genome import add_sample_pipeline

logger = logging.getLogger(__name__)


def run_init_pipeline(samples, collection_dir, temp_dir, args, timing_log):
    """
    Run initial pan-genome analysis.

    Pipeline description: 
        (1) extract gene sequences
        (2) cluster by CD-HIT
        (3) All-agaist-all comparision by BLASTP
        (4) cluster by Markov clustering algorithms - MCL

    Parameters
    ----------
    samples : list
        list of samples
    collection_dir : path
        collection directory
    temp_dir : path
        temporary directory
    args : object
        Command-line input arguments
    timing_log : path
        path of time.log

    Returns
    -------
    clusters : list of list
        list of sequence clusters
    gene_dictionary : dictionary of {gene_id:(tuple)}
        a data structure contains information of each gene.
        Includes: sample, contig, length, name, product
    """
    gene_dictionary = data_preparation.extract_proteins(
        samples,collection_dir,args)

    combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=temp_dir,
        samples=samples,
        timing_log=timing_log)

    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=combined_faa,
        out_dir=temp_dir,
        threads=args.threads,
        timing_log=timing_log)

    blast_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        timing_log=timing_log,
        evalue = args.evalue,
        max_target_seqs=2000,
        threads=args.threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir, 
        args=args)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        timing_log=timing_log)

    clusters = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_clusters,
        mcl_file=mcl_file)

    return clusters, gene_dictionary


def run_add_pipeline(new_samples, old_represent_faa, previous_clusters, 
                     collection_dir, temp_dir, args, timing_log):
    """
    Add new samples to previous collection.

    Pipeline description: 
        (1) extract gene sequences of new samples
        (2) match new sequence with previous clusters by CD-HIT-2D
        (3) cluster the remain sequences by CD-HIT 
        (4) match new sequence with previous clusters by BLASTP
        (5) All-agaist-all comparision by BLASTP
        (6) cluster by Markov clustering algorithms - MCL

    Parameters
    ----------
    new_samples : list
        list of new samples
    old_represent_faa : path
        path of reference_pangenome.fasta
    previous_clusters : list  of []
        each empty list correspond to a previous cluster
    collection_dir : path
        collection directory
    temp_dir : path
        temporary directory
    args : object
        Command-line input arguments
    timing_log : path
        path of time.log

    Returns
    -------
    new_clusters : list of list
        list of new clusters
    gene_dictionary : dictionary of {gene_id:(tuple)}
        a data structure contains information of each gene.
        Includes: sample, contig, length, name, product
    """
    gene_dictionary = data_preparation.extract_proteins(
        new_samples,collection_dir,args)
    
    new_combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=temp_dir,
        samples=new_samples,
        timing_log=timing_log)

    notmatched_faa, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1 = old_represent_faa,
        database_2 = new_combined_faa,
        out_dir = temp_dir,
        threads=args.threads,
        timing_log=timing_log)

    add_sample_pipeline.add_gene_cd_hit_2d(previous_clusters, cd_hit_2d_clusters)

    num_seq = subprocess.run(
        f'grep ">" {notmatched_faa} | wc -l', 
        capture_output=True, text=True, shell=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return None, None

    notmatched_represent_faa, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=notmatched_faa,
        out_dir=temp_dir,
        threads=args.threads,
        timing_log=timing_log)

    blast_1_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = old_represent_faa,
        query_fasta = notmatched_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast1'),
        timing_log=timing_log,
        evalue = args.evalue,
        max_target_seqs=2000,
        threads=args.threads)

    remain_fasta = add_sample_pipeline.add_gene_blast(
        previous_clusters=previous_clusters,
        cd_hit_clusters = cd_hit_clusters,
        blast_result=blast_1_result, 
        fasta_file=notmatched_represent_faa, 
        out_dir=temp_dir,
        args=args)

    num_seq = subprocess.run(
        f'grep ">" {remain_fasta} | wc -l', 
        capture_output=True, text=True, shell=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return None, None
    
    blast_2_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = remain_fasta,
        query_fasta = remain_fasta,
        out_dir = os.path.join(temp_dir, 'blast2'),
        timing_log=timing_log,
        evalue = args.evalue,
        max_target_seqs=2000,
        threads=args.threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_2_result,
        out_dir = temp_dir, 
        args=args)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        timing_log=timing_log)

    new_clusters = main_pipeline.reinflate_clusters(
        cd_hit_clusters = cd_hit_clusters,
        mcl_file=mcl_file)

    return new_clusters, gene_dictionary