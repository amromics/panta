#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import multiprocessing
import logging
import json
import csv
from datetime import datetime
from panta import *
from panta.utils import *
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import tracemalloc
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)


def collect_sample(sample_id_list, args):
    samples = []
    if args.tsv != None:
        with open(args.tsv,'r') as fh:
            csv_reader = csv.reader(fh, delimiter='\t')
            for row in csv_reader:
                gff = row[1]
                if (not gff.endswith('.gff')) and (not gff.endswith('.gff.gz')):
                    raise Exception(f'{gff} should be a gff3 file (file ending with .gff or .gff.gz')
                sample_id = row[0].replace('-','_')#Make sure that - is not part of sample_id

                if sample_id in sample_id_list:
                    logging.info(f'{sample_id} already exists -- skip')
                    continue
                else:
                    sample_id_list.append(sample_id)
                assembly = row[2]
                if row[2] == '':
                    assembly = None
                samples.append({'id':sample_id, 'gff_file':gff, 'assembly':assembly})

    elif args.gff != None:
        gff_list = args.gff
        for gff in gff_list:
            base_name = os.path.basename(gff)
            if gff.endswith('.gff') or gff.endswith('.GFF'):
                sample_id = base_name[:-4]
            elif gff.endswith('.gff.gz') or gff.endswith('.GFF.gz'):
                sample_id = base_name[:-7]
            elif gff.endswith('.gff3') or gff.endswith('.GFF3'):
                sample_id = base_name[:-5]
            elif gff.endswith('.gff3.gz') or gff.endswith('.GFF3.gz'):
                sample_id = base_name[:-8]
            else:
                raise Exception(f'{gff} file should have one of suffices .gff, .gff3, .GFF, .GFF3')

            sample_id = sample_id.replace('-','_')#Make sure that - is not part of sample_id
            if sample_id in sample_id_list:
                logging.info(f'{sample_id} already exists -- skip')
                continue
            else:
                sample_id_list.append(sample_id)
            samples.append({'id':sample_id, 'gff_file':gff, 'assembly':None})
    else:
        raise Exception(f'Please specify -t or -g')

    #samples.sort(key= lambda x:x['id'])
    return samples

def run_main_pipeline_a_diamond_c_mcl(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    
    blast_result = main_pipeline.pairwise_alignment_diamond(
        
        database_fasta = combined_faa,
        query_fasta = combined_faa,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.make_clusters_from_mcl(
        
        mcl_file=mcl_file,
        map_file=combined_faa_map
        )
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples,samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)

    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
    #if os.path.exists(temp_dir):
    #    shutil.rmtree(temp_dir)
def run_main_pipeline_g_mmseq_a_diamond_c_mcl(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_groups)}')  """
    cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_mmseq_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len mmseq_clusters = {len(cd_hit_groups)}') 

    #print(f'Diamond = {args.diamond}')
    blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_groups,
        mcl_file=mcl_file)
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples,samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
        post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,threads=threads)
        post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_diamond_a_diamond_c_mcl(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_clusters)}') """
    diamond_represent_fasta, diamond_groups = main_pipeline.run_diamond_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        cover=args.cov,
        threads=threads)
    logger.info(f'len diamond_clusters = {len(diamond_groups)}') 

    #print(f'Diamond = {args.diamond}')
    blast_result = main_pipeline.pairwise_alignment_diamond(
        #diamond=(args.blast=='diamond'),
        database_fasta = diamond_represent_fasta,
        query_fasta =diamond_represent_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        groups=diamond_groups,
        mcl_file=mcl_file)
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post 
    
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples,samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)
        post_analysis.run_hhm_profile_clusters(annotated_clusters,out_dir,  threads=threads)
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.move(diamond_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(inflated_clusters, open(os.path.join(out_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
    #if os.path.exists(temp_dir):
    #    shutil.rmtree(temp_dir)
def run_main_pipeline_g_cdhit_a_diamond_c_mcl(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_groups)}') 
    """ cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_mmseq_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len mmseq_clusters = {len(cd_hit_groups)}')  """

    #print(f'Diamond = {args.diamond}')
    blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        groups=cd_hit_groups,
        mcl_file=mcl_file)
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples, samples,out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)

    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_c_diamond(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_clusters)}') """
    inflated_clusters = main_pipeline.run_diamond_clustering(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    #logger.info(f'len mmseq_clusters = {len(cd_hit_groups)}') 

   

    


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples, samples,out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)

    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    #json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_diamond_a_diamond_c_diamond(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_clusters)}') """
    inflated_clusters = main_pipeline.run_diamond_clustering_pipeline(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    #logger.info(f'len mmseq_clusters = {len(cd_hit_groups)}') 

   

    


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)

    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    #json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_c_mmseq(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_clusters)}') """
    inflated_clusters = main_pipeline.run_mmseq_clustering(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    #logger.info(f'len mmseq_clusters = {len(cd_hit_groups)}') 

   

    


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples,samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)

    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    #json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_m_ref_g_mmseq_a_diamond_c_mcl(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    ref_db='gene_families/gene_families_diamond_db.dmnd'
    ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    group_seq_dir=os.path.join(out_dir, 'groups')
    os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_groups)}')  """
    mmseq_represent_fasta, mmseq_groups = main_pipeline.run_mmseq_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        group_dir =group_seq_dir,       
        threads=threads)
    not_match_represent_faa,clusters_by_ref,remain_groups= main_pipeline.match_seqs_to_ref(
        out_dir = temp_dir,
        seqs_file = mmseq_represent_fasta,
        
        groups=mmseq_groups,
        refdb=ref_db,
        ref_clusters=ref_clusters,
        
        threads=threads)
    json.dump(clusters_by_ref, open(os.path.join(out_dir, 'clusters_by_ref.json'), 'w'), indent=4, sort_keys=True)

    # mmseq_represent_fasta2, mmseq_groups2 = main_pipeline.run_mmseq_with_map(
    #     faa_file=not_match_faa,
    #     map_file=combined_faa_map,
    #     out_dir=temp_dir,
    #     threads=threads)
    # logger.info(f'len new mmseq_clusters = {len(mmseq_groups2)}') 

    #print(f'Diamond = {args.diamond}')
    blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = not_match_represent_faa,
        query_fasta = not_match_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        groups=remain_groups,
        mcl_file=mcl_file)
    logger.info(f'len new inflated_clusters = {len(inflated_clusters)} len groups = {len(clusters)}')
    json.dump(inflated_clusters, open(os.path.join(out_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len new split_clusters = {len(split_clusters)}')
    json.dump(split_clusters, open(os.path.join(out_dir, 'split_clusters.json'), 'w'), indent=4, sort_keys=True)

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)
    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters0.json'), 'w'), indent=4, sort_keys=True)

    annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    #post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
 
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_m_ref_g_mmseq_a_diamond_c_mcl_nearest_group(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()
    print("using "+str(threads) +" CPU cores")
    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    ref_db='gene_families/gene_families_diamond_db.dmnd'
    ref_cdb='gene_families/db'
    ref_clusters = json.load(open('gene_families/gene_family_splited_clusters.json','r'))
    group_seq_dir=os.path.join(out_dir, 'groups')
    os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_groups)}')  """
    mmseq_represent_fasta, mmseq_groups = main_pipeline.run_mmseq_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        group_dir=group_seq_dir,
        threads=threads)
    not_match_represent_faa,clusters_by_ref,remain_groups= main_pipeline.match_seqs_to_ref_by_nearest_group(
        out_dir = temp_dir,
        seqs_file = mmseq_represent_fasta,
        refcdb=ref_cdb,
        groups=mmseq_groups,
        refdb=ref_db,
        ref_clusters=ref_clusters,       
        threads=threads)
    # mmseq_represent_fasta2, mmseq_groups2 = main_pipeline.run_mmseq_with_map(
    #     faa_file=not_match_faa,
    #     map_file=combined_faa_map,
    #     out_dir=temp_dir,
    #     threads=threads)
    # logger.info(f'len new mmseq_clusters = {len(mmseq_groups2)}') 
    #print(f'Diamond = {args.diamond}')
    blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = not_match_represent_faa,
        query_fasta = not_match_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        inflation=2,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters_by_groups(
        groups=remain_groups,
        mcl_file=mcl_file)
    logger.info(f'len new inflated_clusters = {len(inflated_clusters)} len groups = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs_by_group(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len new split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster_main(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)

    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_m_ref_s_diamond_progressive(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()
    print("using "+str(threads) +" CPU cores")
    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    ref_db='gene_families/gene_families_diamond_db.dmnd'
    ref_cdb='gene_families/db'
    ref_clusters = json.load(open('gene_families/gene_family_splited_clusters.json','r'))
    group_seq_dir=os.path.join(out_dir, 'groups')
    os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)
    for sample in samples:
        list_seqs=parse_sample_to_list_seq(sample,out_dir)


    annotated_clusters = post_analysis.annotate_cluster_main(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    
    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_mmseq_a_diamond_c_mcl_rand_old(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    ref_db='gene_families/gene_families_diamond_db.dmnd'
    ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_groups)}')  """
    #firsly, grouping with i=1 and c=1 to get unique seds
    unique_seds_fasta, unique_groups = main_pipeline.run_mmseq_with_map_unique_seqs(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,      
        threads=threads)
    json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)

    #grouping with 98% to get high similar groups 
    groups_representative_fasta, similar_groups = main_pipeline.run_mmseq_with_map_similar_seqs(
        faa_file=unique_seds_fasta,
        
        out_dir=temp_dir,      
        threads=threads)
    json.dump(similar_groups, open(os.path.join(out_dir, 'similar_groups.json'), 'w'), indent=4, sort_keys=True)

    # mmseq_represent_fasta2, mmseq_groups2 = main_pipeline.run_mmseq_with_map(
    #     faa_file=not_match_faa,
    #     map_file=combined_faa_map,
    #     out_dir=temp_dir,
    #     threads=threads)
    # logger.info(f'len new mmseq_clusters = {len(mmseq_groups2)}') 

    #print(f'Diamond = {args.diamond}')
    blast_result = main_pipeline.pairwise_alignment_partion_diamond(
      
        #database_fasta = groups_representative_fasta,
        input_fasta = groups_representative_fasta,
        out_dir = os.path.join(temp_dir),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters_with_unique_seqs(
        unique_groups=unique_groups,
        similar_groups=similar_groups,
        mcl_file=mcl_file)
    logger.info(f'len new inflated_clusters = {len(inflated_clusters)} len groups = {len(clusters)}')
    json.dump(inflated_clusters, open(os.path.join(out_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)


    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len new split_clusters = {len(split_clusters)}')
    json.dump(split_clusters, open(os.path.join(out_dir, 'split_clusters.json'), 'w'), indent=4, sort_keys=True)

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)
    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)

    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    #post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    #post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    #post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
 
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seds_fasta, main_unique_seds)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(cd_hit_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_mmseq_a_diamond_c_mcl_rand(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    """ cd_hit_represent_fasta, cd_hit_groups = main_pipeline.run_cd_hit_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len cd_hit_clusters = {len(cd_hit_groups)}')  """
    groups_represent_fasta, similar_groups = main_pipeline.run_mmseq_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        threads=threads)
    logger.info(f'len mmseq_clusters = {len(similar_groups)}') 

    #print(f'Diamond = {args.diamond}')
    sourmash_result_file = main_pipeline.pairwise_alignment_sourmash(
      
        #database_fasta = groups_representative_fasta,
        
        query_fasta = groups_represent_fasta,
        out_dir = out_dir,
        ksize=args.ksize,
        similarity=args.similarity,
        diff_len=args.LD,
        threads=threads,
        timing_log=None)
   
    mcl_file = main_pipeline.cluster_with_mcl_from_faiss(
        out_dir = out_dir,
        file_faiss = sourmash_result_file,
        threads=threads,
        inflation=2,
        timing_log=None)
    

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster_include_groups(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    

    output.create_outputs1(annotated_clusters,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    if args.alignment:
        post_analysis.run_gene_alignment(annotated_clusters, samples,samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
        post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,threads=threads)
        post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')

    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)

    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.move(groups_represent_fasta, os.path.join(out_dir, 'representative.fasta'))
    json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_mmseq_a_diamond_c_mcl_connected_componnet(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    #ref_db='gene_families/gene_families_diamond_db.dmnd'
    #ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=out_dir,
        samples=samples)

    
    unique_seds_fasta, unique_groups = main_pipeline.run_mmseq_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    
    #json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    
    

    inflated_clusters, similar_groups, representative_similar_groups,represent_clusters_file=main_pipeline.clustering(
        input_fasta_file=unique_seds_fasta,
        out_dir=temp_dir,
        threads=threads
        )

   
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn)
    json.dump(similar_groups, open(os.path.join(out_dir, 'similar_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'before_annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    annotated_clusters_file=post_analysis.expandClusterMembers(annotated_clusters_file,unique_groups)
    unique_groups=save_unique_seqs(unique_groups,out_dir)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs(annotated_clusters_file,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    #post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    #post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    #post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
    represent_clusters_file=main_pipeline.make_representative_clusters_from_similar_groups(
        annotated_clusters_file=annotated_clusters_file,
        representative_groups=representative_similar_groups,
        out_dir=out_dir
    )
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
   
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seds_fasta, main_unique_seds)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.rmtree(os.path.join(out_dir, 'samples'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_mmseq_a_diamond_c_mcl_core_analysis(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    #ref_db='gene_families/gene_families_diamond_db.dmnd'
    #ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=out_dir,
        samples=samples)

    
    unique_seds_fasta, unique_groups = main_pipeline.run_mmseq_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    
    #json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    
    

    inflated_clusters, similar_groups, representative_similar_groups=main_pipeline.clustering(
        input_fasta_file=unique_seds_fasta,
        out_dir=temp_dir,
        threads=threads
        )

   
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn)
    json.dump(similar_groups, open(os.path.join(out_dir, 'similar_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'before_annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    annotated_clusters_file=post_analysis.expandClusterMembers(annotated_clusters_file,unique_groups)
    unique_groups=save_unique_seqs(unique_groups,out_dir)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs(annotated_clusters_file,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    #post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    #post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    #post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
    represent_clusters_file=main_pipeline.make_representative_clusters_from_similar_groups(
        annotated_clusters_file=annotated_clusters_file,
        representative_groups=representative_similar_groups,
        out_dir=out_dir
    )
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
   
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seds_fasta, main_unique_seds)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.rmtree(os.path.join(out_dir, 'samples'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_mmseq_a_diamond_c_mcl_hash(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    #ref_db='gene_families/gene_families_diamond_db.dmnd'
    #ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=out_dir,
        samples=samples)

    
    unique_seqs_fasta, gene_hash, map_rep_hash= main_pipeline.run_hash_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    
    #json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    
    

    inflated_clusters, similar_groups, representative_similar_groups=main_pipeline.clustering(
        input_fasta_file=unique_seqs_fasta,
        out_dir=temp_dir,
        threads=threads
        )

    group_hash={}
    for g in similar_groups:
        h=map_rep_hash[g]
        group_hash[h]=[h]
        for id in similar_groups[g]:
            group_hash[h].append(map_rep_hash[id])
    json.dump(inflated_clusters, open(os.path.join(out_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters = post_analysis.annotate_cluster_hash(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn,
        gene_hash=gene_hash,
       
        map_gene_hash=map_rep_hash)
    #json.dump(similar_groups, open(os.path.join(out_dir, 'similar_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'before_annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters_file=os.path.join(out_dir, 'clusters.json')
    
    json.dump(annotated_clusters, open(annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    #annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    #annotated_clusters_file=post_analysis.expandClusterMembersByHash(annotated_clusters_file,gene_hash,map_rep_hash)
    json.dump(gene_hash, open(os.path.join(out_dir, 'gene_hash.json'), 'w'), indent=4, sort_keys=True)
    group_hash_file=os.path.join(out_dir, 'group_hash.json')
    
    json.dump(group_hash, open(group_hash_file, 'w'), indent=4, sort_keys=True)
    
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs_from_hash(annotated_clusters_file,gene_hash,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    
    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
    main_similar = os.path.join(out_dir, 'similar.tsv')
    main_similar_fasta=os.path.join(out_dir, 'similar.fasta')
    shutil.move(representative_similar_groups, main_similar_fasta)
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seqs_fasta, main_unique_seds)
    shutil.move(os.path.join(temp_dir, 'filtered_blast_results'), main_similar)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.rmtree(os.path.join(out_dir, 'samples'))
  

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_mmseq_a_diamond_c_mcl_hash_opt(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    hash_dir = os.path.join(out_dir, 'hash')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    #ref_db='gene_families/gene_families_diamond_db.dmnd'
    #ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=out_dir,
        samples=samples)

    
    unique_seqs_fasta, gene_hash, map_rep_hash= main_pipeline.run_hash_unique_seqs_opt(
        faa_file=combined_faa,

        out_dir=temp_dir,      
        hash_dir=hash_dir,
        threads=threads)
    
    #json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    
    

    inflated_clusters, similar_groups, representative_similar_groups=main_pipeline.clustering(
        input_fasta_file=unique_seqs_fasta,
        out_dir=temp_dir,
        threads=threads
        )

    group_hash={}
    for g in similar_groups:
        h=map_rep_hash[g]
        group_hash[h]=[h]
        for id in similar_groups[g]:
            group_hash[h].append(map_rep_hash[id])
    json.dump(inflated_clusters, open(os.path.join(out_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters = post_analysis.annotate_cluster_hash_opt(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn,
       
       
        map_gene_hash=map_rep_hash)
    #json.dump(similar_groups, open(os.path.join(out_dir, 'similar_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'before_annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters_file=os.path.join(out_dir, 'clusters.json')
    
    json.dump(annotated_clusters, open(annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    #annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    #annotated_clusters_file=post_analysis.expandClusterMembersByHash(annotated_clusters_file,gene_hash,map_rep_hash)
    json.dump(gene_hash, open(os.path.join(out_dir, 'gene_hash.json'), 'w'), indent=4, sort_keys=True)
    group_hash_file=os.path.join(out_dir, 'group_hash.json')
    map_rep_hash_file=os.path.join(out_dir, 'map_rep_hash.json')
    json.dump(group_hash, open(group_hash_file, 'w'), indent=4, sort_keys=True)
    json.dump(map_rep_hash, open(map_rep_hash_file, 'w'), indent=4, sort_keys=True)
    
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs_from_hash_opt(annotated_clusters_file,gene_hash,samples,out_dir,hash_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    
    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
    main_similar = os.path.join(out_dir, 'similar.tsv')
    main_similar_fasta=os.path.join(out_dir, 'similar.fasta')
    shutil.move(representative_similar_groups, main_similar_fasta)
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seqs_fasta, main_unique_seds)
    shutil.move(os.path.join(temp_dir, 'filtered_blast_results'), main_similar)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.rmtree(os.path.join(out_dir, 'samples'))
  

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_g_mmseq_a_diamond_c_mcl_hash_opt2(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    #ref_db='gene_families/gene_families_diamond_db.dmnd'
    #ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa,combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=out_dir,
        samples=samples)

    
    unique_seqs_fasta, gene_hash, map_rep_hash= main_pipeline.run_hash_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    
    #json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    
    

    inflated_clusters, similar_groups, representative_similar_groups=main_pipeline.clustering(
        input_fasta_file=unique_seqs_fasta,
        out_dir=temp_dir,
        threads=threads
        )

    group_hash={}
    for g in similar_groups:
        h=map_rep_hash[g]
        group_hash[h]=[h]
        for id in similar_groups[g]:
            group_hash[h].append(map_rep_hash[id])
    json.dump(inflated_clusters, open(os.path.join(out_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    no_annotated_clusters = post_analysis.make_no_annotated_cluster(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn,
        gene_hash=gene_hash,
       
        map_gene_hash=map_rep_hash)
    #json.dump(similar_groups, open(os.path.join(out_dir, 'similar_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'before_annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    clusters_file=os.path.join(out_dir, 'clusters.json')
    
    json.dump(no_annotated_clusters, open(clusters_file, 'w'), indent=4, sort_keys=True)
    #annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    #annotated_clusters_file=post_analysis.expandClusterMembersByHash(annotated_clusters_file,gene_hash,map_rep_hash)
    json.dump(gene_hash, open(os.path.join(out_dir, 'gene_hash.json'), 'w'), indent=4, sort_keys=True)
    group_hash_file=os.path.join(out_dir, 'group_hash.json')
    
    json.dump(group_hash, open(group_hash_file, 'w'), indent=4, sort_keys=True)
    
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #output.create_outputs_from_hash(annotated_clusters_file,gene_hash,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    
    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
    main_similar = os.path.join(out_dir, 'similar.tsv')
    main_similar_fasta=os.path.join(out_dir, 'similar.fasta')
    main_combine_map=os.path.join(out_dir, 'combined.map')
    shutil.move(representative_similar_groups, main_similar_fasta)
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seqs_fasta, main_unique_seds)
    shutil.move(combined_faa_map, main_combine_map)
    shutil.move(os.path.join(temp_dir, 'filtered_blast_results'), main_similar)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.rmtree(os.path.join(out_dir, 'samples'))
  

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_faiss_mmseq_mcl(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    #ref_db='gene_families/gene_families_diamond_db.dmnd'
    #ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=out_dir,
        samples=samples)

    
    # unique_seds_fasta, unique_groups = main_pipeline.run_faiss_unique_seqs(
    #     faa_file=combined_faa,
      
    #     out_dir=temp_dir,      
    #     threads=threads)
    unique_seds_fasta, unique_groups = main_pipeline.run_mmseq_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    #json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    
    

    inflated_clusters, similar_groups, representative_similar_groups,represent_clusters_file=main_pipeline.clustering_faiss(
        input_fasta_file=unique_seds_fasta,
        out_dir=temp_dir,
        threads=threads
        )

   
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'before_annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    annotated_clusters_file=post_analysis.expandClusterMembers(annotated_clusters_file,unique_groups)
    unique_groups=save_unique_seqs(unique_groups,out_dir)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs(annotated_clusters_file,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    #post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    #post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    #post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
    represent_clusters_file=main_pipeline.make_representative_clusters_from_similar_groups(
        annotated_clusters_file=annotated_clusters_file,
        representative_groups=representative_similar_groups,
        out_dir=out_dir
    )
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
   
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seds_fasta, main_unique_seds)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.rmtree(os.path.join(out_dir, 'samples'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main_pipeline_ems_mmseq_mcl(args):
    starttime = datetime.now()

    out_dir = args.outdir
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    temp_dir = os.path.join(out_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    #ref_db='gene_families/gene_families_diamond_db.dmnd'
    #ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    #group_seq_dir=os.path.join(out_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # collect samples
    sample_id_list = []
    samples = collect_sample(sample_id_list, args)
    if len(samples) < 2:
        raise Exception(f'There must be at least 2 samples')

    data_preparation.extract_proteins_tofile(
        samples=samples,
        out_dir=out_dir,
        gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        table=args.table,
        threads=threads)

    # combined_faa = data_preparation.combine_proteins(
    #     out_dir=out_dir,
    #     samples=samples)

    # main_pipeline
    # cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
    #     faa_file=combined_faa,
    #     out_dir=temp_dir,
    #     threads=threads)

    combined_faa = data_preparation.combine_proteins(
        out_dir=out_dir,
        samples=samples)

    
    # unique_seds_fasta, unique_groups = main_pipeline.run_faiss_unique_seqs(
    #     faa_file=combined_faa,
      
    #     out_dir=temp_dir,      
    #     threads=threads)
    unique_seds_fasta, unique_groups = main_pipeline.run_mmseq_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    #json.dump(unique_groups, open(os.path.join(out_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    
    

    inflated_clusters, similar_groups, representative_similar_groups,represent_clusters_file=main_pipeline.clustering_ems(
        input_fasta_file=unique_seds_fasta,
        out_dir=temp_dir,
        threads=threads
        )

   
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'before_annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    annotated_clusters_file=post_analysis.expandClusterMembers(annotated_clusters_file,unique_groups)
    unique_groups=save_unique_seqs(unique_groups,out_dir)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs(annotated_clusters_file,samples,out_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #post_analysis.run_gene_alignment(annotated_clusters, samples, out_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa)
    #post_analysis.create_poa_protein_consensus(annotated_clusters,out_dir,gene_families_dir='gene_families',threads=threads)
    #post_analysis.make_protein_consensus_db(annotated_clusters,out_dir,threads=threads)
    #post_analysis.create_protein_db_for_clusters(annotated_clusters,out_dir,group_seq_dir,threads=threads)
    represent_clusters_file=main_pipeline.make_representative_clusters_from_similar_groups(
        annotated_clusters_file=annotated_clusters_file,
        representative_groups=representative_similar_groups,
        out_dir=out_dir
    )
    
    # output for next run
    #output.export_gene_annotation(gene_annotation, out_dir)
    #json.dump(gene_position, open(os.path.join(out_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)

    main_gene_annotation_fn = os.path.join(out_dir, 'gene_annotation.csv')
    main_gene_position_fn = os.path.join(out_dir, 'gene_position.csv')
    main_unique_seds = os.path.join(out_dir, 'unique_seqs.fasta')
   
    shutil.move(gene_annotation_fn, main_gene_annotation_fn)
    shutil.move(gene_position_fn, main_gene_position_fn)
    shutil.move(unique_seds_fasta, main_unique_seds)
    json.dump(samples, open(os.path.join(out_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    shutil.rmtree(os.path.join(out_dir, 'samples'))
    # json.dump(clusters, open(os.path.join(out_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.copy(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #shutil.move(blast_result, os.path.join(out_dir, 'blast.tsv'))
    #cmd = f'gzip -c {blast_result} > ' + os.path.join(out_dir, 'blast.tsv.gz')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_main(args):
    if args.mode=='a_diamond_c_mcl':
        run_main_pipeline_a_diamond_c_mcl(args)
    if args.mode=='g_cdhit_a_diamond_c_mcl':
        run_main_pipeline_g_cdhit_a_diamond_c_mcl(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl':
        run_main_pipeline_g_mmseq_a_diamond_c_mcl(args)
    if args.mode=='g_diamond_a_diamond_c_mcl':
        run_main_pipeline_g_diamond_a_diamond_c_mcl(args)
    if args.mode=='c_mmseq':
        run_main_pipeline_c_mmseq(args)
    if args.mode=='c_diamond':
        run_main_pipeline_c_diamond(args)
    if args.mode=='g_diamond_a_diamond_c_diamond':
        run_main_pipeline_g_diamond_a_diamond_c_diamond(args)
    if args.mode=='m_ref_g_mmseq_a_diamond_c_mcl':
        run_main_pipeline_m_ref_g_mmseq_a_diamond_c_mcl(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_rand':
        run_main_pipeline_g_mmseq_a_diamond_c_mcl_rand(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_cc':
        run_main_pipeline_g_mmseq_a_diamond_c_mcl_connected_componnet(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_hash':
        run_main_pipeline_g_mmseq_a_diamond_c_mcl_hash(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_hash_opt':
        run_main_pipeline_g_mmseq_a_diamond_c_mcl_hash_opt(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_hash_opt2':
        run_main_pipeline_g_mmseq_a_diamond_c_mcl_hash_opt2(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_core_analysis':    
        run_main_pipeline_g_mmseq_a_diamond_c_mcl_core_analysis(args)
    if args.mode=='g_faiss_a_faiss_c_mcl':
        run_main_pipeline_faiss_mmseq_mcl(args)
    if args.mode=='g_mmseq_a_ems_c_mcl':
        run_main_pipeline_ems_mmseq_mcl(args)
def run_add(args):
    if args.mode=='g_mmseq_a_diamond_c_mcl_core_analysis':
        run_add_sample_pipeline_core_analysis(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_core_analysis2':
        run_add_sample_pipeline_core_analysis2(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_rand':
        run_add_sample_pipeline_rand(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_hash':
        run_add_sample_pipeline_hash(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_hash_opt':
        run_add_sample_pipeline_hash_opt(args)
    if args.mode=='g_mmseq_a_diamond_c_mcl_hash_opt2':
        run_add_sample_pipeline_hash_opt2(args)

def run_add_sample_pipeline(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')


    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))

    old_represent_faa = os.path.join(collection_dir, 'representative.fasta')
    if not os.path.isfile(old_represent_faa):
        raise Exception(f'{old_represent_faa} does not exist')

    old_blast_result = os.path.join(collection_dir, 'blast.tsv')
    if not os.path.isfile(old_blast_result):
        raise Exception(f'{old_blast_result} does not exist')

    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')

    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    new_combined_faa = data_preparation.combine_proteins(
        out_dir=collection_dir,
        samples=new_samples)

    not_match_faa, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1 = old_represent_faa,
        database_2 = new_combined_faa,
        out_dir = temp_dir,
        threads=threads)

    not_match_represent_faa, not_match_clusters = main_pipeline.run_cd_hit(
        faa_file=not_match_faa,
        out_dir=temp_dir,
        threads=threads)

    blast_1_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = old_represent_faa,
        query_fasta = not_match_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast1'),
        evalue = evalue,
        threads=threads
        )

    blast_2_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = not_match_represent_faa,
        query_fasta = not_match_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast2'),
        evalue = evalue,
        threads=threads
        )

    combined_blast_result = add_sample_pipeline.combine_blast_results(
        blast_1=old_blast_result,
        blast_2=blast_1_result,
        blast_3=blast_2_result,
        outdir=temp_dir)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=combined_blast_result,
        #gene_annotation=gene_annotation,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)


    logger.info(f'len cd_hit_2d_clusters = {len(cd_hit_2d_clusters)} len not_match_clusters = {len(not_match_clusters)} len old_clusters = {len(old_clusters)}')
    inflated_clusters, new_clusters = add_sample_pipeline.reinflate_clusters(
        old_clusters=old_clusters,
        cd_hit_2d_clusters=cd_hit_2d_clusters,
        not_match_clusters=not_match_clusters,
        mcl_file=mcl_file
        )
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len new clusters = {len(new_clusters)}')

    # post analysis
    #new_samples.extend(old_samples)
    old_samples.extend(new_samples)
    new_samples = old_samples
    #new_samples.sort(key= lambda x:x['id'])

    split_clusters = post_analysis.split_paralogs(
        #gene_annotation_fn=gene_annotation_fn,
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)

    output.create_outputs(annotated_clusters,new_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #print(annotated_clusters)
    if args.alignment:
        samples_dir = os.path.join(collection_dir, 'samples')
        if not os.path.exists(samples_dir):
            raise Exception(f'{samples_dir} does not exist')
        post_analysis.run_gene_alignment(annotated_clusters, new_samples,old_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)

    # output for next run
    #main_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv.gz')
    #main_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv.gz')

    #Replace the main existing files by the new ones
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)

    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(new_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')

    #if os.path.exists(temp_dir):
    #    shutil.rmtree(temp_dir)
def run_add_sample_pipeline2(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        #shutil.rmtree(temp_dir)
        #os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    group_seq_dir=os.path.join(collection_dir, 'groups')
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'annotated_clusters.json'), 'r'))
    ref_clusters = json.load(open('gene_families/gene_family_clusters.json','r'))
    ref_db='gene_families/gene_families_diamond_db.dmnd'
    ref_cdb='gene_families/db'
    old_concensus_faa = os.path.join(collection_dir, 'consensus.fasta')
    if not os.path.isfile(old_concensus_faa):
        raise Exception(f'{old_concensus_faa} does not exist')
    
    old_consensus_db=None
    if args.blast=='diamond':
        old_consensus_db=os.path.join(collection_dir, 'consenus_diamond_db.dmnd')
    if args.blast=='mmseq':
        old_consensus_db=os.path.join(collection_dir, 'consenus_mmseq_db')
    
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')

    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=collection_dir,
        samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    mmseq_represent_fasta, mmseq_groups = main_pipeline.run_mmseq_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        group_dir=group_seq_dir,
        threads=threads)
    not_match_faa1,old_clusters,mmseq_groups = add_sample_pipeline.extend_by_match_seqs_to_ref(
        seqs_file = mmseq_represent_fasta,
        groups=mmseq_groups,
        old_clusters=old_clusters,
        out_dir = temp_dir,
        refdb=ref_db,
        refcdb=ref_cdb,
        ref_clusters=ref_clusters,
        threads=threads)
    not_match_faa2,mmseq_groups = add_sample_pipeline.match_new_sequence_to_oldcluster(
        new_seqs_file = not_match_faa1,
        old_clusters=old_clusters,
        out_dir = temp_dir,
        groups=mmseq_groups,
        consensusdb=old_consensus_db,
        threads=threads)
    
    # not_match_faa2,old_clusters,mmseq_groups = add_sample_pipeline.extend_by_match_seqs_to_ref(
    #     seqs_file = not_match_faa1,
    #     groups=mmseq_groups,
    #     old_clusters=old_clusters,
    #     out_dir = temp_dir,
    #     refdb=ref_db,
    #     ref_clusters=ref_clusters,
    #     threads=threads)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(not_match_faa2,collection_dir)
    # not_match_represent_faa, not_match_clusters = main_pipeline.run_cd_hit(
    #     faa_file=not_match_faa,
    #     out_dir=temp_dir,
    #     threads=threads)
    # group_represent_fasta2, mmseq_groups2 = main_pipeline.run_mmseq_with_map(
    #     faa_file=combined_faa_file,
    #     map_file=combined_faa_map,
    #     out_dir=temp_dir,
    #     threads=threads)
    # logger.info(f'len mmseq_clusters 2 = {len(mmseq_groups2)}') 
    blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = not_match_faa2,
        query_fasta = not_match_faa2,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters(
        groups=mmseq_groups,
        mcl_file=mcl_file)
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)
    # post analysis
    #new_samples.extend(old_samples)
    annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,old_clusters)
    old_samples.extend(new_samples)
    #new_samples = old_samples
    #new_samples.sort(key= lambda x:x['id'])

    
    json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
   
    output.create_outputs(annotated_clusters,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #print(annotated_clusters)
    if args.alignment:
        samples_dir = os.path.join(collection_dir, 'samples')
        if not os.path.exists(samples_dir):
            raise Exception(f'{samples_dir} does not exist')
        #post_analysis.run_gene_alignment(annotated_clusters, new_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)
        post_analysis.run_gene_alignment(annotated_clusters, old_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa,add=True)
        post_analysis.create_poa_protein_consensus(annotated_clusters,collection_dir,threads=threads)
        post_analysis.make_protein_consensus_db(annotated_clusters,collection_dir,threads=threads)
    # output for next run
    #main_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv.gz')
    #main_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv.gz')

    #Replace the main existing files by the new ones
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)

    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_nearest_group(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        #shutil.rmtree(temp_dir)
        #os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')


    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'annotated_clusters.json'), 'r'))
    group_seq_dir=os.path.join(collection_dir, 'groups')
    cluster_dir=os.path.join(collection_dir, 'clusters')
    ref_clusters = json.load(open('gene_families/gene_family_splited_clusters.json','r'))
    ref_db='gene_families/gene_families_diamond_db.dmnd'
    ref_cdb='gene_families/db'
    group_seq_dir=os.path.join(collection_dir, 'groups')
    
    old_concensus_faa = os.path.join(collection_dir, 'consensus.fasta')
    if not os.path.isfile(old_concensus_faa):
        raise Exception(f'{old_concensus_faa} does not exist')
    
    old_consensus_db=None
    if args.blast=='diamond':
        old_consensus_db=os.path.join(collection_dir, 'consenus_diamond_db.dmnd')
    if args.blast=='mmseq':
        old_consensus_db=os.path.join(collection_dir, 'consenus_mmseq_db')
    
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')

    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=collection_dir,
        samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    mmseq_represent_fasta, mmseq_groups = main_pipeline.run_mmseq_with_map(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,
        group_dir=group_seq_dir,
        threads=threads)
    not_match_faa1,old_clusters,mmseq_groups = add_sample_pipeline.extend_by_match_seqs_to_ref(
        seqs_file = mmseq_represent_fasta,
        groups=mmseq_groups,
        old_clusters=old_clusters,
        out_dir = temp_dir,
        refdb=ref_db,
        refcdb=ref_cdb,
        ref_clusters=ref_clusters,
        threads=threads)
    not_match_faa2,mmseq_groups = add_sample_pipeline.match_new_sequence_to_oldcluster(
        new_seqs_file = not_match_faa1,
        old_clusters=old_clusters,
        out_dir = temp_dir,
        groups=mmseq_groups,
        group_dir=group_seq_dir,
        cluster_dir=cluster_dir,
        consensusdb=old_consensus_db,
        threads=threads)
    
    # not_match_faa2,old_clusters,mmseq_groups = add_sample_pipeline.extend_by_match_seqs_to_ref(
    #     seqs_file = not_match_faa1,
    #     groups=mmseq_groups,
    #     old_clusters=old_clusters,
    #     out_dir = temp_dir,
    #     refdb=ref_db,
    #     ref_clusters=ref_clusters,
    #     threads=threads)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(not_match_faa2,collection_dir)
    # not_match_represent_faa, not_match_clusters = main_pipeline.run_cd_hit(
    #     faa_file=not_match_faa,
    #     out_dir=temp_dir,
    #     threads=threads)
    # group_represent_fasta2, mmseq_groups2 = main_pipeline.run_mmseq_with_map(
    #     faa_file=combined_faa_file,
    #     map_file=combined_faa_map,
    #     out_dir=temp_dir,
    #     threads=threads)
    # logger.info(f'len mmseq_clusters 2 = {len(mmseq_groups2)}') 
    blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = not_match_faa2,
        query_fasta = not_match_faa2,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)

    inflated_clusters, clusters = main_pipeline.reinflate_clusters_by_groups(
        groups=mmseq_groups,
        mcl_file=mcl_file)
    logger.info(f'len inflated_clusters = {len(inflated_clusters)} len clusters = {len(clusters)}')


    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len split_clusters = {len(split_clusters)}')

    annotated_clusters = post_analysis.annotate_cluster_main(
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn)
    # post analysis
    #new_samples.extend(old_samples)
    annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,old_clusters)
    old_samples.extend(new_samples)
    #new_samples = old_samples
    #new_samples.sort(key= lambda x:x['id'])

    
    json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
   
    output.create_outputs(annotated_clusters,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #print(annotated_clusters)
    samples_dir = os.path.join(collection_dir, 'samples')
    if not os.path.exists(samples_dir):
        raise Exception(f'{samples_dir} does not exist')
    #post_analysis.run_gene_alignment(annotated_clusters, new_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)
    post_analysis.run_gene_alignment(annotated_clusters, old_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa,add=True)
    post_analysis.create_poa_protein_consensus(annotated_clusters,collection_dir,gene_families_dir='gene_families',threads=threads)
    post_analysis.make_protein_consensus_db(annotated_clusters,collection_dir,threads=threads)
    post_analysis.create_protein_db_for_clusters(annotated_clusters,collection_dir,group_seq_dir,threads=threads)
    
    # output for next run
    #main_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv.gz')
    #main_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv.gz')

    #Replace the main existing files by the new ones
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)

    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
    #if os.path.exists(temp_dir):
    #    shutil.rmtree(temp_dir)

def build_reference_gene_family_db(args):
    starttime = datetime.now()
    ref_db_dir='gene_families'
    if os.path.exists(ref_db_dir):
        shutil.rmtree(ref_db_dir)
    else:
        os.mkdir(ref_db_dir)
        os.mkdir(ref_db_dir+"/sequences")
        os.mkdir(ref_db_dir+"/db")
    mmseq_represent_fasta= os.path.join(ref_db_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(ref_db_dir, 'mmseq_cluster.tsv')
    cmd = f'mmseqs easy-linclust {args.input} {ref_db_dir}/mmseq {ref_db_dir}/tmp --min-seq-id {args.identity} -c {args.coverage} --threads {args.threads} > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseq')        

    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq2 with {args.identity} identity part 1 -- time taken {elapsed}')

    groups = {}
   
    count = 0
    
    represent_corrected_fasta = os.path.join(ref_db_dir, 'mmseq.fasta')
    map_id_des={}
    unique_gene_name=set()
    num_dup=0
    with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                str_name=line[1:].strip()
                
                geneid=str_name.split(' ')[0].split('|')[1]
                
                des=str_name[str_name.find(' ') + 1:]
                map_id_des[geneid]={'des':des}
                pattern = r'GN=([^ ]+)'


                match = re.search(pattern, des)


                if match:
                 
                    genename = match.group(1)
                    
                else:
                    genename=geneid
                # if genename in unique_gene_name:
                #     suffix=1
                #     new_gene_name=genename+'_'+str(suffix)
                #     num_dup=num_dup+1
                #     while new_gene_name in unique_gene_name:
                #         suffix=suffix+1
                #         new_gene_name=genename+'_'+str(suffix)
                #     map_id_des[geneid]['genename']=   new_gene_name
                #     unique_gene_name.add(new_gene_name)
                #     ofh.write(f'>{new_gene_name}\n')
                # else:
                #     map_id_des[geneid]['genename']= genename
                #     unique_gene_name.add(genename)
                #     ofh.write(f'>{genename}\n')
                if genename in unique_gene_name:
                    num_dup=num_dup+1
                else:
                    unique_gene_name.add(genename)
                map_id_des[geneid]['gene_name']= genename
                unique_gene_name.add(genename)
                ofh.write(f'>{geneid}\n')       
            else:
                ofh.write(line)  
        
    for ref_seq in SeqIO.parse(open(represent_corrected_fasta),'fasta'):    
        SeqIO.write(ref_seq, ref_db_dir+"/sequences/"+ref_seq.id+".fasta", "fasta")
    elapsed = datetime.now() - starttime
    #logging.info(f'Run mmseq with {args.identity} identity part 1 -- time taken {elapsed}')
    c_cursor='0'
    cluster_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=rep_name
                
                groups[c_cursor] = {'gene_id':[]} 
                groups[c_cursor]['representative'] = ref_db_dir+"/sequences/"+member+".fasta"
                groups[c_cursor]['gene_name']=map_id_des[c_cursor]['gene_name']
                cluster_count=cluster_count+1
                
            
            groups[c_cursor]['gene_id'].append(member)                 

           
    logging.info(f'###{args.identity} {args.coverage} {cluster_count} {num_dup}')
    diamond_result=main_pipeline.pairwise_alignment_diamond(represent_corrected_fasta,represent_corrected_fasta,ref_db_dir+'/temp')          
    mcl_file=main_pipeline.cluster_with_mcl(diamond_result,ref_db_dir+'/temp',inflation=2)
    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    clusters=[]
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            cluster={'id':genes[0],'groups':[]}
            for gene in genes:
                #inflated_genes.append(gene)
                
                inflated_genes.extend(groups[gene]['gene_id'])
                cluster['groups'].append(gene)
            clusters.append(cluster) 
            inflated_clusters.append(inflated_genes)
    dict_gene={}
    with open(args.input) as fh:
        for seq in SeqIO.parse(fh, 'fasta'):
            geneid=seq.id.split(' ')[0].split('|')[1]
            dict_gene[geneid]=SeqRecord(seq.seq, id = geneid, description = '')
    
    map_name={}
    big_clusters={}
    for c in clusters:
        big_clusters[c['id']]={'groups':[]}

        big_clusters[c['id']]['gene_name']=map_id_des[c['id']]['gene_name']
        big_clusters[c['id']]['description']=map_id_des[c['id']]['des']
        #clusters[k]['gene_name']=map_id_des[c[0]]['gene_name']
        #clusters[k]['description']=map_id_des[c[0]]['des']
        for g in c['groups']:

            big_clusters[c['id']]['groups'].append(groups[g])
       
    #split clusters by gene name
    splited_clusters={}
    for c in big_clusters.keys():
        #count gene name in clusters:
        unique_gene_in_cluster={}
        set_names=set()
        for g in big_clusters[c]['groups']:
        
            if g['gene_name'] not in unique_gene_in_cluster.keys():
                unique_gene_in_cluster[g['gene_name']]=1
            else:
                unique_gene_in_cluster[g['gene_name']]=1 + unique_gene_in_cluster[g['gene_name']]
        if len(unique_gene_in_cluster.keys())>1:
            #neeed to split
            for gene_name in unique_gene_in_cluster.keys():
                list_groups=[]
                for g in big_clusters[c]['groups']:
                    if g['gene_name']==gene_name:
                        list_groups.append(g)
                newid=list_groups[0]['gene_id'][0]
                splited_clusters[newid]={}
                splited_clusters[newid]['groups']=list_groups
                #for g in list_groups:
                #    splited_clusters[newid]['groups'].append(g)
                splited_clusters[newid]['gene_name']=gene_name
                splited_clusters[newid]['description']=map_id_des[newid]['des']
        else:
           splited_clusters[c]=big_clusters[c]
    pool = multiprocessing.Pool(processes=args.threads)
    results = []
    for c in splited_clusters.keys():
        #make fasta for each cluster
        gene_seq_consensus_file=ref_db_dir+"/temp/"+c+".cons.faa"
        gene_seq_file=ref_db_dir+"/temp/"+c+".faa"
        with open(gene_seq_file, 'wt') as fh:
            for g in splited_clusters[c]['groups']:
                for gene in g['gene_id']:
            
                    SeqIO.write(dict_gene[gene], fh, 'fasta')
        cmd = f"abpoa -c -t BLOSUM62.mtx {gene_seq_file} > {gene_seq_consensus_file} && sed -i 's/^>Consensus_sequence/>{c}/' {gene_seq_consensus_file}"
        results.append(pool.apply_async(run_command,(cmd, None)))
        gene_diamonddb_file=ref_db_dir+"/db/"+c+".db"
        cmd2=f'./diamond makedb --in {gene_seq_file} -d {gene_diamonddb_file} -p 1 --quiet'
        results.append(pool.apply_async(run_command,(cmd2, None)))
    pool.close()
    pool.join()

    represent_clusters_file=os.path.join(ref_db_dir, 'ref_consensus.fasta')
    with open(represent_clusters_file, 'wt') as fh:
        for c in clusters:
            consensus_gene_file=ref_db_dir+"/temp/"+c['id']+".cons.faa"
            with open(consensus_gene_file) as fi:
                for seq in SeqIO.parse(fi, 'fasta'):
                    SeqIO.write(seq, fh, 'fasta')
    diamond_db = os.path.join(ref_db_dir, 'gene_families_diamond_db')

    cmd = f'./diamond makedb --in {represent_clusters_file} -d {diamond_db} -p {args.threads} --quiet'

    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond makedb for consensus file')

    diamond_group_db = os.path.join(ref_db_dir, 'gene_groups_diamond_db')

    cmd = f'./diamond makedb --in {represent_corrected_fasta} -d {diamond_group_db} -p {args.threads} --quiet'

    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond makedb for groups')
    json.dump(groups, open(os.path.join(ref_db_dir, 'gene_family_groups.json'), 'w'), indent=4, sort_keys=True)
    
    
    json.dump(big_clusters, open(os.path.join(ref_db_dir, 'gene_family_clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(splited_clusters, open(os.path.join(ref_db_dir, 'gene_family_splited_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    elapsed = datetime.now() - starttime
    #logging.info(f'Run make gene family db with {args.identity} identity -- time taken {str(elapsed)}')
def build_reference_gene_family_db_hmmer(args):
    starttime = datetime.now()
    ref_db_dir='hmmer_gene_families'
    if os.path.exists(ref_db_dir):
        shutil.rmtree(ref_db_dir)
    else:
        os.mkdir(ref_db_dir)
        os.mkdir(ref_db_dir+"/sequences")
        os.mkdir(ref_db_dir+"/db")
    mmseq_represent_fasta= os.path.join(ref_db_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(ref_db_dir, 'mmseq_cluster.tsv')
    cmd = f'mmseqs easy-linclust {args.input} {ref_db_dir}/mmseq {ref_db_dir}/tmp --min-seq-id {args.identity} -c {args.coverage} --threads {args.threads} > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseq')        

    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq2 with {args.identity} identity part 1 -- time taken {elapsed}')

    groups = {}
   
    count = 0
    
    represent_corrected_fasta = os.path.join(ref_db_dir, 'mmseq.fasta')
    map_id_des={}
    unique_gene_name=set()
    num_dup=0
    with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                str_name=line[1:].strip()
                
                geneid=str_name.split(' ')[0].split('|')[1]
                
                des=str_name[str_name.find(' ') + 1:]
                map_id_des[geneid]={'des':des}
                pattern = r'GN=([^ ]+)'


                match = re.search(pattern, des)


                if match:
                 
                    genename = match.group(1)
                    
                else:
                    genename=geneid
                # if genename in unique_gene_name:
                #     suffix=1
                #     new_gene_name=genename+'_'+str(suffix)
                #     num_dup=num_dup+1
                #     while new_gene_name in unique_gene_name:
                #         suffix=suffix+1
                #         new_gene_name=genename+'_'+str(suffix)
                #     map_id_des[geneid]['genename']=   new_gene_name
                #     unique_gene_name.add(new_gene_name)
                #     ofh.write(f'>{new_gene_name}\n')
                # else:
                #     map_id_des[geneid]['genename']= genename
                #     unique_gene_name.add(genename)
                #     ofh.write(f'>{genename}\n')
                if genename in unique_gene_name:
                    num_dup=num_dup+1
                else:
                    unique_gene_name.add(genename)
                map_id_des[geneid]['gene_name']= genename
                unique_gene_name.add(genename)
                ofh.write(f'>{geneid}\n')       
            else:
                ofh.write(line)  
        
    for ref_seq in SeqIO.parse(open(represent_corrected_fasta),'fasta'):    
        SeqIO.write(ref_seq, ref_db_dir+"/sequences/"+ref_seq.id+".fasta", "fasta")
    elapsed = datetime.now() - starttime
    #logging.info(f'Run mmseq with {args.identity} identity part 1 -- time taken {elapsed}')
    c_cursor='0'
    cluster_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=rep_name
                
                groups[c_cursor] = {'gene_id':[]} 
                groups[c_cursor]['representative'] = ref_db_dir+"/sequences/"+member+".fasta"
                groups[c_cursor]['gene_name']=map_id_des[c_cursor]['gene_name']
                cluster_count=cluster_count+1
                
            
            groups[c_cursor]['gene_id'].append(member)                 

           
    logging.info(f'###{args.identity} {args.coverage} {cluster_count} {num_dup}')
    diamond_result=main_pipeline.pairwise_alignment_diamond(represent_corrected_fasta,represent_corrected_fasta,ref_db_dir+'/temp')          
    mcl_file=main_pipeline.cluster_with_mcl(diamond_result,ref_db_dir+'/temp',inflation=2)
    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    clusters=[]
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            cluster={'id':genes[0],'groups':[]}
            for gene in genes:
                #inflated_genes.append(gene)
                
                inflated_genes.extend(groups[gene]['gene_id'])
                cluster['groups'].append(gene)
            clusters.append(cluster) 
            inflated_clusters.append(inflated_genes)
    dict_gene={}
    with open(args.input) as fh:
        for seq in SeqIO.parse(fh, 'fasta'):
            geneid=seq.id.split(' ')[0].split('|')[1]
            dict_gene[geneid]=SeqRecord(seq.seq, id = geneid, description = '')
    
    map_name={}
    big_clusters={}
    for c in clusters:
        big_clusters[c['id']]={'groups':[]}

        big_clusters[c['id']]['gene_name']=map_id_des[c['id']]['gene_name']
        big_clusters[c['id']]['description']=map_id_des[c['id']]['des']
        #clusters[k]['gene_name']=map_id_des[c[0]]['gene_name']
        #clusters[k]['description']=map_id_des[c[0]]['des']
        for g in c['groups']:

            big_clusters[c['id']]['groups'].append(groups[g])
       
    #split clusters by gene name
    
    pool = multiprocessing.Pool(processes=args.threads)
    results = []
    for c in big_clusters.keys():
        #make fasta for each cluster
        gene_seq_aln_file=ref_db_dir+"/temp/"+c+".aln"
        gene_seq_file=ref_db_dir+"/temp/"+c+".faa"
        with open(gene_seq_file, 'wt') as fh:
            for g in big_clusters[c]['groups']:
                for gene in g['gene_id']:
            
                    SeqIO.write(dict_gene[gene], fh, 'fasta')
        cmd = f"abpoa {gene_seq_file} -r1 > {gene_seq_aln_file}"
        results.append(pool.apply_async(run_command,(cmd, None)))
        gene_hmmer_file=ref_db_dir+"/db/"+c+".hmm"
        cmd2=f'hmmbuild --cpu 1 --amino {gene_hmmer_file} {gene_seq_aln_file}'
        results.append(pool.apply_async(run_command,(cmd2, None)))
    pool.close()
    pool.join()

    hmmer_clusters_file=os.path.join(ref_db_dir, 'ref.hmm')
    for c in clusters:
        gene_hmmer_file=ref_db_dir+"/db/"+c+".hmm"
        cmd= f"cat {hmmer_clusters_file} {gene_hmmer_file} > {hmmer_clusters_file}"
        ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond makedb for consensus file')

    

    #ret = os.system(cmd)
    

    
    json.dump(groups, open(os.path.join(ref_db_dir, 'gene_family_groups.json'), 'w'), indent=4, sort_keys=True)
    
    
    json.dump(big_clusters, open(os.path.join(ref_db_dir, 'gene_family_clusters.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(splited_clusters, open(os.path.join(ref_db_dir, 'gene_family_splited_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    elapsed = datetime.now() - starttime    
def run_add_sample_pipeline_rand(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'annotated_clusters.json'), 'r'))
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    old_unique_groups=json.load(open(os.path.join(collection_dir, 'unique_groups.json'), 'r'))
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')

    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=collection_dir,
        samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    new_unique_seds_fasta, new_unique_groups = main_pipeline.run_mmseq_with_map_unique_seqs(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,        
        threads=threads)
    #json.dump(new_unique_groups, open(os.path.join(collection_dir, 'new_unique_groups.json'), 'w'), indent=4, sort_keys=True)
    concat_unique_seq_fasta=concat2fasta(old_unique_seqs,new_unique_seds_fasta,os.path.join(temp_dir,"concat_unique_seqs.fasta"))
    combined_unique_seqs_fasta, combined_unique_group=main_pipeline.run_mmseq_unique_seqs(
        faa_file=concat_unique_seq_fasta,
        out_dir=temp_dir,        
        threads=threads)
    #json.dump(combined_unique_group, open(os.path.join(collection_dir, 'combined_unique_group.json'), 'w'), indent=4, sort_keys=True)
    flat_unique_seq=add_sample_pipeline.flat_combined_unique_seq(new_unique_groups,old_unique_groups,combined_unique_group)
    json.dump(flat_unique_seq, open(os.path.join(collection_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    groups_representative_fasta, similar_groups = main_pipeline.run_mmseq_with_map_similar_seqs(
        faa_file=combined_unique_seqs_fasta,
        
        out_dir=temp_dir,      
        threads=threads)
    blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = groups_representative_fasta,
        query_fasta = groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)
    #json.dump(similar_groups, open(os.path.join(collection_dir, 'similar_groups.json'), 'w'), indent=4, sort_keys=True)
    
    inflated_clusters, clusters, unique_seq_by_similar_group = add_sample_pipeline.reinflate_clusters_with_hirachical_group(
        old_unique_groups=old_unique_groups,
        new_unique_groups=new_unique_groups,
        combined_unique_group=combined_unique_group,
        similar_groups=similar_groups,
        mcl_file=mcl_file)
    logger.info(f'len new inflated_clusters = {len(inflated_clusters)} len groups = {len(clusters)}')
    #merge clusters
    
    #json.dump(inflated_clusters, open(os.path.join(collection_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len new split_clusters = {len(split_clusters)}')
    json.dump(split_clusters, open(os.path.join(collection_dir, 'split_clusters.json'), 'w'), indent=4, sort_keys=True)

    annotated_clusters = post_analysis.annotate_cluster_progressive_unique_seq(
        unique_seq_by_similar_group=unique_seq_by_similar_group,
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn,
        out_dir = temp_dir
        )
    json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    logger.info(f'samples after extend = {len(old_samples)}')
    output.create_outputs1(annotated_clusters,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #print(annotated_clusters)
    if args.alignment:
        samples_dir = os.path.join(collection_dir, 'samples')
        if not os.path.exists(samples_dir):
            raise Exception(f'{samples_dir} does not exist')
        #post_analysis.run_gene_alignment(annotated_clusters, new_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)
        post_analysis.run_gene_alignment(annotated_clusters, old_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa,add=True)
        post_analysis.create_poa_protein_consensus(annotated_clusters,collection_dir,threads=threads)
        post_analysis.make_protein_consensus_db(annotated_clusters,collection_dir,threads=threads)
    # output for next run
    #main_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv.gz')
    #main_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv.gz')

    #Replace the main existing files by the new ones
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
    shutil.copy(combined_unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_cc(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)
    #tracemalloc.start()
    mem_usage = mem_report(0, "begin")
    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    total_used_mem=0
    timing_log= os.path.join(collection_dir, 'cmds.log')
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    mem_usage = mem_report(mem_usage, "after old_samples")
    #old_clusters = json.load(open(os.path.join(collection_dir, 'annotated_clusters.json'), 'r'))
    old_clusters =os.path.join(collection_dir, 'clusters.json')
    mem_usage = mem_report(mem_usage, "after old_clusters")
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    #old_unique_groups=json.load(open(os.path.join(collection_dir, 'unique_groups.json'), 'r'))
    old_unique_groups=os.path.join(collection_dir, 'unique_groups.json')
    mem_usage = mem_report(mem_usage, "after old_unique_groups")
    old_representative_clusters=os.path.join(collection_dir, 'representative_clusters.fasta')
    #print("mem after reload")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after loading")
    #check size
    total_used_mem=total_used_mem+sys.getsizeof(old_samples)
    logging.info("size of old_samples: "+str(sys.getsizeof(old_samples))+" , total is "+str(total_used_mem))
    
    #total_used_mem=total_used_mem+sys.getsizeof(old_unique_groups)
    #logging.info("size of old_unique_groups: "+str(sys.getsizeof(old_unique_groups)/1E9)+" GB , total is "+str(total_used_mem))
    #total_used_mem=total_used_mem+sys.getsizeof(old_clusters)
    #logging.info("size of old_clusters: "+str(sys.getsizeof(old_unique_groups)/1E9)+" GB , total is "+str(total_used_mem))
    
   # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    total_used_mem=total_used_mem+sys.getsizeof(new_samples)
    logging.info("size of new_samples: "+str(sys.getsizeof(new_samples))+" , total is "+str(total_used_mem))
    
    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa = data_preparation.combine_proteins(
        out_dir=collection_dir,
        samples=new_samples)
    #print("mem after combine")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after combine")
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    new_unique_seds_fasta, new_unique_groups = main_pipeline.run_mmseq_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads,
        timing_log=timing_log
        )
    json.dump(new_unique_groups, open(os.path.join(collection_dir, 'new_unique_groups1.json'), 'w'), indent=4, sort_keys=True)
    total_used_mem=total_used_mem+sys.getsizeof(new_unique_groups)
    logging.info("size of new_samples: "+str(sys.getsizeof(new_unique_groups))+" , total is "+str(total_used_mem))
    #print("mem after run_mmseq_unique_seqs")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after run_mmseq_unique_seqs")
    _old_clusters= json.load(open(old_clusters, 'r'))
    count_num_seq_in_cluster=0
    for c in _old_clusters.keys():
        list_useq=read_array(_old_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in old cluster {count_num_seq_in_cluster} 1")
        
    old_clusters, old_unique_groups,new_unique_groups,new_unique_seds_fasta=main_pipeline.identical_matching_old_clusters(
        old_clusters_file=old_clusters,  
        old_unique_groups_file=old_unique_groups,
        old_unique_sequence=old_unique_seqs,
        new_unique_groups=new_unique_groups,
        new_unique_seqs=new_unique_seds_fasta,
        out_dir=temp_dir,
        root_dir=collection_dir,
        threads=threads,
        timing_log=timing_log)
    _old_clusters= json.load(open(old_clusters, 'r'))
    count_num_seq_in_cluster=0
    for c in _old_clusters.keys():
        list_useq=read_array(_old_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in old cluster {count_num_seq_in_cluster} 2")
    #json.dump(new_unique_groups, open(os.path.join(collection_dir, 'new_unique_groups2.json'), 'w'), indent=4, sort_keys=True)
    #print("mem after identical_matching_old_clusters")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after identical_matching_old_clusters")
    unmatch_combined_faa_file,updated_clusters=main_pipeline.simple_match_new_seq_to_old_clusters(
        old_clusters_file=old_clusters,
        representative_clusters_file=old_representative_clusters,
        new_sequences_file=new_unique_seds_fasta,
        old_sequences_file=old_unique_seqs,
        out_dir=temp_dir,
        threads=threads
    )
    # unmatch_combined_faa_file,updated_clusters=main_pipeline.hc_match_new_seq_to_old_clusters(
    #     old_clusters_file=old_clusters,
    #     representative_clusters_file=old_representative_clusters,
    #     new_sequences_file=new_unique_seds_fasta,
    #     old_sequences_file=old_unique_seqs,
    #     out_dir=temp_dir,
    #     threads=threads
    # )
    total_used_mem=total_used_mem+sys.getsizeof(updated_clusters)
    logging.info("size of updated_clusters: "+str(sys.getsizeof(updated_clusters)/1E9)+" GB , total is "+str(total_used_mem))
    #print("mem after simple_match_new_seq_to_old_clusters")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after simple_match_new_seq_to_old_clusters")
    _old_clusters= json.load(open(old_clusters, 'r'))
    count_num_seq_in_cluster=0
    for c in _old_clusters.keys():
        list_useq=read_array(_old_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in old cluster {count_num_seq_in_cluster} 3")
    new_inflated_clusters, similar_groups, representative_similar_groups,new_representative_clusters=main_pipeline.clustering(
        input_fasta_file=unmatch_combined_faa_file,
        out_dir=collection_dir,
        threads=threads,
        timing_log=timing_log
        )
    total_used_mem=total_used_mem+sys.getsizeof(new_inflated_clusters)
    logging.info("size of new_inflated_clusters: "+str(sys.getsizeof(new_inflated_clusters))+" , total is "+str(total_used_mem))
    #print("mem after clustering")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after clustering")
    
    new_annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=new_inflated_clusters,
        gene_annotation_fn=gene_annotation_fn)
    total_used_mem=total_used_mem+sys.getsizeof(new_annotated_clusters)
    logging.info("size of new_annotated_clusters: "+str(sys.getsizeof(new_annotated_clusters))+" , total is "+str(total_used_mem))
    #print("mem after annotate_cluster")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after annotate_cluster")

    merge_clusters=post_analysis.merge_new_cluster_to_old_clusters(
        old_annotated_clusters_file=updated_clusters,
        new_annotated_clusters=new_annotated_clusters,
        root_dir=collection_dir
        )
    _old_clusters= json.load(open(old_clusters, 'r'))
    count_num_seq_in_cluster=0
    for c in _old_clusters.keys():
        list_useq=read_array(_old_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in old cluster {count_num_seq_in_cluster} 4")
    #print("mem after merge_new_cluster_to_old_clusters")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after merge_new_cluster_to_old_clusters")
    
    merge_clusters=post_analysis.expandClusterMembers(
        merge_clusters,
        new_unique_groups)
    total_used_mem=total_used_mem+sys.getsizeof(merge_clusters)
    logging.info("size of merge_clusters: "+str(sys.getsizeof(merge_clusters)/1E9)+" GB, total is "+str(total_used_mem))
    #print("mem after expandClusterMembers")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after expandClusterMembers")
    _old_clusters= json.load(open(old_clusters, 'r'))
    count_num_seq_in_cluster=0
    for c in _old_clusters.keys():
        list_useq=read_array(_old_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in old cluster {count_num_seq_in_cluster} 5")
    count_num_seq_in_cluster=0
    _merge_clusters= json.load(open(merge_clusters, 'r'))
    for c in _merge_clusters.keys():
        list_useq=read_array(_merge_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in merge cluster {count_num_seq_in_cluster} 0")
    annotated_clusters=main_pipeline.refine_updated_clusters(
        updated_clusters_file=merge_clusters,
        new_annotated_clusters=new_annotated_clusters,
        old_representatve=old_representative_clusters,
        new_representative=new_representative_clusters,
        gene_annotation_fn=gene_annotation_fn,
        unique_sequences=old_unique_seqs,
        out_dir=temp_dir,
        root_dir=collection_dir,
        threads=threads)
    total_used_mem=total_used_mem+sys.getsizeof(annotated_clusters)
    logging.info("size of annotated_clusters: "+str(sys.getsizeof(annotated_clusters)/1E9)+" GB , total is "+str(total_used_mem))
    #print("mem after refine_updated_clusters")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after refine_updated_clusters")
    _old_clusters= json.load(open(old_clusters, 'r'))
    count_num_seq_in_cluster=0
    for c in _old_clusters.keys():
        list_useq=read_array(_old_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in old cluster {count_num_seq_in_cluster} 6")
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    logger.info(f'samples after extend = {len(old_samples)}')
    output.create_outputs(annotated_clusters,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #print(annotated_clusters)
    #print("mem after create_outputs")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after create_outputs")
    # if args.alignment:
    #     samples_dir = os.path.join(collection_dir, 'samples')
    #     if not os.path.exists(samples_dir):
    #         raise Exception(f'{samples_dir} does not exist')
    #     #post_analysis.run_gene_alignment(annotated_clusters, new_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)
    #     post_analysis.run_gene_alignment(annotated_clusters, old_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa,add=True)
    #     post_analysis.create_poa_protein_consensus(annotated_clusters,collection_dir,threads=threads)
    #     post_analysis.make_protein_consensus_db(annotated_clusters,collection_dir,threads=threads)
    # # output for next run
    #main_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv.gz')
    #main_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv.gz')
    represent_clusters_file=main_pipeline.make_representative_clusters(annotated_clusters,old_unique_seqs,collection_dir)
    #print("mem after make_representative_clusters")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after make_representative_clusters")
    #Replace the main existing files by the new ones
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
    #shutil.copy(combined_unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(old_unique_groups, open(os.path.join(collection_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)

    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    #print("mem after all")
    #print(tracemalloc.get_traced_memory())
    #tracemalloc.stop()
    mem_usage = mem_report(mem_usage, "after all")
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_core_analysis(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)
    #tracemalloc.start()
    mem_usage = mem_report(0, "begin")
    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_presentation_tab = os.path.join(collection_dir, 'gene_presence_absence.Rtab')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    total_used_mem=0
    timing_log= os.path.join(collection_dir, 'cmds.log')
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    mem_usage = mem_report(mem_usage, "after old_samples")
    #old_clusters = json.load(open(os.path.join(collection_dir, 'annotated_clusters.json'), 'r'))
    old_clusters =os.path.join(collection_dir, 'clusters.json')
    mem_usage = mem_report(mem_usage, "after old_clusters")
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    #old_unique_groups=json.load(open(os.path.join(collection_dir, 'unique_groups.json'), 'r'))
    old_unique_groups=os.path.join(collection_dir, 'unique_groups.json')
    mem_usage = mem_report(mem_usage, "after old_unique_groups")
    old_representative_clusters=os.path.join(collection_dir, 'representative_clusters.fasta')
    #print("mem after reload")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after loading")
    #check size
    total_used_mem=total_used_mem+sys.getsizeof(old_samples)
    logging.info("size of old_samples: "+str(sys.getsizeof(old_samples))+" , total is "+str(total_used_mem))
    
    #total_used_mem=total_used_mem+sys.getsizeof(old_unique_groups)
    #logging.info("size of old_unique_groups: "+str(sys.getsizeof(old_unique_groups)/1E9)+" GB , total is "+str(total_used_mem))
    #total_used_mem=total_used_mem+sys.getsizeof(old_clusters)
    #logging.info("size of old_clusters: "+str(sys.getsizeof(old_unique_groups)/1E9)+" GB , total is "+str(total_used_mem))
    
   # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    total_used_mem=total_used_mem+sys.getsizeof(new_samples)
    logging.info("size of new_samples: "+str(sys.getsizeof(new_samples))+" , total is "+str(total_used_mem))
    
    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa = data_preparation.combine_proteins(
        out_dir=collection_dir,
        samples=new_samples)
    #print("mem after combine")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after combine")
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    new_unique_seds_fasta, new_unique_groups = main_pipeline.run_mmseq_unique_seqs(
        faa_file=combined_faa,
      
        out_dir=temp_dir,      
        threads=threads,
        timing_log=timing_log
        )
    json.dump(new_unique_groups, open(os.path.join(collection_dir, 'new_unique_groups1.json'), 'w'), indent=4, sort_keys=True)
    total_used_mem=total_used_mem+sys.getsizeof(new_unique_groups)
    logging.info("size of new_samples: "+str(sys.getsizeof(new_unique_groups))+" , total is "+str(total_used_mem))
    #print("mem after run_mmseq_unique_seqs")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after run_mmseq_unique_seqs")
    _old_clusters= json.load(open(old_clusters, 'r'))
    count_num_seq_in_cluster=0
    for c in _old_clusters.keys():
        list_useq=read_array(_old_clusters[c]['gene_id'])
        count_num_seq_in_cluster=count_num_seq_in_cluster+len(list_useq)
    logging.info(f"number seq in old cluster {count_num_seq_in_cluster} 1")
        
    
    new_inflated_clusters, similar_groups, representative_similar_groups=main_pipeline.clustering(
        input_fasta_file=new_unique_seds_fasta,
        out_dir=collection_dir,
        threads=threads,
        timing_log=timing_log
        )
    total_used_mem=total_used_mem+sys.getsizeof(new_inflated_clusters)
    logging.info("size of new_inflated_clusters: "+str(sys.getsizeof(new_inflated_clusters))+" , total is "+str(total_used_mem))
    #print("mem after clustering")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after clustering")
    
    new_annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=new_inflated_clusters,
        gene_annotation_fn=gene_annotation_fn)
    total_used_mem=total_used_mem+sys.getsizeof(new_annotated_clusters)
    logging.info("size of new_annotated_clusters: "+str(sys.getsizeof(new_annotated_clusters))+" , total is "+str(total_used_mem))
    #print("mem after annotate_cluster")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after annotate_cluster")
    #new_annotated_clusters,old_annotated_clusters_file,old_unique_seqs,
    # new_unique_seqs,gene_present_tab,root_dir,temp_dir
    unique_seqs_fasta,merge_clusters=main_pipeline.merge_new_cluster_to_old_clusters_by_core_analysis2(
        new_annotated_clusters=new_annotated_clusters,
        old_annotated_clusters_file=old_clusters,
        old_unique_seqs=old_unique_seqs,
        new_unique_seqs=new_unique_seds_fasta,
        gene_present_tab=gene_presentation_tab,
        root_dir=collection_dir,
        temp_dir=temp_dir,
        threads=threads,
        timing_log=timing_log
        )
    
   
    #clusters_file, old_groups,new_groups,old_seqs_fasta,new_seqs_fasta,root_dir,temp_dir
    unique_seqs_fasta=main_pipeline.reduce_unique_seqs_and_expand_seq_ids(
        clusters_file=old_clusters,
        old_groups_file=old_unique_groups,
        new_groups=new_unique_groups,
        old_seqs_fasta=old_unique_seqs,
        new_seqs_fasta=new_unique_seds_fasta,
        root_dir=collection_dir,
        temp_dir=temp_dir,
        threads=threads,
        timing_log=timing_log
    )
   
    
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    logger.info(f'samples after extend = {len(old_samples)}')
    output.create_outputs(old_clusters,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #print(annotated_clusters)
    #print("mem after create_outputs")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after create_outputs")
    # if args.alignment:
    #     samples_dir = os.path.join(collection_dir, 'samples')
    #     if not os.path.exists(samples_dir):
    #         raise Exception(f'{samples_dir} does not exist')
    #     #post_analysis.run_gene_alignment(annotated_clusters, new_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)
    #     post_analysis.run_gene_alignment(annotated_clusters, old_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa,add=True)
    #     post_analysis.create_poa_protein_consensus(annotated_clusters,collection_dir,threads=threads)
    #     post_analysis.make_protein_consensus_db(annotated_clusters,collection_dir,threads=threads)
    # # output for next run
    #main_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv.gz')
    #main_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv.gz')
    represent_clusters_file=main_pipeline.make_representative_clusters(old_clusters,old_unique_seqs,collection_dir)
    #print("mem after make_representative_clusters")
    #print(tracemalloc.get_traced_memory())
    mem_usage = mem_report(mem_usage, "after make_representative_clusters")
    #Replace the main existing files by the new ones
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
    shutil.copy(unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(old_unique_groups, open(os.path.join(collection_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)

    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    #print("mem after all")
    #print(tracemalloc.get_traced_memory())
    #tracemalloc.stop()
    mem_usage = mem_report(mem_usage, "after all")
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_core_analysis2(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')

    
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    old_unique_groups=json.load(open(os.path.join(collection_dir, 'unique_groups.json'), 'r'))
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')

    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
        out_dir=collection_dir,
        samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    new_unique_seds_fasta, new_unique_groups = main_pipeline.run_mmseq_with_map_unique_seqs(
        faa_file=combined_faa,
        map_file=combined_faa_map,
        out_dir=temp_dir,        
        threads=threads)
    #json.dump(new_unique_groups, open(os.path.join(collection_dir, 'new_unique_groups.json'), 'w'), indent=4, sort_keys=True)
    concat_unique_seq_fasta=concat2fasta(old_unique_seqs,new_unique_seds_fasta,os.path.join(temp_dir,"concat_unique_seqs.fasta"))
    combined_unique_seqs_fasta, combined_unique_group=main_pipeline.run_mmseq_unique_seqs(
        faa_file=concat_unique_seq_fasta,
        out_dir=temp_dir,        
        threads=threads)
    #json.dump(combined_unique_group, open(os.path.join(collection_dir, 'combined_unique_group.json'), 'w'), indent=4, sort_keys=True)
    flat_unique_seq=add_sample_pipeline.flat_combined_unique_seq(new_unique_groups,old_unique_groups,combined_unique_group)
    json.dump(flat_unique_seq, open(os.path.join(collection_dir, 'unique_groups.json'), 'w'), indent=4, sort_keys=True)
    groups_representative_fasta, similar_groups = main_pipeline.run_mmseq_with_map_similar_seqs(
        faa_file=combined_unique_seqs_fasta,
        
        out_dir=temp_dir,      
        threads=threads)
    blast_result = main_pipeline.pairwise_alignment_partion_diamond(
      
        #database_fasta = groups_representative_fasta,
        input_fasta = groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result,
        threads=threads)
    #json.dump(similar_groups, open(os.path.join(collection_dir, 'similar_groups.json'), 'w'), indent=4, sort_keys=True)
    
    inflated_clusters, clusters, unique_seq_by_similar_group = add_sample_pipeline.reinflate_clusters_with_hirachical_group(
        old_unique_groups=old_unique_groups,
        new_unique_groups=new_unique_groups,
        combined_unique_group=combined_unique_group,
        similar_groups=similar_groups,
        mcl_file=mcl_file)
    logger.info(f'len new inflated_clusters = {len(inflated_clusters)} len groups = {len(clusters)}')
    #merge clusters
    
    #json.dump(inflated_clusters, open(os.path.join(collection_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    split_clusters = post_analysis.split_paralogs(
        gene_position_fn=gene_position_fn,
        unsplit_clusters= inflated_clusters,
        dontsplit=args.dont_split
        )
    logger.info(f'len new split_clusters = {len(split_clusters)}')
    json.dump(split_clusters, open(os.path.join(collection_dir, 'split_clusters.json'), 'w'), indent=4, sort_keys=True)

    annotated_clusters = post_analysis.annotate_cluster_progressive_unique_seq(
        unique_seq_by_similar_group=unique_seq_by_similar_group,
        unlabeled_clusters=split_clusters,
        gene_annotation_fn=gene_annotation_fn,
        out_dir = temp_dir
        )
    json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    logger.info(f'samples after extend = {len(old_samples)}')
    output.create_outputs(annotated_clusters,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    #json.dump(annotated_clusters, open(os.path.join(collection_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #print(annotated_clusters)
    if args.alignment:
        samples_dir = os.path.join(collection_dir, 'samples')
        if not os.path.exists(samples_dir):
            raise Exception(f'{samples_dir} does not exist')
        #post_analysis.run_gene_alignment(annotated_clusters, new_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads)
        post_analysis.run_gene_alignment(annotated_clusters, old_samples, collection_dir, args.alignment, coverage_threshold=args.ratio_coverage, threads=threads,poa=args.poa,add=True)
        post_analysis.create_poa_protein_consensus(annotated_clusters,collection_dir,threads=threads)
        post_analysis.make_protein_consensus_db(annotated_clusters,collection_dir,threads=threads)
    # output for next run
    #main_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv.gz')
    #main_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv.gz')

    #Replace the main existing files by the new ones
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
    shutil.copy(combined_unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_hash(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    gene_hash= json.load(open(os.path.join(collection_dir, 'gene_hash.json'), 'r'))
    group_hash= json.load(open(os.path.join(collection_dir, 'group_hash.json'), 'r'))
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    old_filtered_blast_result=os.path.join(collection_dir, 'similar.tsv')
    old_similar_seqs=os.path.join(collection_dir, 'similar.fasta')
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')

    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa = data_preparation.combine_proteins(
        out_dir=collection_dir,
        samples=new_samples)
    #combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
    #    out_dir=collection_dir,
    #    samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    #shutil.rmtree(os.path.join(collection_dir,"samples"))
    gene_hash,remain_new_faa=main_pipeline.add_hash_unique_seqs(
        old_gene_hash=gene_hash,
        faa_file=combined_faa, 
        out_dir=temp_dir,
        )
    new_unique_seqs_fasta, new_gene_hash, new_map_rep_hash= main_pipeline.run_hash_unique_seqs(
        faa_file=remain_new_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    
    #merge gene_hash and map gene hash
    map_rep_hash={}
    for h in gene_hash:
        map_rep_hash[gene_hash[h][0]]=h
    
    for h in new_gene_hash:
        gene_hash[str(h)]=new_gene_hash[h]
    for id in new_map_rep_hash:
        map_rep_hash[id]=new_map_rep_hash[id]
    logging.info(f'size gene_hash after extend: {len(gene_hash)}')
    
    #merge unique fasta
    combined_unique_seqs_fasta=concat2fasta(old_unique_seqs,new_unique_seqs_fasta,os.path.join(temp_dir,"concat_unique_seqs.fasta"))
    
    new_groups_representative_fasta, new_similar_groups = main_pipeline.run_mmseq_with_map_similar_seqs(
        faa_file=new_unique_seqs_fasta,
        
        out_dir=temp_dir,      
        threads=threads)
    similar_groups={}
    for  h in group_hash:
        first_id=gene_hash[h][0]
        similar_groups[first_id]=[]
        for hh in group_hash[h]:
           #print(type(h))
            #print(type(hh))
            similar_groups[first_id].append(gene_hash[str(hh)][0])
        if similar_groups[first_id][0]==first_id:
            similar_groups[first_id].pop(0)
    logging.info(f'old similar groups: {len(similar_groups)}')
    
    for g in new_similar_groups:
        h=new_map_rep_hash[g]
        group_hash[h]=[h]
        for id in new_similar_groups[g]:
            group_hash[h].append(map_rep_hash[id])
        similar_groups[g]=new_similar_groups[g]
    logging.info(f'combined similar groups: {len(similar_groups)}')
    
    new_blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = new_groups_representative_fasta,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    new_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=new_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    concat_filter_blast_result=appendTextfile(old_filtered_blast_result,new_filtered_blast_result)
    pairwise_blast_result = main_pipeline.pairwise_alignment_diamond_split_db(
      
        database_fasta = old_similar_seqs,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)
    pairwise_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=pairwise_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    concat_filter_blast_result=appendTextfile(concat_filter_blast_result,pairwise_filtered_blast_result)
        
    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = concat_filter_blast_result,
        threads=threads)
    
   
    
    #json.dump(similar_groups, open(os.path.join(collection_dir, 'similar_groups.json'), 'w'), indent=4, sort_keys=True)
    
    inflated_clusters, groups = main_pipeline.reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)  
    #json.dump(inflated_clusters, open(os.path.join(collection_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    annotated_clusters = post_analysis.annotate_cluster_hash(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn,
        gene_hash=gene_hash,
       
        map_gene_hash=map_rep_hash)
    
    annotated_clusters_file=os.path.join(collection_dir, 'clusters.json')
    
    json.dump(annotated_clusters, open(annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    #annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    #annotated_clusters_file=post_analysis.expandClusterMembersByHash(annotated_clusters_file,gene_hash,map_rep_hash)
    json.dump(gene_hash, open(os.path.join(collection_dir, 'gene_hash.json'), 'w'), indent=4, sort_keys=True)
    group_hash_file=os.path.join(collection_dir, 'group_hash.json')
    
    json.dump(group_hash, open(group_hash_file, 'w'), indent=4, sort_keys=True)
    
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs_from_hash(annotated_clusters_file,gene_hash,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    combined_similar_group_seqs_fasta=appendTextfile(old_similar_seqs,new_groups_representative_fasta)
 
    main_similar = os.path.join(collection_dir, 'similar.tsv')
   
    shutil.move(concat_filter_blast_result, main_similar)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    
    logger.info(f'samples after extend = {len(old_samples)}')
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
       
    shutil.copy(combined_unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_hash2(args):
    starttime = datetime.now()
    mem_usage = mem_report(0, "start run_add_sample_pipeline_hash2")
    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    gene_hash= json.load(open(os.path.join(collection_dir, 'gene_hash.json'), 'r'))
    group_hash= json.load(open(os.path.join(collection_dir, 'group_hash.json'), 'r'))
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    old_filtered_blast_result=os.path.join(collection_dir, 'similar.tsv')
    old_similar_seqs=os.path.join(collection_dir, 'similar.fasta')
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 1: after load old data")
    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 2: after extract_proteins_tofile")
    
    combined_faa = data_preparation.combine_proteins(
        out_dir=collection_dir,
        samples=new_samples)
    #combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
    #    out_dir=collection_dir,
    #    samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    #shutil.rmtree(os.path.join(collection_dir,"samples"))'
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 3: after combine_proteins")
    
    gene_hash,remain_new_faa=main_pipeline.add_hash_unique_seqs(
        old_gene_hash=gene_hash,
        faa_file=combined_faa, 
        out_dir=temp_dir,
        )
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 4: after add_hash_unique_seqs")
    new_unique_seqs_fasta, new_gene_hash, new_map_rep_hash= main_pipeline.run_hash_unique_seqs(
        faa_file=remain_new_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 5: after run_hash_unique_seqs")
    
    matched_groups,unmatched_u_seqs=main_pipeline.diamond_cd_hit_2d_split(old_similar_seqs,new_unique_seqs_fasta,
        out_dir = os.path.join(temp_dir, 'merge'),
        evalue = args.evalue,
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 6: after diamond_cd_hit_2d_split")
    
    #merge gene_hash and map gene hash
    map_rep_hash={}
    for h in gene_hash:
        map_rep_hash[gene_hash[h][0]]=h
    
    for h in new_gene_hash:
        gene_hash[str(h)]=new_gene_hash[h]
    for id in new_map_rep_hash:
        map_rep_hash[id]=new_map_rep_hash[id]
    logging.info(f'size gene_hash after extend: {len(gene_hash)}')
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 7: after merge rep_hash gene_hash")
    
    #merge unique fasta
    combined_unique_seqs_fasta=concat2fasta(old_unique_seqs,new_unique_seqs_fasta,os.path.join(temp_dir,"concat_unique_seqs.fasta"))
    #merge group_gene
    for g in matched_groups:
        h=map_rep_hash[g]
        for i in matched_groups[g]:
            group_hash[h].append(map_rep_hash[i])
    
    new_groups_representative_fasta, new_similar_groups = main_pipeline.run_mmseq_with_map_similar_seqs(
        faa_file=unmatched_u_seqs,
        
        out_dir=temp_dir,      
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 8: after run_mmseq_with_map_similar_seqs")
    
     
    similar_groups={}
    for  h in group_hash:
        first_id=gene_hash[h][0]
        similar_groups[first_id]=[]
        for hh in group_hash[h]:
           #print(type(h))
            #print(type(hh))
            similar_groups[first_id].append(gene_hash[str(hh)][0])
        if similar_groups[first_id][0]==first_id:
            similar_groups[first_id].pop(0)
    logging.info(f'old similar groups: {len(similar_groups)}')
    
    for g in new_similar_groups:
        h=new_map_rep_hash[g]
        group_hash[h]=[h]
        for id in new_similar_groups[g]:
            group_hash[h].append(map_rep_hash[id])
        similar_groups[g]=new_similar_groups[g]
    logging.info(f'combined similar groups: {len(similar_groups)}')
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 9: after merged similar_groups")
    
    new_blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = new_groups_representative_fasta,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 10: after pairwise_alignment_diamond")
    
    new_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=new_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 11: after filter_blast_result")
    
    concat_filter_blast_result=appendTextfile(old_filtered_blast_result,new_filtered_blast_result)
    pairwise_blast_result = main_pipeline.pairwise_alignment_diamond_split_db(
      
        database_fasta = old_similar_seqs,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 12: after pairwise_alignment_diamond_split_db")
    
    pairwise_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=pairwise_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    
    concat_filter_blast_result=appendTextfile(concat_filter_blast_result,pairwise_filtered_blast_result)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 13: after filter_blast_result")
        
    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = concat_filter_blast_result,
        threads=threads)
    
   
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 14: after cluster_with_mcl")
    
    #json.dump(similar_groups, open(os.path.join(collection_dir, 'similar_groups.json'), 'w'), indent=4, sort_keys=True)
    
    inflated_clusters, groups = main_pipeline.reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)  
    #json.dump(inflated_clusters, open(os.path.join(collection_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 15: after reinflate_clusters")
    
    annotated_clusters = post_analysis.annotate_cluster_hash(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn,
        gene_hash=gene_hash,
       
        map_gene_hash=map_rep_hash)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 16: after annotate_cluster_hash")
    
    annotated_clusters_file=os.path.join(collection_dir, 'clusters.json')
    
    json.dump(annotated_clusters, open(annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    #annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    #annotated_clusters_file=post_analysis.expandClusterMembersByHash(annotated_clusters_file,gene_hash,map_rep_hash)
    json.dump(gene_hash, open(os.path.join(collection_dir, 'gene_hash.json'), 'w'), indent=4, sort_keys=True)
    group_hash_file=os.path.join(collection_dir, 'group_hash.json')
    
    json.dump(group_hash, open(group_hash_file, 'w'), indent=4, sort_keys=True)
    
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs_from_hash(annotated_clusters_file,gene_hash,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 17: after create_outputs_from_hash")
    
    combined_similar_group_seqs_fasta=appendTextfile(old_similar_seqs,new_groups_representative_fasta)
 
    main_similar = os.path.join(collection_dir, 'similar.tsv')
   
    shutil.move(concat_filter_blast_result, main_similar)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    
    logger.info(f'samples after extend = {len(old_samples)}')
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
       
    shutil.copy(combined_unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 18: after all")
    
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_hash_opt(args):
    starttime = datetime.now()
    mem_usage = mem_report(0, "start run_add_sample_pipeline_hash2")
    collection_dir = args.collection_dir
    hash_dir=os.path.join(collection_dir, 'hash')
    

    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    gene_hash= json.load(open(os.path.join(collection_dir, 'gene_hash.json'), 'r'))
    group_hash= json.load(open(os.path.join(collection_dir, 'group_hash.json'), 'r'))
    map_rep_hash=json.load(open(os.path.join(collection_dir, 'map_rep_hash.json'), 'r'))
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    old_filtered_blast_result=os.path.join(collection_dir, 'similar.tsv')
    old_similar_seqs=os.path.join(collection_dir, 'similar.fasta')
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 1: after load old data")
    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 2: after extract_proteins_tofile")
    
    combined_faa = data_preparation.combine_proteins(
        out_dir=collection_dir,
        samples=new_samples)
    #combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
    #    out_dir=collection_dir,
    #    samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    #shutil.rmtree(os.path.join(collection_dir,"samples"))'
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 3: after combine_proteins")
    
    gene_hash,remain_new_faa=main_pipeline.add_hash_unique_seqs_opt(
        old_gene_hash=gene_hash,
        faa_file=combined_faa, 
        hash_dir=hash_dir,
        out_dir=temp_dir,
        )
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 4: after add_hash_unique_seqs")
    new_unique_seqs_fasta, new_gene_hash, new_map_rep_hash= main_pipeline.run_hash_unique_seqs_opt(
        faa_file=remain_new_faa,
      
        out_dir=temp_dir,     
        hash_dir=hash_dir, 
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 5: after run_hash_unique_seqs")
    
    matched_groups,unmatched_u_seqs=main_pipeline.diamond_cd_hit_2d_split(old_similar_seqs,new_unique_seqs_fasta,
        out_dir = os.path.join(temp_dir, 'merge'),
        evalue = args.evalue,
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 6: after diamond_cd_hit_2d_split")
    
    #merge gene_hash and map gene hash
    
    
    for h in new_gene_hash:
        gene_hash[str(h)]=[]
    for id in new_map_rep_hash:
        map_rep_hash[id]=new_map_rep_hash[id]
    map_hash_rep={}
    for id in map_rep_hash:
        map_hash_rep[str(map_rep_hash[id])]=id
    logging.info(f'size gene_hash after extend: {len(gene_hash)}')
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 7: after merge rep_hash gene_hash")
    
    #merge unique fasta
    combined_unique_seqs_fasta=concat2fasta(old_unique_seqs,new_unique_seqs_fasta,os.path.join(temp_dir,"concat_unique_seqs.fasta"))
    #merge group_gene
    for g in matched_groups:
        h=map_rep_hash[g]
        for i in matched_groups[g]:
            group_hash[h].append(map_rep_hash[i])
    
    new_groups_representative_fasta, new_similar_groups = main_pipeline.run_mmseq_with_map_similar_seqs(
        faa_file=unmatched_u_seqs,
        
        out_dir=temp_dir,      
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 8: after run_mmseq_with_map_similar_seqs")
    
     
    similar_groups={}
    for  h in group_hash:
        first_id=map_hash_rep[h]
        similar_groups[first_id]=[]
        for hh in group_hash[h]:
           #print(type(h))
            #print(type(hh))
            similar_groups[first_id].append(map_hash_rep[str(hh)])
        if similar_groups[first_id][0]==first_id:
            similar_groups[first_id].pop(0)
    logging.info(f'old similar groups: {len(similar_groups)}')
    
    for g in new_similar_groups:
        h=new_map_rep_hash[g]
        group_hash[h]=[h]
        for id in new_similar_groups[g]:
            group_hash[h].append(map_rep_hash[id])
        similar_groups[g]=new_similar_groups[g]
    logging.info(f'combined similar groups: {len(similar_groups)}')
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 9: after merged similar_groups")
    
    new_blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = new_groups_representative_fasta,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 10: after pairwise_alignment_diamond")
    
    new_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=new_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 11: after filter_blast_result")
    
    concat_filter_blast_result=appendTextfile(old_filtered_blast_result,new_filtered_blast_result)
    pairwise_blast_result = main_pipeline.pairwise_alignment_diamond_split_db(
      
        database_fasta = old_similar_seqs,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 12: after pairwise_alignment_diamond_split_db")
    
    pairwise_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=pairwise_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    
    concat_filter_blast_result=appendTextfile(concat_filter_blast_result,pairwise_filtered_blast_result)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 13: after filter_blast_result")
        
    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = concat_filter_blast_result,
        threads=threads)
    
   
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 14: after cluster_with_mcl")
    
    #json.dump(similar_groups, open(os.path.join(collection_dir, 'similar_groups.json'), 'w'), indent=4, sort_keys=True)
    
    inflated_clusters, groups = main_pipeline.reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)  
    #json.dump(inflated_clusters, open(os.path.join(collection_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 15: after reinflate_clusters")
    
    annotated_clusters = post_analysis.annotate_cluster_hash_opt(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn,
       
       
        map_gene_hash=map_rep_hash)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 16: after annotate_cluster_hash")
    
    annotated_clusters_file=os.path.join(collection_dir, 'clusters.json')
    
    json.dump(annotated_clusters, open(annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    #annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    #annotated_clusters_file=post_analysis.expandClusterMembersByHash(annotated_clusters_file,gene_hash,map_rep_hash)
    json.dump(gene_hash, open(os.path.join(collection_dir, 'gene_hash.json'), 'w'), indent=4, sort_keys=True)
    group_hash_file=os.path.join(collection_dir, 'group_hash.json')
    
    json.dump(group_hash, open(group_hash_file, 'w'), indent=4, sort_keys=True)
    map_rep_hash_file=os.path.join(collection_dir, 'map_rep_hash.json')
    json.dump(map_rep_hash, open(map_rep_hash_file, 'w'), indent=4, sort_keys=True)
    
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    output.create_outputs_from_hash_opt(annotated_clusters_file,gene_hash,old_samples,collection_dir,hash_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 17: after create_outputs_from_hash")
    
    combined_similar_group_seqs_fasta=appendTextfile(old_similar_seqs,new_groups_representative_fasta)
 
    main_similar = os.path.join(collection_dir, 'similar.tsv')
   
    shutil.move(concat_filter_blast_result, main_similar)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    
    logger.info(f'samples after extend = {len(old_samples)}')
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
       
    shutil.copy(combined_unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    #shutil.move(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)
    mem_usage = mem_report(mem_usage, "run_add_sample_pipeline_hash2-checkpoint 18: after all")
    
    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def run_add_sample_pipeline_hash_opt2(args):
    starttime = datetime.now()

    collection_dir = args.collection_dir
    if not os.path.exists(collection_dir):
        raise Exception(f'{collection_dir} does not exist')
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    diamond=(args.blast=='diamond')

    identity = args.identity
    evalue = args.evalue


    temp_dir = os.path.join(collection_dir, 'temp')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
        pass
    else:
        os.makedirs(temp_dir)

    gene_annotation_fn = os.path.join(temp_dir, 'gene_annotation.csv')
    gene_position_fn = os.path.join(temp_dir, 'gene_position.csv')
    combined_map= os.path.join(collection_dir, 'combined.map')
    gene_hash= json.load(open(os.path.join(collection_dir, 'gene_hash.json'), 'r'))
    group_hash= json.load(open(os.path.join(collection_dir, 'group_hash.json'), 'r'))
    #os.makedirs(group_seq_dir)
    # Check required files
    existing_gene_annotation_fn = os.path.join(collection_dir, 'gene_annotation.csv')
    if not os.path.isfile(existing_gene_annotation_fn):
        raise Exception(f'{existing_gene_annotation_fn} does not exist')
    #gene_annotation = output.import_gene_annotation(gene_annotation_file)

    existing_gene_position_fn = os.path.join(collection_dir, 'gene_position.csv')
    #gene_position = json.load(open(os.path.join(collection_dir, 'gene_position.json'), 'r'))

    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    #old_clusters = json.load(open(os.path.join(collection_dir, 'clusters.json'), 'r'))
    old_unique_seqs=os.path.join(collection_dir, 'unique_seqs.fasta')
    old_filtered_blast_result=os.path.join(collection_dir, 'similar.tsv')
    old_similar_seqs=os.path.join(collection_dir, 'similar.fasta')
    # collect new samples
    sample_id_list = [sample['id'] for sample in old_samples]
    new_samples = collect_sample(sample_id_list, args)
    if len(new_samples) == 0:
        raise Exception(f'There must be at least one new sample')

    # data preparation
    data_preparation.extract_proteins_tofile(
        samples=new_samples,
        out_dir=collection_dir,
        gene_annotation_fn = gene_annotation_fn,
        gene_position_fn = gene_position_fn,
        table=args.table,
        existing_gene_annotation_fn=existing_gene_annotation_fn,
        existing_gene_position_fn=existing_gene_position_fn,
        threads=threads,
        )
    combined_faa,combined_faa_map = data_preparation.add_combine_proteins_with_maps(
        out_dir=collection_dir,
        samples=new_samples,
        combined_map=combined_map
        )
    #combined_faa, combined_faa_map = data_preparation.combine_proteins_with_maps(
    #    out_dir=collection_dir,
    #    samples=new_samples)
    #combined_faa_file, combined_faa_map=data_preparation.make_combine_maps(new_combined_faa,collection_dir)
    #shutil.rmtree(os.path.join(collection_dir,"samples"))
    gene_hash,remain_new_faa=main_pipeline.add_hash_unique_seqs(
        old_gene_hash=gene_hash,
        faa_file=combined_faa, 
        out_dir=temp_dir,
        )
    new_unique_seqs_fasta, new_gene_hash, new_map_rep_hash= main_pipeline.run_hash_unique_seqs(
        faa_file=remain_new_faa,
      
        out_dir=temp_dir,      
        threads=threads)
    
    #merge gene_hash and map gene hash
    map_rep_hash={}
    for h in gene_hash:
        map_rep_hash[gene_hash[h][0]]=h
    
    for h in new_gene_hash:
        gene_hash[str(h)]=new_gene_hash[h]
    for id in new_map_rep_hash:
        map_rep_hash[id]=new_map_rep_hash[id]
    logging.info(f'size gene_hash after extend: {len(gene_hash)}')
    
    #merge unique fasta
    combined_unique_seqs_fasta=concat2fasta(old_unique_seqs,new_unique_seqs_fasta,os.path.join(temp_dir,"concat_unique_seqs.fasta"))
    
    new_groups_representative_fasta, new_similar_groups = main_pipeline.run_mmseq_with_map_similar_seqs(
        faa_file=new_unique_seqs_fasta,
        
        out_dir=temp_dir,      
        threads=threads)
    similar_groups={}
    for  h in group_hash:
        first_id=gene_hash[h][0]
        similar_groups[first_id]=[]
        for hh in group_hash[h]:
           #print(type(h))
            #print(type(hh))
            similar_groups[first_id].append(gene_hash[str(hh)][0])
        if similar_groups[first_id][0]==first_id:
            similar_groups[first_id].pop(0)
    logging.info(f'old similar groups: {len(similar_groups)}')
    
    for g in new_similar_groups:
        h=new_map_rep_hash[g]
        group_hash[h]=[h]
        for id in new_similar_groups[g]:
            group_hash[h].append(map_rep_hash[id])
        similar_groups[g]=new_similar_groups[g]
    logging.info(f'combined similar groups: {len(similar_groups)}')
    
    new_blast_result = main_pipeline.pairwise_alignment_diamond(
      
        database_fasta = new_groups_representative_fasta,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    new_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=new_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    concat_filter_blast_result=appendTextfile(old_filtered_blast_result,new_filtered_blast_result)
    pairwise_blast_result = main_pipeline.pairwise_alignment_diamond_split_db(
      
        database_fasta = old_similar_seqs,
        query_fasta = new_groups_representative_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)
    pairwise_filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=pairwise_blast_result,
        out_dir = temp_dir,
        identity=args.identity,
        length_difference=args.LD,
        alignment_coverage_short=args.AS,
        alignment_coverage_long=args.AL)
    concat_filter_blast_result=appendTextfile(concat_filter_blast_result,pairwise_filtered_blast_result)
        
    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = concat_filter_blast_result,
        threads=threads)
    
   
    
    #json.dump(similar_groups, open(os.path.join(collection_dir, 'similar_groups.json'), 'w'), indent=4, sort_keys=True)
    
    inflated_clusters, groups = main_pipeline.reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)  
    #json.dump(inflated_clusters, open(os.path.join(collection_dir, 'inflated_clusters.json'), 'w'), indent=4, sort_keys=True)
    
    no_annotated_clusters = post_analysis.make_no_annotated_cluster(
        unlabeled_clusters=inflated_clusters,
        gene_annotation_fn=gene_annotation_fn,
        gene_hash=gene_hash,
       
        map_gene_hash=map_rep_hash)
    
    no_annotated_clusters_file=os.path.join(collection_dir, 'clusters.json')
    
    json.dump(no_annotated_clusters, open(no_annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    #annotated_clusters_file=save_clusters(annotated_clusters,out_dir)
    #annotated_clusters_file=post_analysis.expandClusterMembersByHash(annotated_clusters_file,gene_hash,map_rep_hash)
    json.dump(gene_hash, open(os.path.join(collection_dir, 'gene_hash.json'), 'w'), indent=4, sort_keys=True)
    group_hash_file=os.path.join(collection_dir, 'group_hash.json')
    
    json.dump(group_hash, open(group_hash_file, 'w'), indent=4, sort_keys=True)
    
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    old_samples.extend(new_samples)
    #annotated_clusters=post_analysis.merge_new_cluster_to_old_clusters(annotated_clusters,clusters_by_ref)
    #json.dump(annotated_clusters, open(os.path.join(out_dir, 'annotated_clusters.json'), 'w'), indent=4, sort_keys=True)
    #output.create_outputs_from_hash(annotated_clusters_file,gene_hash,old_samples,collection_dir,t_core=args.core,t_soft=args.soft,t_shell=args.shell)
    combined_similar_group_seqs_fasta=appendTextfile(old_similar_seqs,new_groups_representative_fasta)
 
    main_similar = os.path.join(collection_dir, 'similar.tsv')
   
    shutil.move(concat_filter_blast_result, main_similar)
    shutil.rmtree(os.path.join(collection_dir, 'samples'))
    
    logger.info(f'samples after extend = {len(old_samples)}')
    shutil.copy(gene_annotation_fn, existing_gene_annotation_fn)
    shutil.copy(gene_position_fn, existing_gene_position_fn)
       
    shutil.copy(combined_unique_seqs_fasta, old_unique_seqs)
    #output.export_gene_annotation(gene_annotation, collection_dir)
    #json.dump(gene_position, open(os.path.join(collection_dir, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    
    json.dump(old_samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    #add_sample_pipeline.combine_representative(not_match_represent_faa, old_represent_faa, collection_dir)
    #json.dump(new_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    shutil.move(combined_faa_map, combined_map)
    #shutil.copy(combined_blast_result, os.path.join(collection_dir, 'blast.tsv'))
    #cmd = f'gzip -c {combined_blast_result} > ' + os.path.join(collection_dir, 'blast.tsv.gz')
    #cmd = f'mv {combined_blast_result}  ' + os.path.join(collection_dir, 'blast.tsv')
    #os.system(cmd)

    elapsed = datetime.now() - starttime
    logging.info(f'Done -- time cli {str(elapsed)}')
    logging.info(f'Done -- time taken {str(elapsed)}')
def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    main_cmd = subparsers.add_parser(
        'main',
        description='Main pipeline: run pan-genome analysis for the first time',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    main_cmd.set_defaults(func=run_main)
    main_cmd.add_argument('-m', '--mode', help='mode running', required=True, type=str)
    main_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    main_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    main_cmd.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    main_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    main_cmd.add_argument('-b', '--blast', help='method for all-against-all alignment', default='diamond', action='store', choices=['diamond', 'blast'])
    main_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.70, type=float)
    main_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0.70, type=float)
    main_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    main_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0, type=float)
    main_cmd.add_argument('-c','--cov' ,help='coverage of grouping step', default=90, type=float)
    main_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    main_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    main_cmd.add_argument('--table', help='codon table', default=11, type=int)
    main_cmd.add_argument('-a', '--alignment', help='run alignment for each gene cluster', default=None, choices=['nucleotide', 'protein'])
    main_cmd.add_argument('-r', '--ratio-coverage', help='Ratio of coverage to align', default=0.0, type=float)
    main_cmd.add_argument('--poa', help='Alignment with POA', default=False, action='store_true')
    main_cmd.add_argument('--core', help='Percentage of core genes', default=0.99, type=float)
    main_cmd.add_argument('--soft', help='Percentage of soft core genes', default=0.95, type=float)
    main_cmd.add_argument('--shell', help='Percentage of shell genes', default=0.15, type=float)
    main_cmd.add_argument('--ksize', help='K-mer length', default=3, type=int)
    main_cmd.add_argument('--similarity', help='similarity threadshold', default=0.70, type=float)
    #main_cmd.add_argument('--diff_len', help='diff len threadshold', default=0.30, type=float)

    add_cmd = subparsers.add_parser(
        'add',
        description='Add pipeline: add sample into previous collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=run_add)
    add_cmd.add_argument('-m', '--mode', help='mode running', required=True, type=str)
    add_cmd.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    add_cmd.add_argument('-f', '--tsv', help='tsv input file',default=None, type=str)
    add_cmd.add_argument('-c', '--collection-dir', help='previous collection directory', required=True, type=str)
    add_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    add_cmd.add_argument('-b', '--blast', help='method for all-against-all alignment', default='diamond', action='store', choices=['diamond', 'blast'])
    add_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=0.70, type=float)
    add_cmd.add_argument('--LD', help='length difference cutoff between two sequences', default=0.70, type=float)
    add_cmd.add_argument('--AL', help='alignment coverage for the longer sequence', default=0, type=float)
    add_cmd.add_argument('--AS', help='alignment coverage for the shorter sequence', default=0, type=float)
    add_cmd.add_argument('-e', '--evalue', help='Blast evalue', default=1E-6, type=float)
    add_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--table', help='codon table', default=11, type=int)
    add_cmd.add_argument('-a', '--alignment', help='run alignment for each gene cluster', default=None, choices=['nucleotide', 'protein'])
    add_cmd.add_argument('-r', '--ratio-coverage', help='Ratio of coverage to align', default=0.0, type=float)
    add_cmd.add_argument('--poa', help='Alignment with POA', default=False, action='store_true')
   
    add_cmd.add_argument('--core', help='Percentage of core genes', default=0.99, type=float)
    add_cmd.add_argument('--soft', help='Percentage of soft core genes', default=0.95, type=float)
    add_cmd.add_argument('--shell', help='Percentage of shell genes', default=0.15, type=float)
    
    ref_cmd = subparsers.add_parser(
        'build',
        description='Build reference gene famyly clusters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ref_cmd.set_defaults(func=build_reference_gene_family_db_hmmer)
    ref_cmd.add_argument('-i', '--input', help='reference gene sequences in fasta format',required=True, type=str)
    ref_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=1, type=int)
    ref_cmd.add_argument('-d', '--identity', help='identity', default=0.7, type=float)
    ref_cmd.add_argument('-c', '--coverage', help='coverage', default=0.7, type=float)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
