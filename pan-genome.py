import argparse
import os
import shutil
import logging
import json
from glob import glob
from pan_genome import *

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)

def run_main_pipeline(args):
    collection_dir = args.collection_dir
    threads = args.threads
    dontsplit = args.dont_split
    fasta = args.fasta
    diamond = args.diamond
    identity = args.identity
    number = args.number
    
    samples = []
    fasta_ext = ('.fasta', '.fna', 'ffn')
    for path in args.gff_files:
        dir_name = os.path.dirname(path)
        base_name = os.path.basename(path)
        sample_id = base_name.rsplit('.', 1)[0]
        sample = {'id':sample_id, 'gff_file':path}
        if fasta == True:
            try:
                fasta_file = [f for f in glob(dir_name+'/'+sample_id+'*') if f.endswith(fasta_ext)][0]
                sample['fasta_file'] = fasta_file
            except:
                raise Exception(f'The corresponding fasta file of {sample_id} does not exist')
        samples.append(sample)
    samples.sort(key= lambda x:x['id'])

    temp_dir = os.path.join(collection_dir, 'temp')
    timing_log = os.path.join(collection_dir, 'time.log')
    if not os.path.exists(collection_dir):
        os.makedirs(collection_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # data preparation
    gene_annotation = {}
    data_preparation.extract_proteins(
        samples=samples,
        out_dir=collection_dir,
        gene_annotation = gene_annotation,
        fasta=fasta,
        threads=threads
        )

    subset = samples[0:number]
    remain = samples[number:]
    
    # run a subset of collection
    subset_dir = os.path.join(temp_dir, 'subset')
    if not os.path.exists(subset_dir):
        os.makedirs(subset_dir)

    subset_combined_faa = data_preparation.combine_proteins(
        out_dir=subset_dir,
        samples=subset,
        timing_log=timing_log)

    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=subset_combined_faa,
        out_dir=subset_dir,
        threads=threads, 
        timing_log=timing_log)

    subset_blast_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(subset_dir, 'blast'),
        identity=identity,
        threads=threads, 
        timing_log=timing_log
        )

    subset_mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = subset_dir,
        blast_result = subset_blast_result,
        threads=threads, 
        timing_log=timing_log)



    # run the remain of collection
    if len(remain) == 0:
        inflated_clusters = main_pipeline.reinflate_clusters(
            cd_hit_clusters=cd_hit_clusters,
            mcl_file=subset_mcl_file
        )
    else:
        remain_dir = os.path.join(temp_dir, 'remain')
        if not os.path.exists(remain_dir):
            os.makedirs(remain_dir)

        remain_combined_faa = data_preparation.combine_proteins(
            out_dir=remain_dir,
            samples=remain,
            timing_log=timing_log)

        not_match_fasta, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
            database_1 = cd_hit_represent_fasta,
            database_2 = remain_combined_faa,
            out_dir = remain_dir,
            threads=threads, 
            timing_log=timing_log)

        not_match_represent_fasta, not_match_clusters = main_pipeline.run_cd_hit(
            faa_file=not_match_fasta,
            out_dir=remain_dir,
            threads=threads, 
            timing_log=timing_log)

        blast_1_result = main_pipeline.pairwise_alignment(
            diamond=diamond,
            database_fasta = cd_hit_represent_fasta,
            query_fasta = not_match_represent_fasta,
            out_dir = os.path.join(remain_dir, 'blast1'),
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )

        blast_2_result = main_pipeline.pairwise_alignment(
            diamond=diamond,
            database_fasta = not_match_represent_fasta,
            query_fasta = not_match_represent_fasta,
            out_dir = os.path.join(remain_dir, 'blast2'),
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )

        combined_blast_result = add_sample_pipeline.combine_blast_results(
            blast_1=subset_blast_result, 
            blast_2=blast_1_result, 
            blast_3=blast_2_result, 
            outdir=remain_dir)

        remain_mcl_file = main_pipeline.cluster_with_mcl(
            out_dir = remain_dir,
            blast_result = combined_blast_result,
            threads=threads, 
            timing_log=timing_log)

        inflated_clusters = add_sample_pipeline.reinflate_clusters(
            cd_hit_clusters=cd_hit_clusters, 
            cd_hit_2d_clusters=cd_hit_2d_clusters,
            not_match_clusters=not_match_clusters,
            mcl_file=remain_mcl_file
            )

    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_annotation=gene_annotation,
        unsplit_clusters= inflated_clusters,
        dontsplit=dontsplit
        )
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters, 
        gene_annotation=gene_annotation)

    # output
    spreadsheet_file = output.create_spreadsheet(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=collection_dir
    )
    rtab_file = output.create_rtab(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=collection_dir
    )
    summary_file = output.create_summary(
        rtab_file=rtab_file, 
        out_dir=collection_dir
    )

    # output for next run
    protein_database = data_preparation.combine_proteins(
        out_dir=collection_dir,
        samples=samples,
        timing_log=timing_log)

    json.dump(gene_annotation, open(os.path.join(collection_dir, 'gene_annotation.json'), 'w'), indent=4, sort_keys=True)
    json.dump(inflated_clusters, open(os.path.join(collection_dir, 'unsplit_clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)
    

def run_add_sample_pipeline(args):
    collection_dir = args.collection_dir
    threads = args.threads
    dontsplit = args.dont_split
    fasta = args.fasta
    diamond = args.diamond
    identity = args.identity
    
    samples = []
    fasta_ext = ('.fasta', '.fna', 'ffn')
    for path in args.gff_files:
        dir_name = os.path.dirname(path)
        base_name = os.path.basename(path)
        sample_id = base_name.split('.')[0]
        sample = {'id':sample_id, 'gff_file':path}
        if fasta == True:
            try:
                fasta_file = [f for f in glob(dir_name+'/'+sample_id+'.*') if f.endswith(fasta_ext)][0]
                sample['fasta_file'] = fasta_file
            except:
                raise Exception(f'The corresponding fasta file of {sample_id} does not exist')
        samples.append(sample)
    samples.sort(key= lambda x:x['id'])
    
    temp_dir = os.path.join(collection_dir, 'temp')
    timing_log = os.path.join(collection_dir, 'time.log')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
    else:
        os.makedirs(temp_dir)
    
    # Check needed files
    representative_fasta = os.path.join(collection_dir, 'representative.fasta')
    if not os.path.isfile(representative_fasta):
        raise Exception(f'{representative_fasta} is not exist')
    
    protein_database = os.path.join(collection_dir, 'combined.faa')
    if not os.path.isfile(protein_database):
        raise Exception(f'{protein_database} is not exist')
    
    gene_annotation = json.load(open(os.path.join(collection_dir, 'gene_annotation.json'), 'r'))
    old_clusters = json.load(open(os.path.join(collection_dir, 'unsplit_clusters.json'), 'r'))

    # data preparation
    data_preparation.extract_proteins(
        samples=samples,
        out_dir=collection_dir,
        gene_annotation = gene_annotation,
        fasta=fasta,
        threads=threads
        )
    combined_faa_file = data_preparation.combine_proteins(
        out_dir=temp_dir,
        samples=samples,
        timing_log=timing_log)

    # add sample pipeline
    not_match_fasta, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1 = representative_fasta,
        database_2 = combined_faa_file,
        out_dir = temp_dir,
        identity=identity,
        threads=threads, 
        timing_log=timing_log)
    
    blast_1_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = representative_fasta,
        query_fasta = not_match_fasta,
        out_dir = os.path.join(temp_dir, 'blast1'),
        identity=identity,
        threads=threads, 
        timing_log=timing_log
        )
    
    blast_remain_fasta = add_sample_pipeline.filter_fasta(
        blast_result = blast_1_result, 
        fasta_file = not_match_fasta, 
        out_dir = temp_dir
        )

    blast_2_result = main_pipeline.pairwise_alignment(
        diamond=diamond,
        database_fasta = blast_remain_fasta,
        query_fasta = blast_remain_fasta,
        out_dir = os.path.join(temp_dir, 'blast2'),
        identity=identity,
        threads=threads, 
        timing_log=timing_log
        )

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = blast_2_result,
        threads=threads, 
        timing_log=timing_log)

    inflated_clusters = add_sample_pipeline.reinflate_clusters(
        old_clusters=old_clusters, 
        cd_hit_2d_clusters=cd_hit_2d_clusters, 
        blast_1_result_file=blast_1_result, 
        mcl_clusters=mcl_file
        )

    # post analysis
    split_clusters = post_analysis.split_paralogs(
        gene_annotation=gene_annotation,
        unsplit_clusters= inflated_clusters,
        dontsplit=dontsplit
        )
    annotated_clusters = post_analysis.annotate_cluster(
        unlabeled_clusters=split_clusters, 
        gene_annotation=gene_annotation)
    
    # output
    old_samples = json.load(open(os.path.join(collection_dir, 'samples.json'), 'r'))
    samples.extend(old_samples)
    samples.sort(key= lambda x:x['id'])
    
    spreadsheet_file = output.create_spreadsheet(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=collection_dir
    )
    rtab_file = output.create_rtab(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=collection_dir
    )
    summary_file = output.create_summary(
        rtab_file=rtab_file, 
        out_dir=collection_dir
    )

    # output for next run
    os.system(f'cat {combined_faa_file} >> {protein_database}')

    representative_fasta = output.create_representative_fasta(
        clusters = inflated_clusters,
        gene_annotation = gene_annotation,
        faa_fasta = protein_database,
        out_dir = collection_dir
        )

    json.dump(gene_annotation, open(os.path.join(collection_dir, 'gene_annotation.json'), 'w'), indent=4, sort_keys=True)
    json.dump(inflated_clusters, open(os.path.join(collection_dir, 'unsplit_clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    
    main_cmd = subparsers.add_parser(
        'main',
        description='Main pipeline: run pan-genome analysis for the first time',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    main_cmd.set_defaults(func=run_main_pipeline)
    main_cmd.add_argument('gff_files', help='a.gff b.gff ... (*.gff)', type=str, nargs='+')
    main_cmd.add_argument('-c', '--collection-dir', help='collection directory', required=True, type=str)
    main_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    main_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    main_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=95, type=float)
    main_cmd.add_argument('-f', '--fasta', help='fasta files are seperated from gff files. (fasta file must have the same name, be in the same folder of coresponding gff file, and have one of following extension: .fasta .fna .fnn)', default=False, action='store_true')
    main_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)
    main_cmd.add_argument('-n', '--number', help='number of samples which are analysed first', default=50, type=int)

    add_cmd = subparsers.add_parser(
        'add',
        description='Add pipeline: add sample into previous collection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=run_add_sample_pipeline)
    add_cmd.add_argument('gff_files', help='a.gff b.gff ... (*.gff)', type=str, nargs='+')
    add_cmd.add_argument('-c', '--collection-dir', help='previous collection directory', required=True, type=str)
    add_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    add_cmd.add_argument('-d', '--diamond', help='use Diamond for all-agaist-all alignment instead of Blastp', default=False, action='store_true')
    add_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=95, type=float)
    add_cmd.add_argument('-f', '--fasta', help='fasta files are seperated from gff files. (fasta file must have the same name, be in the same folder of coresponding gff file, and have one of following extension: .fasta .fna .fnn)', default=False, action='store_true')
    add_cmd.add_argument('-t', '--threads', help='number of threads to use, 0 for all', default=0, type=int)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()