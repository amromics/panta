import argparse
import os
import shutil
import logging
import json
from glob import glob
from pan_genome import *

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s : %(message)s',
                    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)

def run_main_pipeline(args):
    pan_genome_folder = args.out_dir
    timing_log = args.time_log
    threads = args.threads
    overwrite = args.over_write
    dontsplit = args.dont_split
    fasta = args.fasta
    diamond = args.diamond
    identity = args.identity
    
    samples = []
    fasta_ext = ('.fasta', '.fna', 'ffn')
    for path in args.inputs:
        dir_name = os.path.dirname(path)
        base_name = os.path.basename(path)
        sample_id = base_name.split('.')[0]
        sample = {'id':sample_id, 'gff_file':path}
        if fasta == True:
            try:
                fasta_file = [f for f in glob(dir_name+'/'+sample_id+'*') if f.endswith(fasta_ext)][0]
                sample['fasta_file'] = fasta_file
            except:
                raise Exception(f'The corresponding fasta file of {sample_id} does not exist')
        samples.append(sample)
    
    temp_dir = os.path.join(pan_genome_folder, 'temp')

    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    # Check if pan-genome has run
    pan_genome_output = os.path.join(pan_genome_folder,'summary_statistics.txt')
    if os.path.isfile(pan_genome_output) and (not overwrite):
        return

    # data preparation
    gene_annotation = {}
    gene_position = {}
    data_preparation.extract_proteins(
        samples=samples,
        out_dir=pan_genome_folder,
        gene_annotation = gene_annotation,
        gene_position = gene_position, 
        timing_log=timing_log,
        fasta=fasta
        )
    combined_faa_file = data_preparation.combine_proteins(
        out_dir=temp_dir,
        samples=samples,
        timing_log=timing_log)

    # create protein database
    protein_database = os.path.join(pan_genome_folder, 'protein_database')
    shutil.copyfile(combined_faa_file, protein_database)

    # main pipeline
    cd_hit_represent_fasta, excluded_cluster, cd_hit_clusters = main_pipeline.run_cd_hit_iterative(
        faa_file=combined_faa_file,
        samples=samples,
        out_dir=temp_dir, 
        threads=threads, 
        timing_log=timing_log)

    if diamond == False:
        blast_result = main_pipeline.all_against_all_blast(
            out_dir = os.path.join(temp_dir, 'blast'),
            database_fasta = cd_hit_represent_fasta,
            query_fasta = cd_hit_represent_fasta,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )
    else:
        blast_result = main_pipeline.run_diamond(
            out_dir = os.path.join(temp_dir, 'blast'),
            database_fasta = cd_hit_represent_fasta,
            query_fasta = cd_hit_represent_fasta,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_results = blast_result,
        threads=threads, 
        timing_log=timing_log)

    inflated_clusters = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_clusters,
        mcl_file=mcl_file,
        excluded_cluster=excluded_cluster
    )
    # post analysis
    if dontsplit == False:
        split_clusters = post_analysis.split_paralogs(
            gene_annotation=gene_annotation,
            gene_position=gene_position,
            unsplit_clusters= inflated_clusters
            )
        annotated_clusters = post_analysis.annotate_cluster(
            unlabeled_clusters=split_clusters, 
            gene_annotation=gene_annotation)
    else: 
        annotated_clusters = post_analysis.annotate_cluster(
            unlabeled_clusters=inflated_clusters, 
            gene_annotation=gene_annotation)

    # output
    spreadsheet_file = output.create_spreadsheet(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=pan_genome_folder
    )
    rtab_file = output.create_rtab(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=pan_genome_folder
    )
    summary_file = output.create_summary(
        rtab_file=rtab_file, 
        out_dir=pan_genome_folder
    )
    representative_fasta = output.create_representative_fasta(
        clusters = inflated_clusters,
        gene_annotation = gene_annotation,
        faa_fasta = protein_database,
        out_dir = pan_genome_folder
        )

    json.dump(gene_annotation, open(os.path.join(pan_genome_folder, 'gene_annotation.json'), 'w'), indent=4, sort_keys=True)
    json.dump(gene_position, open(os.path.join(pan_genome_folder, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(inflated_clusters, open(os.path.join(pan_genome_folder, 'unsplit_clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(samples, open(os.path.join(pan_genome_folder, 'samples.json'), 'w'), indent=4, sort_keys=True)
    

def run_add_sample_pipeline(args):
    pan_genome_folder = args.collection_dir
    timing_log = args.time_log
    threads = args.threads
    dontsplit = args.dont_split
    fasta = args.fasta
    diamond = args.diamond
    identity = args.identity
    
    samples = []
    fasta_ext = ('.fasta', '.fna', 'ffn')
    for path in args.inputs:
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
    
    temp_dir = os.path.join(pan_genome_folder, 'temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    # Check if last collection exist
    representative_fasta = os.path.join(pan_genome_folder, 'representative.fasta')
    if not os.path.isfile(representative_fasta):
        raise Exception(f'{representative_fasta} is not exist')
    protein_database = os.path.join(pan_genome_folder, 'protein_database')
    if not os.path.isfile(protein_database):
        raise Exception(f'{protein_database} is not exist')
    gene_annotation = json.load(open(os.path.join(pan_genome_folder, 'gene_annotation.json'), 'r'))
    gene_position = json.load(open(os.path.join(pan_genome_folder, 'gene_position.json'), 'r'))
    old_clusters = json.load(open(os.path.join(pan_genome_folder, 'unsplit_clusters.json'), 'r'))

    # data preparation
    data_preparation.extract_proteins(
        samples=samples,
        out_dir=pan_genome_folder,
        gene_annotation = gene_annotation,
        gene_position = gene_position, 
        timing_log=timing_log,
        fasta=fasta
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
    
    if diamond == False:
        blast_1_result = main_pipeline.all_against_all_blast(
            out_dir = os.path.join(temp_dir, 'blast1'),
            database_fasta = representative_fasta,
            query_fasta = not_match_fasta,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )
    else:
        blast_1_result = main_pipeline.run_diamond(
            out_dir = os.path.join(temp_dir, 'blast1'),
            database_fasta = representative_fasta,
            query_fasta = not_match_fasta,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )
    
    blast_remain_fasta = add_sample_pipeline.filter_fasta(
        blast_result = blast_1_result, 
        fasta_file = not_match_fasta, 
        out_dir = temp_dir
        )
    
    if diamond == False:
        blast_2_result = main_pipeline.all_against_all_blast(
            out_dir = os.path.join(temp_dir, 'blast2'),
            database_fasta = blast_remain_fasta,
            query_fasta = blast_remain_fasta,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )
    else:
        blast_2_result = main_pipeline.run_diamond(
            out_dir = os.path.join(temp_dir, 'blast2'),
            database_fasta = blast_remain_fasta,
            query_fasta = blast_remain_fasta,
            identity=identity,
            threads=threads, 
            timing_log=timing_log
            )

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_results = blast_2_result,
        threads=threads, 
        timing_log=timing_log)

    inflated_clusters = add_sample_pipeline.reinflate_clusters(
        old_clusters=old_clusters, 
        cd_hit_2d_clusters=cd_hit_2d_clusters, 
        blast_1_result_file=blast_1_result, 
        mcl_clusters=mcl_file
        )

    # post analysis
    if dontsplit == False:
        split_clusters = post_analysis.split_paralogs(
            gene_annotation=gene_annotation,
            gene_position=gene_position,
            unsplit_clusters= inflated_clusters
            )
        annotated_clusters = post_analysis.annotate_cluster(
            unlabeled_clusters=split_clusters, 
            gene_annotation=gene_annotation)
    else: 
        annotated_clusters = post_analysis.annotate_cluster(
            unlabeled_clusters=inflated_clusters, 
            gene_annotation=gene_annotation)

    # output
    old_samples = json.load(open(os.path.join(pan_genome_folder, 'samples.json'), 'r'))
    samples.extend(old_samples)

    spreadsheet_file = output.create_spreadsheet(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=pan_genome_folder
    )
    rtab_file = output.create_rtab(
        annotated_clusters=annotated_clusters, 
        gene_annotation=gene_annotation,
        samples=samples,
        out_dir=pan_genome_folder
    )
    summary_file = output.create_summary(
        rtab_file=rtab_file, 
        out_dir=pan_genome_folder
    )

    # add new protein to existed protein database
    os.system(f'cat {combined_faa_file} >> {protein_database}')

    representative_fasta = output.create_representative_fasta(
        clusters = inflated_clusters,
        gene_annotation = gene_annotation,
        faa_fasta = protein_database,
        out_dir = pan_genome_folder
        )

    json.dump(gene_annotation, open(os.path.join(pan_genome_folder, 'gene_annotation.json'), 'w'), indent=4, sort_keys=True)
    json.dump(gene_position, open(os.path.join(pan_genome_folder, 'gene_position.json'), 'w'), indent=4, sort_keys=True)
    json.dump(inflated_clusters, open(os.path.join(pan_genome_folder, 'unsplit_clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(samples, open(os.path.join(pan_genome_folder, 'samples.json'), 'w'), indent=4, sort_keys=True)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    
    main_cmd = subparsers.add_parser(
        'main',
        description='main pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    main_cmd.set_defaults(func=run_main_pipeline)
    #main_cmd.set_defaults(func=run_main_pipeline_2)
    main_cmd.add_argument('inputs', help='Input files', type=str, nargs='+')
    main_cmd.add_argument('-o', '--out_dir', help='Output directory', required=True, type=str)
    main_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    main_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)
    main_cmd.add_argument('-w', '--over-write', help='over write', default=False, action='store_true')
    main_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    main_cmd.add_argument('-f', '--fasta', help='input fasta file with gff file', default=False, action='store_true')
    main_cmd.add_argument('-d', '--diamond', help='use diamond instead of blastp', default=False, action='store_true')
    main_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=95, type=float)

    add_cmd = subparsers.add_parser(
        'add',
        description='add sample pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_cmd.set_defaults(func=run_add_sample_pipeline)
    add_cmd.add_argument('inputs', help='Input files', type=str, nargs='+')
    add_cmd.add_argument('-c', '--collection-dir', help='Collection directory', required=True, type=str)
    add_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    add_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)
    add_cmd.add_argument('-s', '--dont-split', help='dont split paralog clusters', default=False, action='store_true')
    add_cmd.add_argument('-f', '--fasta', help='input fasta file with gff file', default=False, action='store_true')
    add_cmd.add_argument('-d', '--diamond', help='use diamond instead of blastp', default=False, action='store_true')
    add_cmd.add_argument('-i', '--identity', help='minimum percentage identity', default=95, type=float)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()