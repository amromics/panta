import os
import subprocess
import shutil
import logging
import json
from pan_genome import *
import pan_genome.post_analysis as pa

logger = logging.getLogger(__name__)

def add_sample(new_samples, old_represent_faa, old_clusters, gene_dictionary, temp_dir, collection_dir, threads, args):
    new_combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=temp_dir,
        samples=new_samples)

    unmatched_faa, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1 = old_represent_faa,
        database_2 = new_combined_faa,
        out_dir = temp_dir,
        threads=threads)

    gene_to_cluster, old_clusters = add_sample_pipeline.add_gene_cd_hit_2d(old_clusters, cd_hit_2d_clusters)

    num_seq = subprocess.run(f'grep ">" {unmatched_faa} | wc -l', shell=True, capture_output=True, text=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return old_clusters, new_combined_faa

    unmatched_represent_faa, unmatched_clusters = main_pipeline.run_cd_hit(
        faa_file=unmatched_faa,
        out_dir=temp_dir,
        threads=threads)

    blast_1_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = old_represent_faa,
        query_fasta = unmatched_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast1'),
        evalue = args.evalue,
        threads=threads)

    remain_fasta, old_clusters = add_sample_pipeline.add_gene_blast(
        old_clusters=old_clusters,
        gene_to_cluster=gene_to_cluster,
        unmatched_clusters = unmatched_clusters,
        blast_result=blast_1_result, 
        fasta_file=unmatched_represent_faa, 
        out_dir=temp_dir,
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    num_seq = subprocess.run(f'grep ">" {remain_fasta} | wc -l', shell=True, capture_output=True, text=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return old_clusters, new_combined_faa
    
    blast_2_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = remain_fasta,
        query_fasta = remain_fasta,
        out_dir = os.path.join(temp_dir, 'blast2'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_2_result,
        out_dir = temp_dir, 
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result)

    new_clusters, new_represent_list = add_sample_pipeline.create_new_clusters(
        unmatched_clusters = unmatched_clusters,
        mcl_file=mcl_file,
        gene_dictionary = gene_dictionary
        )

    # create representative
    new_represent_fasta = os.path.join(temp_dir, 'representative.fasta')
    utils.create_fasta_include(
        fasta_file_list=[unmatched_represent_faa], 
        include_list=new_represent_list, 
        output_file=new_represent_fasta
        ) 


    return new_clusters, new_represent_fasta


def run_main_pipeline(samples, gene_dictionary, gene_position, collection_dir, temp_dir, db_dir, args, anno, threads):

    data_preparation.extract_proteins(samples,collection_dir,gene_dictionary,gene_position,args.table,anno,threads)

    # split collection
    number = args.number
    if number == 0:
        subset = samples[0:]
        remain = []
    else:
        subset = samples[0:number]
        remain = samples[number:]

    # run a subset of collection
    subset_dir = os.path.join(temp_dir, 'subset')
    if not os.path.exists(subset_dir):
        os.makedirs(subset_dir)
    
    subset_combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=subset_dir,
        samples=subset)

    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=subset_combined_faa,
        out_dir=subset_dir,
        threads=threads)

    blast_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(subset_dir, 'blast'),
        evalue = args.evalue,
        threads=threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = subset_dir, 
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    subset_mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = subset_dir,
        blast_result = filtered_blast_result)

    subset_inflated_clusters, subset_represent_list = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_clusters,
        mcl_file=subset_mcl_file,
        gene_dictionary=gene_dictionary)

    # create representative
    subset_represent_fasta = os.path.join(subset_dir, 'representative.fasta')
    utils.create_fasta_include(
        fasta_file_list=[cd_hit_represent_fasta], 
        include_list=subset_represent_list, 
        output_file=subset_represent_fasta
        )


    # run the remain of collection
    if len(remain) == 0:
        final_clusters = subset_inflated_clusters
        final_represent_fasta = os.path.join(collection_dir, 'representative.fasta')
        shutil.move(subset_represent_fasta,final_represent_fasta)
    else:
        remain_dir = os.path.join(temp_dir, 'remain')
        if not os.path.exists(remain_dir):
            os.makedirs(remain_dir)

        remain_clusters, remain_represent_fasta = add_sample(
            new_samples=remain, 
            old_represent_faa=subset_represent_fasta,
            old_clusters=subset_inflated_clusters,
            gene_dictionary=gene_dictionary,
            temp_dir=remain_dir, 
            collection_dir=collection_dir, 
            threads=threads, 
            args=args)

        subset_inflated_clusters.extend(remain_clusters)
        final_clusters = subset_inflated_clusters
        final_represent_fasta = os.path.join(collection_dir, 'representative.fasta')
        os.system(f'cat {subset_represent_fasta} {remain_represent_fasta} > {final_represent_fasta}')


    # split_clusters = post_analysis.split_paralogs(
    #     gene_dictionary=gene_dictionary,
    #     gene_position=gene_position,
    #     unsplit_clusters= inflated_clusters,
    #     split=args.split
    #     )
    
 
    if anno == False:
        clusters_annotation = post_analysis.annotate_cluster(
            unlabeled_clusters=final_clusters, 
            gene_dictionary=gene_dictionary,
            start = 1)
    else:
        clusters_annotation = annotate.annotate_cluster(
            unlabeled_clusters=final_clusters,
            rep_fasta = final_represent_fasta,
            temp_dir=temp_dir,
            db_dir = db_dir,
            threads = threads,
            start = 1
            )

    output.create_spreadsheet(final_clusters, clusters_annotation, gene_dictionary, samples, collection_dir)
    rtab_file = output.create_rtab(final_clusters, clusters_annotation, gene_dictionary,samples,collection_dir)
    output.create_summary(rtab_file, collection_dir)

    output.write_gene_dictionary(gene_dictionary, collection_dir, 'w')
    output.write_gene_position(gene_position, collection_dir, 'w')
    json.dump(samples, open(os.path.join(collection_dir, 'samples.json'), 'w'), indent=4, sort_keys=True)

    json.dump(final_clusters, open(os.path.join(collection_dir, 'clusters.json'), 'w'), indent=4, sort_keys=True)
    json.dump(clusters_annotation, open(os.path.join(collection_dir, 'clusters_annotation.json'), 'w'), indent=4, sort_keys=True)

    return final_clusters


def run_gene_alignment(annotated_clusters, gene_dictionary, samples, collection_dir, alignment, threads):
    
    if alignment == None:
        return

    gene_to_cluster_name = {}
    pan_ref_list = set()

    clusters_dir = os.path.join(collection_dir, 'clusters')
    if os.path.exists(clusters_dir):
        shutil.rmtree(clusters_dir)
        os.mkdir(clusters_dir)
    else:
        os.mkdir(clusters_dir)

    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(collection_dir, 'clusters', cluster_name)
        if not os.path.exists(cluster_dir):
            os.mkdir(cluster_dir)
        length_max = 0
        representative = None
        for gene in annotated_clusters[cluster_name]['gene_id']:
            gene_to_cluster_name[gene] = cluster_name
            length = gene_dictionary[gene][2]
            if length > length_max:
                representative = gene
                length_max = length
        pan_ref_list.add(representative)

    if alignment == 'protein':
        pa.create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, collection_dir)
        pa.create_pro_file_for_each_cluster(samples, gene_to_cluster_name, collection_dir)
        pa.run_mafft_protein_alignment(annotated_clusters, collection_dir, threads)
        pa.create_nucleotide_alignment(annotated_clusters, collection_dir)
    if alignment == 'nucleotide':
        pa.create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, collection_dir)
        pa.run_mafft_nucleotide_alignment(annotated_clusters, collection_dir, threads)
    
    pa.create_core_gene_alignment(annotated_clusters, gene_dictionary,samples,collection_dir)