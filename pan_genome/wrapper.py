import os
import subprocess
import shutil
import logging
import json
from pan_genome import *
import pan_genome.post_analysis as pa

logger = logging.getLogger(__name__)

def add_sample(new_samples, old_represent_faa, old_clusters, gene_to_old_cluster, collection_dir, temp_dir, args):
    
    gene_dictionary = data_preparation.extract_proteins(new_samples,collection_dir,args)
    
    new_combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=temp_dir,
        samples=new_samples)

    unmatched_faa, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1 = old_represent_faa,
        database_2 = new_combined_faa,
        out_dir = temp_dir,
        threads=args.threads)

    add_sample_pipeline.add_gene_cd_hit_2d(old_clusters, cd_hit_2d_clusters, gene_to_old_cluster)

    num_seq = subprocess.run(f'grep ">" {unmatched_faa} | wc -l', shell=True, capture_output=True, text=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return None, None

    unmatched_represent_faa, unmatched_clusters = main_pipeline.run_cd_hit(
        faa_file=unmatched_faa,
        out_dir=temp_dir,
        threads=args.threads)

    blast_1_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = old_represent_faa,
        query_fasta = unmatched_represent_faa,
        out_dir = os.path.join(temp_dir, 'blast1'),
        evalue = args.evalue,
        max_target_seqs=2000,
        threads=args.threads)

    remain_fasta, old_clusters = add_sample_pipeline.add_gene_blast(
        old_clusters=old_clusters,
        gene_to_cluster=gene_to_old_cluster,
        unmatched_clusters = unmatched_clusters,
        blast_result=blast_1_result, 
        fasta_file=unmatched_represent_faa, 
        out_dir=temp_dir,
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    num_seq = subprocess.run(f'grep ">" {remain_fasta} | wc -l', shell=True, capture_output=True, text=True)
    if int(num_seq.stdout.rstrip()) == 0:
        return None, None
    
    blast_2_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = remain_fasta,
        query_fasta = remain_fasta,
        out_dir = os.path.join(temp_dir, 'blast2'),
        evalue = args.evalue,
        max_target_seqs=2000,
        threads=args.threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_2_result,
        out_dir = temp_dir, 
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result)

    new_clusters, gene_to_new_cluster = main_pipeline.reinflate_clusters(
        cd_hit_clusters = unmatched_clusters,
        mcl_file=mcl_file,
        gene_dictionary = gene_dictionary
        )

    return new_clusters, gene_to_new_cluster, unmatched_represent_faa, gene_dictionary


def run_main_pipeline(samples, collection_dir, temp_dir, db_dir, args):
    
    gene_dictionary = data_preparation.extract_proteins(samples,collection_dir,args)

    combined_faa = data_preparation.combine_proteins(
        collection_dir= collection_dir, 
        out_dir=temp_dir,
        samples=samples)

    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=combined_faa,
        out_dir=temp_dir,
        threads=args.threads)

    blast_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta = cd_hit_represent_fasta,
        query_fasta = cd_hit_represent_fasta,
        out_dir = os.path.join(temp_dir, 'blast'),
        evalue = args.evalue,
        max_target_seqs=2000,
        threads=args.threads)

    filtered_blast_result = main_pipeline.filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir, 
        identity=args.identity, LD=args.LD, AS=args.AS, AL=args.AL)

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir = temp_dir,
        blast_result = filtered_blast_result)

    clusters, gene_to_cluster = main_pipeline.reinflate_clusters(
        cd_hit_clusters=cd_hit_clusters,
        mcl_file=mcl_file,
        gene_dictionary=gene_dictionary)
 
    if args.fasta == None:
        clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=clusters, 
            gene_dictionary=gene_dictionary,
            start = 1)
    else:
        # create representative
        represent_fasta = os.path.join(temp_dir, 'representative.fasta')
        utils.create_fasta_include(
            fasta_file_list=[cd_hit_represent_fasta], 
            include_list=gene_to_cluster, 
            output_file=represent_fasta
            )
        
        clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=clusters,
            rep_fasta = represent_fasta,
            temp_dir=temp_dir,
            db_dir = db_dir,
            threads = args.threads,
            start = 1
            )

    output.create_representative_fasta(
        gene_to_cluster=gene_to_cluster, 
        clusters_annotation=clusters_annotation, 
        source_fasta=cd_hit_represent_fasta, 
        out_fasta=os.path.join(collection_dir, 'representative.fasta'), 
        mode='w')

    output.create_spreadsheet(clusters, clusters_annotation, gene_dictionary, samples, collection_dir)
    rtab_file = output.create_rtab(clusters, clusters_annotation, gene_dictionary,samples,collection_dir)
    output.create_summary(rtab_file, collection_dir)


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