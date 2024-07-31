import os
import multiprocessing
import logging
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from datetime import datetime
import pandas as pd
from collections import defaultdict

from panta.utils import *

logger = logging.getLogger(__name__)

def find_paralogs(cluster):#, gene_annotation_dict):
    """
    Find paralogs genes in a cluster. Paralogs genes are
    identified as the genes come from the same sample, that have the
    fewest genes, but greater than 1, in the cluster

    """
    samples = {}
    for gene_id in cluster:
        sample_id,_ = get_seq_ids(gene_id)# gene_annotation_dict[gene_id]['sample_id']
        samples.setdefault(sample_id, []).append(gene_id)

    # pick paralogs with the smallest number of genes
    smallest_number = 1000000
    paralog_genes = None
    for sample_id in samples:
        genes = samples[sample_id]
        count = len(genes)
        if count > 1 and count < smallest_number:
            paralog_genes = genes
            smallest_number = count
    return paralog_genes


# def get_neighbour_genes(gene_annotation_fn, gene_position_fn):
#     gene_neighbour_dict = {}
#     chunksize=100000
#     # Read in gene position
#     # gene_position = {}

#     # with gzip.open(gene_position_fn, 'rt') as gp_fp:
#     #     for line in gp_fp.readlines():
#     #         toks = line.strip().split(',')
#     #         gene_position[(toks[0], toks[1])] = toks[2:]

#     df_it = pd.read_csv(gene_annotation_fn, na_filter= False, index_col='gene_id', usecols=['gene_id', 'gene_index'], chunksize=chunksize)
#     gene_index = 0
#     seq_id = None
#     genes_of_contig = []

#     for chunk_df in df_it:
#         gene_anno_dict = chunk_df.to_dict('index')
#         for gene_id in gene_anno_dict:
#             sample_id, m_seq_id = get_seq_ids(gene_id)
#             if m_seq_id != seq_id:
#                 #when we see a new sequence
#                 for j, _gene_id in enumerate(genes_of_contig):
#                     _start = max(0, j - 5)
#                     gene_neighbour_dict[_gene_id] = genes_of_contig[_start:j] + genes_of_contig[j+1:j+6]
#                 seq_id = m_seq_id
#                 print(f"AAAAA {seq_id}")
#                 gene_neighbour_dict = []
#             genes_of_contig.append(gene_id)
#     #The last sequence
#     for j, _gene_id in enumerate(genes_of_contig):
#         _start = max(0, j - 5)
#         gene_neighbour_dict[_gene_id] = genes_of_contig[_start:j] + genes_of_contig[j+1:j+6]

#             # genes_of_contig = gene_position[(sample_id,seq_id)]
#             # index = genes_of_contig.index(gene_id)
#             # if index != gene_anno_dict[gene_id]['gene_index']:
#             #     gene_index = gene_anno_dict[gene_id]['gene_index']
#             #     raise Exception(f'Index {index} vs {gene_index}')
#             # pre_index = index - 5
#             # post_index = index + 6
#             # if pre_index < 0:
#             #     pre_index = 0
#             # length_of_contig = len(genes_of_contig)
#             # if post_index >= length_of_contig:
#             #     post_index = length_of_contig
#             # neighbour_genes = genes_of_contig[pre_index:index] + genes_of_contig[index+1:post_index]
#             # gene_neighbour_dict[gene_id] = neighbour_genes
#     #TODO check size
#     count_elements = sum([len(gene_neighbour_dict[gene_id]) for gene_id in gene_neighbour_dict])
#     logger.info(f'Size of gene_neighbour_dict = {len(gene_neighbour_dict)} no elements = {count_elements}')

#     return gene_neighbour_dict


def create_orthologs(cluster, paralog_genes, gene_position, gene_to_cluster_index):
    # find cluster indices of all the neighbour genes of each paralog gene
    cluster_indices_around_paralogs = []
    for p in paralog_genes:
        sample_id, seq_id = get_seq_ids(p)
        genes_in_seq = gene_position[(sample_id, seq_id)]
        gene_index = genes_in_seq.index(p)
        neighbours_of_p = genes_in_seq[max(0, gene_index - 5):gene_index] + genes_in_seq[gene_index+1:gene_index+6]
        #neighbours_of_p = gene_neighbour_dict[p]
        cluster_indices_around_p = set()
        for neighbour_gene in neighbours_of_p:
            try:
                cluster_index = gene_to_cluster_index[neighbour_gene]
                cluster_indices_around_p.add(cluster_index)
            except:
                continue
        cluster_indices_around_paralogs.append(cluster_indices_around_p)
    # create data structure to hold new clusters
    new_clusters = [[p] for p in paralog_genes]
    new_clusters.append([]) # extra "leftovers" list to gather genes that don't share CGN with any paralog gene

    # add other members of the cluster to their closest match
    for g in cluster:
        if g in paralog_genes:
            continue

        sample_id, seq_id = get_seq_ids(g)
        genes_in_seq = gene_position[(sample_id, seq_id)]
        gene_index = genes_in_seq.index(g)
        neighbour_genes_of_g = genes_in_seq[max(0, gene_index - 5):gene_index] + genes_in_seq[gene_index+1:gene_index+6]
        #neighbour_genes_of_g = gene_neighbour_dict[g]
        if len(neighbour_genes_of_g) == 0:
            new_clusters[-1].append(g)
            continue

        # find paralog gene which is closest match with g
        best_score = 0
        best_score_index = -1  # -1 is the index of "leftovers" list
        for p,_ in enumerate(paralog_genes):
            cluster_indices_around_p = cluster_indices_around_paralogs[p]
            score_of_p = 0
            for neighbour_gene in neighbour_genes_of_g:
                try:
                    cluster_index = gene_to_cluster_index[neighbour_gene]
                except:
                    continue
                if cluster_index in cluster_indices_around_p:
                    score_of_p += 1
            score_of_p = score_of_p / len(neighbour_genes_of_g)
            if score_of_p > best_score:
                best_score = score_of_p
                best_score_index = p

        new_clusters[best_score_index].append(g)
    # check for "leftovers", remove if absent
    if len(new_clusters[-1]) == 0:
        del new_clusters[-1]

    return new_clusters

def split_paralogs(gene_position_fn, unsplit_clusters, dontsplit):
    if dontsplit == True:
        return unsplit_clusters

    starttime = datetime.now()

    #mem_usage = mem_report(0, "split_paralog0")
    # Read in gene position
    gene_position = {}
    with open(gene_position_fn) as gp_fp:
        for line in gp_fp.readlines():
            toks = line.strip().split(',')
            gene_position[(toks[0], toks[1])] = toks[2:]

    #mem_usage = mem_report(mem_usage, "split_paralog0")
    #TODO: gene_position is memory intensive

    clusters_not_paralogs = set()
    # run iteratively
    out_clusters = unsplit_clusters
    logger.info(f'Number of clusters before spliting {len(out_clusters)}')
    for i in range(100000):
        stime = datetime.now()
        in_clusters = out_clusters
        out_clusters = []
        any_paralogs = 0
        # convert in_clusters so we can find the cluster index of gene
        gene_to_cluster_index = {gene:index for index, genes in enumerate(in_clusters) for gene in genes}

        split_count = 0
        for cluster in in_clusters:
            if len(cluster) == 1:
                out_clusters.append(cluster)
                continue
            first_gene = cluster[0]
            if first_gene in clusters_not_paralogs:
                out_clusters.append(cluster)
                continue

            # check paralogs
            paralog_genes = find_paralogs(cluster)#, gene_annotation_dict)
            # if paralog_genes !=None:
            #     logger.info(f'cluster representative seq {cluster[0]} has '+str(len(paralog_genes)))
            #     for g in paralog_genes:
            #         logger.info(f'{g} ')

            if paralog_genes == None:
                clusters_not_paralogs.add(first_gene)
                out_clusters.append(cluster)
                continue

            # split paralogs
            split_count += 1
            orthologs_clusters = create_orthologs(cluster, paralog_genes, gene_position, gene_to_cluster_index)
            out_clusters.extend(orthologs_clusters)
            any_paralogs = 1

        elapsed = datetime.now() - stime
        logger.info(f'Split paralogs iterate {i} -- count = {split_count}  time taken {str(elapsed)}')
        # check if next iteration is required
        #mem_usage = mem_report(mem_usage, "split_paralog2")

        if any_paralogs == 0:
            break
    split_clusters = out_clusters
    logger.info(f'Number of clusters after spliting {len(out_clusters)}')
    elapsed = datetime.now() - starttime
    logging.info(f'Split paralogs -- time taken {str(elapsed)}')
    #mem_usage = mem_report(mem_usage, "split_paralog3")

    return split_clusters


def annotate_cluster(unlabeled_clusters, gene_annotation_fn):
    starttime = datetime.now()
    clusters = {'groups_' + str(i) : cluster for i, cluster in enumerate(unlabeled_clusters)}
    chunksize = 50000

    annotated_clusters = {}
    suffix = 1
    #gene_annotation_dict = read_csv_to_dict(gene_annotation_fn, 'gene_id', ['gene_name','gene_product'])

    geneid_2_cluster = {}
    for cluster_name in clusters:
        for gene_id in clusters[cluster_name]:
            geneid_2_cluster[gene_id] = cluster_name
    #TODO: we can swap clusters and geneid_2_cluser to save memory
    annotate_cluster_name = defaultdict(dict)
    annotate_cluster_product = defaultdict(set)
    annotate_cluster_rep = {}
    annotate_cluster_len = {} #cluster_id -> (min, max, mean, number)


    df_it = pd.read_csv(gene_annotation_fn, na_filter= False, index_col='gene_id', usecols=['gene_id', 'gene_name','gene_product', 'length'], chunksize=chunksize)
    for chunk_df in df_it:
        gene_anno_dict = chunk_df.to_dict('index')
        for gene_id in gene_anno_dict:
            if gene_id in geneid_2_cluster:
                cluster_name = geneid_2_cluster[gene_id]
                gene_name = gene_anno_dict[gene_id]['gene_name']
                gene_product = gene_anno_dict[gene_id]['gene_product']

                if gene_name:
                    cluster_name_dict = annotate_cluster_name[cluster_name]
                    cluster_name_dict[gene_name] =  cluster_name_dict.setdefault(gene_name, 0) + 1
                if gene_product:
                    annotate_cluster_product[cluster_name].add(gene_product)
                gene_length = gene_anno_dict[gene_id]['length']
                if cluster_name not in annotate_cluster_rep:
                    annotate_cluster_rep[cluster_name] = (gene_id, gene_length)
                    annotate_cluster_len[cluster_name] = (gene_length, gene_length, gene_length,1)
                else:
                    if gene_length > annotate_cluster_rep[cluster_name][1]:
                        annotate_cluster_rep[cluster_name] = (gene_id, gene_length)
                    lens = annotate_cluster_len[cluster_name]
                    annotate_cluster_len[cluster_name] = (
                        min(lens[0], gene_length),
                        max(lens[1], gene_length),
                        (lens[2] * lens[3] + gene_length) / (lens[3] + 1),
                        (lens[3] + 1))
    del geneid_2_cluster
    for cluster_name in clusters:
        gene_id_list = clusters[cluster_name]
        cluster_new_name = cluster_name
        cluster_name_dict = annotate_cluster_name[cluster_name]
        if cluster_name_dict:
            cluster_new_name = max(cluster_name_dict, key=cluster_name_dict.get)

        cluster_product = ', '.join(annotate_cluster_product[cluster_name])
        if not cluster_product:
            cluster_product = 'unknown'

        # check if cluster_new_name already exists
        if cluster_new_name in annotated_clusters:
            cluster_new_name += '_{:05d}'.format(suffix)
            suffix += 1
        lens = annotate_cluster_len[cluster_name]
        annotated_clusters[cluster_new_name] = {
            'gene_id':gene_id_list,
            'product':cluster_product,
            'representative': annotate_cluster_rep[cluster_name][0],
            'min_length': lens[0],
            'max_length': lens[1],
            'mean_length': lens[2],
            'size': lens[3]
            }
    ################
    # for cluster_name in clusters:
    #     cluster_new_name = cluster_name
    #     cluster_product = None
    #     gene_name_count = {}
    #     max_number = 0
    #     gene_id_list = clusters[cluster_name] #TODO: check if set is better than list
    #     for gene_id in gene_id_list:
    #         gene_name = gene_annotation_dict[gene_id]['gene_name']
    #         gene_product = gene_annotation_dict[gene_id]['gene_product']
    #         if gene_name:
    #             gene_name_count[gene_name] = gene_name_count.get(gene_name, 0) + 1
    #             if gene_name_count[gene_name] > max_number:
    #                 cluster_new_name = gene_name
    #                 max_number = gene_name_count[gene_name]
    #                 if gene_product:
    #                     cluster_product = gene_product

    #     if cluster_product == None:
    #         cluster_product =[] #TODO: check if set is better than list
    #         for gene_id in gene_id_list:
    #             gene_name = gene_annotation_dict[gene_id]['gene_name']
    #             gene_product = gene_annotation_dict[gene_id]['gene_product']

    #             if gene_product:
    #                 if gene_product not in cluster_product:
    #                     cluster_product.append(gene_product)
    #         if len(cluster_product) > 0:
    #             cluster_product = ', '.join(cluster_product)
    #         else:
    #             cluster_product = 'unknown'
    #     # check if cluster_new_name already exists
    #     if cluster_new_name in annotated_clusters:
    #         cluster_new_name += '_{:05d}'.format(suffix)
    #         suffix += 1
    #    annotated_clusters[cluster_new_name] = {'gene_id':gene_id_list, 'product':cluster_product}

    elapsed = datetime.now() - starttime
    logging.info(f'Annotate clusters -- time taken {str(elapsed)}')
    return annotated_clusters


def create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, out_dir, add=False):
    starttime = datetime.now()
    clusters_dir = os.path.join(out_dir, 'clusters')
    pan_ref_file = os.path.join(out_dir, 'pan_genome_reference.fna')
    with open(pan_ref_file, 'w') as ref_fh:
        for sample in samples:
            sample_id = sample['id']
            fna_file = os.path.join(out_dir, 'samples', sample_id, sample_id + '.fna')
            with open(fna_file) as fna_fh:
                for seq in SeqIO.parse(fna_fh, 'fasta'):
                    seq_id = seq.id
                    if seq_id in gene_to_cluster_name:
                        cluster_name = gene_to_cluster_name[seq_id]
                        seq_fasta = SeqIO.FastaIO.as_fasta(seq)
                        if add:
                            gene_add_file_in_cluster=os.path.join(clusters_dir, cluster_name, cluster_name + '.add.fna')
                            with open(gene_add_file_in_cluster, 'a') as out_fh:
                                out_fh.write(seq_fasta)
                        gene_file_in_cluster=os.path.join(clusters_dir, cluster_name, cluster_name + '.fna')
                
                        with open(gene_file_in_cluster, 'a') as out_fh:
                            out_fh.write(seq_fasta)
                        if seq_id in pan_ref_list:
                            ref_fh.write(seq_fasta)
                    else:
                        #logger.error(f'Gene {seq_id} not in clusters? why')
                        pass

    elapsed = datetime.now() - starttime
    logging.info(f'Create nucleotide sequence file for each gene cluster -- time taken {str(elapsed)}')

def create_pro_file_for_each_cluster(samples, gene_to_cluster_name, out_dir,add=False):
    starttime = datetime.now()

    for sample in samples:
        sample_id = sample['id']
   
        faa_file = os.path.join(out_dir, 'samples', sample_id, sample_id + '.faa')
        with open(faa_file) as faa_fh:
            for seq in SeqIO.parse(faa_fh, 'fasta'):
                if seq.id not in gene_to_cluster_name:
                    #logger.error(f'Gene {seq.id} not in clusters? why')
                    continue
                cluster_name = gene_to_cluster_name[seq.id]
                if add:
                    gene_add_file_in_cluster=os.path.join(out_dir, 'clusters', cluster_name, cluster_name + '.add.faa')
                    with open(gene_add_file_in_cluster, 'a') as out_fh:
                        out_fh.write(SeqIO.FastaIO.as_fasta(seq))
                gene_file_in_cluster=os.path.join(out_dir, 'clusters', cluster_name, cluster_name + '.faa')
                with open(gene_file_in_cluster, 'a') as out_fh:
                    out_fh.write(SeqIO.FastaIO.as_fasta(seq))

    elapsed = datetime.now() - starttime
    logging.info(f'Create protein sequence file for each gene cluster -- time taken {str(elapsed)}')


def run_mafft_protein_alignment(annotated_clusters, out_dir, threads=1):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, 'clusters')

    #cmds_file = os.path.join(clusters_dir,"pro_align_cmds")
    pool = multiprocessing.Pool(processes=threads)
    results = []
    #with open(cmds_file,'w') as cmds:
    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)
        gene_aln_file = os.path.join(cluster_dir, cluster_name + '.faa.aln')
        gene_seq_file = os.path.join(cluster_dir, cluster_name + '.faa')
        if not os.path.isfile(gene_seq_file):
            logger.info('{} does not exist'.format(gene_seq_file))
            continue

        #check if the alignment has been done
        if os.path.isfile(gene_aln_file):
            inset = set()
            with open(gene_seq_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    inset.add(seq.id)

            outset = set()
            with open(gene_aln_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    outset.add(seq.id)

            if inset == outset:
                logger.info(f'Aligment for {cluster_name} exists, skip realignment')
                continue

        cmd = f"mafft --auto --quiet --thread 1 {gene_seq_file} > {gene_aln_file}"
        #cmds.write(cmd + '\n')
        results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running maftt')

    elapsed = datetime.now() - starttime
    logging.info(f'Run protein alignment -- time taken {str(elapsed)}')

def run_poa_protein_alignment(annotated_clusters, out_dir, threads=1):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, 'clusters')

    #cmds_file = os.path.join(clusters_dir,"pro_align_cmds")
    pool = multiprocessing.Pool(processes=threads)
    results = []
    #with open(cmds_file,'w') as cmds:
    cmd=''
    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)
        gene_aln_file = os.path.join(cluster_dir, cluster_name + '.faa.aln')
        gene_graph_file = os.path.join(cluster_dir, cluster_name + '.faa.gfa')
        gene_seq_file = os.path.join(cluster_dir, cluster_name + '.faa')
        if not os.path.isfile(gene_seq_file):
            logger.info('{} does not exist'.format(gene_seq_file))
            continue

        #check if the alignment has been done or need to progressive
        addset=set()
        if os.path.isfile(gene_aln_file):
            inset = set()
            with open(gene_seq_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    inset.add(seq.id)

            outset = set()
            with open(gene_aln_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    outset.add(seq.id)

            if inset == outset:
                logger.info(f'Aligment for {cluster_name} exists, skip realignment')
                continue
            else:
                addset=inset.difference(outset)
            print('len of '+gene_aln_file+' is '+str(len(outset))+', len of '+gene_seq_file+' is '+str(len(inset)))
        if len(addset) > 0 and len(outset) > 1:
            #progressive
            gene_seq_file = os.path.join(cluster_dir, cluster_name + '.add.faa')
            cmd= f"abpoa -c -i {gene_aln_file} {gene_seq_file} -r1 > {gene_aln_file}" 
        else:    
            cmd = f"abpoa  {gene_seq_file} -c -r1 > {gene_aln_file}"
        #cmds.write(cmd + '\n')
        results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running poa')

    elapsed = datetime.now() - starttime
    logging.info(f'Run protein POA alignment -- time taken {str(elapsed)}')

def run_hhm_profile_clusters(annotated_clusters, out_dir, threads=1):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, 'clusters')

    #cmds_file = os.path.join(clusters_dir,"pro_align_cmds")
    pool = multiprocessing.Pool(processes=threads)
    results = []
    #with open(cmds_file,'w') as cmds:
    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)
        gene_aln_file = os.path.join(cluster_dir, cluster_name + '.faa.aln')
        gene_hmm_file = os.path.join(cluster_dir, cluster_name + '.hmm')
        gene_seq_file = os.path.join(cluster_dir, cluster_name + '.faa')
        if not os.path.isfile(gene_aln_file):
            logger.info('{} does not exist'.format(gene_aln_file))
            continue
        if not os.path.isfile(gene_seq_file):
            logger.info('{} does not exist'.format(gene_seq_file))
            continue
        #check if the alignment has been done
        if os.path.isfile(gene_aln_file):
            inset = set()
            with open(gene_seq_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    inset.add(seq.id)

            outset = set()
            with open(gene_aln_file) as fh:
                for seq in Se-iqIO.parse(fh, 'fasta'):
                    outset.add(seq.id)

            if inset == outset and os.path.isfile(gene_hmm_file) :
                logger.info(f'HMM profile for {cluster_name} exists, skip hmm build profile')
                continue

        cmd = f"hmmbuild --amino  --cpu 1  {gene_hmm_file} {gene_aln_file}"
        #cmds.write(cmd + '\n')
        results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running poa')

    elapsed = datetime.now() - starttime
    logging.info(f'Run protein hmmbuild -- time taken {str(elapsed)}')

def run_mafft_nucleotide_alignment(annotated_clusters, out_dir, threads=1):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, 'clusters')

    #cmds_file = os.path.join(clusters_dir,"nu_align_cmds")
    pool = multiprocessing.Pool(processes=threads)
    results = []
    #with open(cmds_file,'w') as cmds:
    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)
        gene_aln_file = os.path.join(cluster_dir, cluster_name + '.fna.aln')
        gene_seq_file = os.path.join(cluster_dir, cluster_name + '.fna')
        if not os.path.isfile(gene_seq_file):
            logger.info('{} does not exist'.format(gene_aln_file))
            continue

        #check if the alignment has been done
        addset=set()
        if os.path.isfile(gene_aln_file):
            inset = set()
            with open(gene_seq_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    inset.add(seq.id)

            outset = set()
            with open(gene_aln_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    outset.add(seq.id)

            if inset == outset:
                logger.info(f'Aligment for {cluster_name} exists, skip realignment')
                continue
            else:
                addset=inset.difference(outset)
        if len(addset)>0 and len(outset)>1:
            gene_seq_file=os.path.join(cluster_dir, cluster_name + '.add.faa')
            cmd= f"mafft --add {gene_seq_file} --reorder {gene_aln_file} > {gene_aln_file}"
        else:
            cmd = f"mafft --auto --quiet --thread 1 {gene_seq_file} > {gene_aln_file}"
            #cmd += f' && rm {gene_seq_file}'
        #cmds.write(cmd + '\n')
        results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running maftt 2')

    elapsed = datetime.now() - starttime
    logging.info(f'Run nucleotide alignment -- time taken {str(elapsed)}')
def run_poa_nucleotide_alignment(annotated_clusters, out_dir, threads=1):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, 'clusters')

    #cmds_file = os.path.join(clusters_dir,"nu_align_cmds")
    pool = multiprocessing.Pool(processes=threads)
    results = []
    #with open(cmds_file,'w') as cmds:
    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)
        gene_aln_file = os.path.join(cluster_dir, cluster_name + '.fna.aln')
        gene_seq_file = os.path.join(cluster_dir, cluster_name + '.fna')
        if not os.path.isfile(gene_seq_file):
            logger.info('{} does not exist'.format(gene_aln_file))
            continue

        #check if the alignment has been done
        addset=set()
        if os.path.isfile(gene_aln_file):
            inset = set()
            with open(gene_seq_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    inset.add(seq.id)

            outset = set()
            with open(gene_aln_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    outset.add(seq.id)

            if inset == outset:
                logger.info(f'Aligment for {cluster_name} exists, skip realignment')
                continue
            else:
                addset=inset.difference(outset)
        if len(addset) > 0 and len(outset) > 1:
            #progressive
            gene_seq_file = os.path.join(cluster_dir, cluster_name + '.add.fna')
            cmd= f"abpoa -i {gene_aln_file} {gene_seq_file} -r1 > {gene_aln_file}" 
        else: 
            cmd = f"abpoa {gene_seq_file} -r1 > {gene_aln_file}"
        #cmd += f' && rm {gene_seq_file}'
        #cmds.write(cmd + '\n')
        results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running maftt 2')

    elapsed = datetime.now() - starttime
    logging.info(f'Run nucleotide alignment by poa -- time taken {str(elapsed)}')

def create_nucleotide_alignment(annotated_clusters, out_dir):
    starttime = datetime.now()

    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(out_dir, 'clusters', cluster_name)

        protein_aln_file = os.path.join(cluster_dir, cluster_name + '.faa.aln')
        if not os.path.isfile(protein_aln_file):
            logger.info('{} does not exist'.format(protein_aln_file))
            continue
        protein_dict = {}
        with open(protein_aln_file, 'rt') as fh:
            for seq_record in SeqIO.parse(fh, 'fasta'):
                protein_dict[seq_record.id] = str(seq_record.seq)

        nucleotide_seq_file = os.path.join(cluster_dir, cluster_name + '.fna')
        if not os.path.isfile(nucleotide_seq_file):
            logger.info('{} does not exist'.format(nucleotide_seq_file))
            continue
        nucleotide_dict = {}
        for seq_record in SeqIO.parse(nucleotide_seq_file, 'fasta'):
            nucleotide_dict[seq_record.id] = str(seq_record.seq)

        nucleotide_aln_file = os.path.join(cluster_dir, cluster_name + '.fna.aln')
        with open(nucleotide_aln_file, 'wt') as fh:
            for seq_id in protein_dict.keys():
                protein = protein_dict[seq_id]
                nucleotide = nucleotide_dict[seq_id]
                result = ''
                codon_pos = 0
                for c in protein:
                    if c == '-':
                        result += '---'
                    else:
                        result += nucleotide[codon_pos * 3: codon_pos * 3 + 3]
                        codon_pos += 1
                new_record = SeqRecord(Seq(result), id = seq_id, description = '')
                SeqIO.write(new_record, fh, 'fasta')

        #os.remove(nucleotide_seq_file)
        #os.remove(os.path.join(cluster_dir, cluster_name + '.faa'))

    elapsed = datetime.now() - starttime
    logging.info(f'Create  nucleotide alignment -- time taken {str(elapsed)}')

def create_core_gene_alignment(annotated_clusters,
                            #gene_annotation_dict,
                            samples, out_dir):
    starttime = datetime.now()
    clusters_dir = os.path.join(out_dir, 'clusters')
    seq_dict = {}
    for sample in samples:
        sample_id = sample['id']
        seq_dict[sample_id] = ''

    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)

        sample_list = set()
        skip = False
        for gene_id in annotated_clusters[cluster_name]['gene_id']:
            sample_id, _ = get_seq_ids(gene_id)
            #sample_id = gene_annotation_dict[gene_id]['sample_id']

            #length = gene_annotation_dict[gene_id]['length']
            if sample_id not in sample_list:
                sample_list.add(sample_id)
            else:
                skip = True
        if len(sample_list) < len(samples):
            skip = True

        if skip == True:
            continue

        nucleotide_aln_file = os.path.join(cluster_dir, cluster_name + '.fna.aln')
        if not os.path.isfile(nucleotide_aln_file):
            logger.info('{} does not exist'.format(nucleotide_aln_file))
            continue

        cluster_dict = {}
        with open(nucleotide_aln_file, 'rt') as fh:
            for seq_record in SeqIO.parse(fh, 'fasta'):
                sample_id,_ = get_seq_ids(seq_record.id)
                #sample_id = gene_annotation_dict[seq_record.id]['sample_id']
                cluster_dict[sample_id] = str(seq_record.seq)

        for sample_id in cluster_dict:
            seq_dict[sample_id] += cluster_dict[sample_id]

    core_gene_aln_file = os.path.join(out_dir, 'core_gene_alignment.aln.gz')
    with gzip.open(core_gene_aln_file, 'wt') as fh:
        for sample in seq_dict:
            new_record = SeqRecord(Seq(seq_dict[sample]), id = sample, description = '')
            SeqIO.write(new_record, fh, 'fasta')

    elapsed = datetime.now() - starttime
    logging.info(f'Create core gene alignment -- time taken {str(elapsed)}')

def run_gene_alignment(annotated_clusters, samples,oldsamples, collection_dir, alignment, coverage_threshold=0.0, threads=1,poa=False, add=False):
    count_threshold = int(len(samples) * coverage_threshold)

    #Need to have at least 2 genes before we can do alignment
    if count_threshold < 2:
        count_threshold = 2

    gene_to_cluster_name = {}
    pan_ref_list = set()

    clusters_dir = os.path.join(collection_dir, 'clusters')
    # if os.path.exists(clusters_dir):
    #     shutil.rmtree(clusters_dir)
    #Create a fresh
    if not os.path.exists(clusters_dir):
        os.mkdir(clusters_dir)

    clusters_to_align = []
    #gene_annotation_dict = read_csv_to_dict(gene_annotation_fn, 'gene_id', ['sample_id','length'])
    for cluster_name in annotated_clusters:
        # if len(annotated_clusters[cluster_name]['gene_id']) < count_threshold:
        #     continue
        clusters_to_align.append(cluster_name)

        cluster_dir = os.path.join(collection_dir, 'clusters', cluster_name)
        if not os.path.exists(cluster_dir):
            os.mkdir(cluster_dir)
        #length_max = 0
        #representative = None
        for gene_id in annotated_clusters[cluster_name]['gene_id']:
            gene_to_cluster_name[gene_id] = cluster_name
            #sample_id = gene_annotation_dict[gene_id]['sample_id']
            #length = gene_annotation_dict[gene_id]['length']
            #if length > length_max:
            #    representative = gene_id
            #    length_max = length
        pan_ref_list.add(annotated_clusters[cluster_name]['representative'])

        #Remove folders that are not in clusters to align (in case the clusters changes)
    existing_folders = os.listdir(clusters_dir)
    for folder_name in existing_folders:
        if folder_name not in clusters_to_align:
            folder_path = os.path.join(clusters_dir,folder_name)
            if os.path.isdir(folder_path):
                shutil.rmtree(folder_path)
                logger.info(f'Clean up cluster {folder_name}')            

    if 'protein' == alignment:
        create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, collection_dir,add=add)
        create_pro_file_for_each_cluster(samples, gene_to_cluster_name, collection_dir,add=add)
        if poa==True:
            run_poa_protein_alignment(clusters_to_align, collection_dir, threads=threads)
        else:
            run_mafft_protein_alignment(clusters_to_align, collection_dir, threads=threads)
        
        #run_hhm_profile_clusters(clusters_to_align,collection_dir,  threads=threads)
        create_nucleotide_alignment(clusters_to_align, collection_dir)
        

    if 'nucleotide' == alignment:
        create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, collection_dir)
        if poa==True:
            run_poa_nucleotide_alignment(clusters_to_align, collection_dir, threads=threads)
        else:
            run_mafft_nucleotide_alignment(clusters_to_align, collection_dir, threads=threads)

    create_core_gene_alignment(annotated_clusters,oldsamples,collection_dir)

def create_poa_protein_consensus(annotated_clusters, out_dir,threads=1):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, 'clusters')

    #cmds_file = os.path.join(clusters_dir,"pro_align_cmds")
    pool = multiprocessing.Pool(processes=threads)
    results = []
    #with open(cmds_file,'w') as cmds:
    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)
        gene_aln_file = os.path.join(cluster_dir, cluster_name + '.faa.aln')
        #gene_hmm_file = os.path.join(cluster_dir, cluster_name + '.hmm')
        gene_seq_file = os.path.join(cluster_dir, cluster_name + '.faa')
        
        gene_seq_consensus_file = os.path.join(cluster_dir, cluster_name + '.consensus.faa')
        #if not os.path.isfile(gene_aln_file):
        #    logger.info('{} does not exist'.format(gene_aln_file))
        #    continue
        if not os.path.isfile(gene_seq_file):
            logger.info('{} does not exist'.format(gene_seq_file))
            continue
        #check if the alignment has been done
        if os.path.isfile(gene_aln_file):
            inset = set()
            with open(gene_seq_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    inset.add(seq.id)

            outset = set()
            with open(gene_aln_file) as fh:
                for seq in SeqIO.parse(fh, 'fasta'):
                    outset.add(seq.id)

            if inset == outset and os.path.isfile(gene_seq_consensus_file) :
                logger.info(f'Consensus for {cluster_name} exists, skip consensus by poa')
                continue

        cmd = f"abpoa -c -t BLOSUM62.mtx {gene_seq_file} > {gene_seq_consensus_file} && sed -i 's/^>Consensus_sequence/>{cluster_name}/' {gene_seq_consensus_file}"
        results.append(pool.apply_async(run_command,(cmd, None)))
        
        annotated_clusters[cluster_name]['representative']=gene_seq_consensus_file
        #cmds.write(cmd + '\n')
        
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running poa consensus')
    
    elapsed = datetime.now() - starttime
    logging.info(f'Run create protein consensus by poa -- time taken {str(elapsed)}')
def make_protein_consensus_db(annotated_clusters, out_dir,method='diamond',threads=1):
    consensus_fasta = os.path.join(out_dir, 'consensus.fasta')
    for cluster_name in annotated_clusters:
        if os.path.isfile(annotated_clusters[cluster_name]["representative"]):
            cmd = f'cat {annotated_clusters[cluster_name]["representative"]} >> {consensus_fasta}'
            ret = run_command(cmd)
            if ret != 0:
                logger.info(f'Concat consensus of {cluster_name} failed!')
    if method=='diamond':
        diamond_db = os.path.join(out_dir, 'consenus_diamond_db')
        cmd = f'./diamond makedb --in {consensus_fasta} -d {diamond_db} -p {threads} --quiet'
        #ret = os.system(cmd)
        ret = run_command(cmd)
        if ret != 0:
            raise Exception('Error running diamond makedb for consensus file')
    if method=='mmseq':
        mmseq_db = os.path.join(out_dir, 'consensus_mmseq_db')
        cmd = f'mmseqs createdb {consensus_fasta} {mmseq_db}'
        #ret = os.system(cmd)
        ret = run_command(cmd)
        if ret != 0:
            raise Exception('Error running mmseqs createdb consensus file')
def merge_new_cluster_to_old_clusters(new_annotated_clusters,old_annotated_clusters):
    count={}
    for clustername in new_annotated_clusters:
        if clustername in old_annotated_clusters:
            #search max suffix
            suffix=1
            while clustername+'_'+str(suffix) in old_annotated_clusters:
                suffix=suffix+1
            
            old_annotated_clusters[clustername+'_'+str(suffix)]=new_annotated_clusters[clustername]
    return old_annotated_clusters