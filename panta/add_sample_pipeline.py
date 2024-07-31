import os
import logging
import copy
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from panta.utils import run_command, parse_cluster_file

logger = logging.getLogger(__name__)

def run_cd_hit_2d(database_1, database_2, out_dir, threads=4):
    starttime = datetime.now()

    not_match_fasta = os.path.join(out_dir, 'cd-hit-2d.fasta')
    cd_hit_cluster_file = not_match_fasta + '.clstr'
    
    cmd = f'cd-hit-2d -i {database_1} -i2 {database_2} -o {not_match_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    #cdhit_log = f'time_cdhit2d_{datetime.timestamp(starttime)}.log'
    #ret = run_command(cmd, cdhit_log)
    ret = os.system(cmd)
    if ret != 0:
        raise Exception('Error running cd-hit-2d')

    clusters = parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT-2D with 98% identity -- time taken {str(elapsed)}')
    return not_match_fasta, clusters


def combine_blast_results(blast_1, blast_2, blast_3, outdir):
    combined_blast_results = os.path.join(outdir, 'combined_blast_results')

    #os.system(f'cat {blast_1}  > {combined_blast_results}')
    os.system(f'cat {blast_1} {blast_2} {blast_3} > {combined_blast_results}')
    os.remove(blast_2)
    os.remove(blast_3)
    

    return combined_blast_results


def combine_representative(new, old, out_dir):
    temp_file = os.path.join(out_dir, 'representative_temp')
    out_file = os.path.join(out_dir, 'representative.fasta')
    os.system(f'cat {old} {new} > {temp_file}')

    os.replace(temp_file, out_file)


def reinflate_clusters(old_clusters, cd_hit_2d_clusters, not_match_clusters, mcl_file):
    starttime = datetime.now()

    # clusters for next run
    new_clusters = copy.deepcopy(old_clusters)
    for gene in new_clusters:
        new_clusters[gene].extend(cd_hit_2d_clusters[gene])
    new_clusters.update(copy.deepcopy(not_match_clusters))

    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                if gene in old_clusters:
                    inflated_genes.extend(old_clusters[gene])
                    inflated_genes.extend(cd_hit_2d_clusters[gene])
                    del old_clusters[gene]
                    del cd_hit_2d_clusters[gene]
                if gene in not_match_clusters:
                    inflated_genes.extend(not_match_clusters[gene])
                    del not_match_clusters[gene]
            inflated_clusters.append(inflated_genes)
    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in old_clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(old_clusters[gene])
        inflated_genes.extend(cd_hit_2d_clusters[gene])
        inflated_clusters.append(inflated_genes)

    for gene in not_match_clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(not_match_clusters[gene])

    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate clusters -- time taken {str(elapsed)}')
    return inflated_clusters, new_clusters

def match_new_sequence(new_seqs_file,old_clusters,out_dir,consensusdb,method='diamond',evalue=1E-6,threads=1):
    #blast with diamond
    starttime = datetime.now()
    if method=='diamond':
        diamond_result = os.path.join(out_dir, 'diamond.tsv')
        cmd = f'./diamond blastp -q {new_seqs_file} -d {consensusdb} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamond_result}'
        ret = run_command(cmd)
        if ret != 0:
            raise Exception('Error running diamond blastp')
        top_match={}
        for line in open(diamond_result, 'r'):
            cells = line.rstrip().split('\t')            

            sid=cells[0]
            clustername=cells[1]
            pident = float(cells[2]) / 100
            alignment_length = int(cells[3]) # * 3
            qlen = int(cells[12])# * 3 + 3
            slen = int(cells[13])# * 3 + 3
            short_seq = min(qlen, slen)
            long_seq = max(qlen, slen)
            len_diff = short_seq / long_seq
            align_short = alignment_length / short_seq
            align_long = alignment_length / long_seq
            if pident <= 0.7 or len_diff <= 0.7 or align_short <= 0.7 or align_long <= 0.7:
                continue
            if not cells[0] in top_match.keys():
                top_match[cells[0]]={'len':qlen,'cluster':clustername,'ident':pident} 
            if top_match[cells[0]]['ident']<pident:
                top_match[cells[0]]['cluster']=clustername
                top_match[cells[0]]['ident']=pident
        for seq in top_match.keys():
            old_clusters[top_match[seq]['cluster']]['gene_id'].append(seq)
            new_size=old_clusters[top_match[seq]['cluster']]['size']+1
            old_clusters[top_match[seq]['cluster']]['size']=new_size
            old_clusters[top_match[seq]['cluster']]['mean_length']=float((int(old_clusters[top_match[seq]['cluster']]['mean_length'])*(new_size-1)+int(top_match[seq]['len'])))/new_size
            old_clusters[top_match[seq]['cluster']]['max_length']=max(int(old_clusters[top_match[seq]['cluster']]['max_length']),int(top_match[seq]['len']))
            old_clusters[top_match[seq]['cluster']]['min_length']=min(int(old_clusters[top_match[seq]['cluster']]['min_length']),int(top_match[seq]['len']))
        un_match_combined_faa_file = os.path.join(out_dir, 'temp', 'unmatch_combined.faa')
        with open(un_match_combined_faa_file, 'w') as fh, open(new_seqs_file, 'rt') as fi:
            for newseq in SeqIO.parse(fi,'fasta'):
                if not newseq.id in top_match.keys():
                    SeqIO.write(newseq,fh,'fasta')
        elapsed = datetime.now() - starttime
        logging.info(f'Matching new seqs to existed clusters -- time taken {str(elapsed)}')
        return un_match_combined_faa_file
    else:
        return None
 
            

        

