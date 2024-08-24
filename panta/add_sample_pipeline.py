import os
import logging
import copy
from datetime import datetime
import multiprocessing

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

def match_new_sequence_to_oldcluster(new_seqs_file,old_clusters,out_dir,groups,group_dir,consensusdb,method='diamond',evalue=1E-6,threads=1):
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
            # if not cells[0] in top_match.keys():
            #     top_match[cells[0]]={'len':qlen,'cluster':clustername,'ident':pident} 
            # if top_match[cells[0]]['ident']<pident:
            #     top_match[cells[0]]['cluster']=clustername
            #     top_match[cells[0]]['ident']=pident
            if not sid in top_match.keys():
                top_match[sid]=[]
            top_match[sid].append({'len':qlen,'cluster':clustername,'ident':pident} )
        temseqdir=os.path.join(out_dir, 'temp_seqs')
        temmatchingdir=os.path.join(out_dir, 'temp_matching')
        if not os.path.exists(temseqdir):
            os.mkdir(temseqdir)
        if not os.path.exists(temmatchingdir):
            os.mkdir(temmatchingdir)
        with open(new_seqs_file, 'r') as fh:
            for seq in SeqIO.parse(fh,'fasta'):
                temp_seq= os.path.join(temseqdir, seq.id+".faa")
                with open(temp_seq, 'w') as fo:
                    SeqIO.write(seq,fo,'fasta')
        #blast with groups in matched clusters
        pool = multiprocessing.Pool(processes=threads)
        results = []
        for sid in top_match.keys():
            os.mkdir(temmatchingdir+"/"+sid)
            if len(top_match[sid])<=1:
                continue
            for c in top_match[sid]:
                diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
                cmd = f'./diamond blastp -q {os.path.join(temseqdir,sid+".faa")} -d {"clusters/"+c["cluster"]+"/"+c["cluster"]+".db.dmnd"} -p 1 --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamondout}'
                results.append(pool.apply_async(run_command,(cmd, None),error_callback=custom_error_callback))
        pool.close()
        pool.join()
        # for result in results:
        #     if result.get() != 0:
        #         #print(result)
        #         raise Exception('Error running diamond with ref clusters')
        for sid in top_match.keys():
            neartest_cluster=top_match[sid][0]
            best_ident=0
            if len(top_match[sid])>1:
                for c in top_match[sid]:
                    diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
                    #read and note max identity group
                    for line in open(diamondout, 'r'):
                        cells = line.rstrip().split('\t')            

                        sid=cells[0]
                        groupname=cells[1]
                        pident = float(cells[2]) / 100
                        if pident>best_ident:
                            best_ident=pident
                            neartest_cluster=c
            nearest_cluster_name=neartest_cluster['cluster']
           
            old_clusters[nearest_cluster_name]['groups'].append(groups[sid])
            new_size=old_clusters[nearest_cluster_name]['size']+1+len(groups[sid]['gene_id'])
            old_clusters[nearest_cluster_name]['size']=new_size
            old_clusters[nearest_cluster_name]['mean_length']=float((int(old_clusters[nearest_cluster_name]['mean_length'])*(new_size-1-len(groups[sid]['gene_id']))+int(neartest_cluster['len'])))/(new_size-len(groups[sid]['gene_id']))
            old_clusters[nearest_cluster_name]['max_length']=max(int(old_clusters[nearest_cluster_name]['max_length']),int(neartest_cluster['len']))
            old_clusters[nearest_cluster_name]['min_length']=min(int(old_clusters[nearest_cluster_name]['min_length']),int(neartest_cluster['len']))
            del groups[sid]
        un_match_combined_faa_file = os.path.join(out_dir,  'unmatch_combined.faa')
        count_unmatched=0
        with open(un_match_combined_faa_file, 'w') as fh, open(new_seqs_file, 'rt') as fi:
            for newseq in SeqIO.parse(fi,'fasta'):
                if not newseq.id in top_match.keys():
                    count_unmatched=count_unmatched+1
                    newseq.description=''
                    newseq.name=''
                    SeqIO.write(newseq,fh,'fasta')
        elapsed = datetime.now() - starttime
        logging.info(f'Matching {len(top_match.keys())} new groups to {len(old_clusters.keys())} existed clusters, remain {count_unmatched} groups -- time taken {str(elapsed)}')
        return un_match_combined_faa_file,groups
    else:
        return None
def custom_error_callback(error):
	print(f'Got error: {error}')
def extend_by_match_seqs_to_ref(seqs_file,groups, old_clusters, refdb, refcdb, ref_clusters,out_dir, threads=1, evalue=1E-6):
    starttime = datetime.now()
    diamond_result = os.path.join(out_dir, 'extend_ref_diamond.tsv')
    cmd = f'./diamond blastp -q {seqs_file} -d {refdb} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamond_result}'
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
        if not sid in top_match.keys():
            top_match[sid]=[]
        top_match[sid].append({'len':qlen,'cluster':clustername,'ident':pident} )
    un_match_combined_faa_file = os.path.join(out_dir, 'unmatch_combined2.faa')
    
    temseqdir=os.path.join(out_dir, 'temp_seqs')
    temmatchingdir=os.path.join(out_dir, 'temp_matching')
    if not os.path.exists(temseqdir):
        os.mkdir(temseqdir)
    if not os.path.exists(temmatchingdir):
        os.mkdir(temmatchingdir)
    with open(seqs_file, 'r') as fh:
        for seq in SeqIO.parse(fh,'fasta'):
            temp_seq= os.path.join(temseqdir, seq.id+".faa")
            if os.path.exists(temp_seq):
                continue
            with open(temp_seq, 'w') as fo:
                SeqIO.write(seq,fo,'fasta')
        #blast with groups in matched clusters
    pool = multiprocessing.Pool(processes=threads)
    results = []
    for sid in top_match.keys():
        if not os.path.exists(temmatchingdir+"/"+sid):
            os.mkdir(temmatchingdir+"/"+sid)
        if len(top_match[sid])<=1:
            continue
        for c in top_match[sid]:
            diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
            cmd = f'./diamond blastp --quiet  -q {os.path.join(temseqdir,sid+".faa")} -d {refcdb+"/"+c["cluster"]+".db.dmnd"} -p 1 --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamondout}'
            results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running diamond with ref clusters')
    gene_map = {}
    count = 0
    
    
    for sid in top_match.keys():
        neartest_cluster=top_match[sid][0]
        best_ident=0
        if len(top_match[sid])>1:
            for c in top_match[sid]:
                diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
                #read and note max identity group
                for line in open(diamondout, 'r'):
                    cells = line.rstrip().split('\t')            

                    #sid=cells[0]
                    groupname=cells[1]
                    pident = float(cells[2]) / 100
                    if pident>best_ident:
                        best_ident=pident
                        neartest_cluster=c
        nearest_cluster_name=neartest_cluster['cluster']
        if not nearest_cluster_name in old_clusters.keys():
            logging.info(f'New cluster to extends:{nearest_cluster_name}')
            old_clusters[nearest_cluster_name]={}
            old_clusters[nearest_cluster_name]['groups']=[]
            old_clusters[nearest_cluster_name]['max_length']=0
            old_clusters[nearest_cluster_name]['mean_length']=0
            old_clusters[nearest_cluster_name]['min_length']=1E6
            old_clusters[nearest_cluster_name]['gene_name']=ref_clusters[nearest_cluster_name]['gene_name']
            old_clusters[nearest_cluster_name]['product']=ref_clusters[nearest_cluster_name]['description']
            old_clusters[nearest_cluster_name]['representative']=''
            old_clusters[nearest_cluster_name]['source']='reference'
            old_clusters[nearest_cluster_name]['size']=0
        old_clusters[nearest_cluster_name]['groups'].append(groups[sid])
        #old_clusters[nearest_cluster_name]['gene_id'].extend(groups[sid])
        #TODO: need to recalculate 
        new_size=old_clusters[nearest_cluster_name]['size']+1+len(groups[sid]['gene_id'])
        old_clusters[nearest_cluster_name]['size']=new_size
        old_clusters[nearest_cluster_name]['mean_length']=float((int(old_clusters[nearest_cluster_name]['mean_length'])*(new_size-len(groups[sid]['gene_id']))+int(neartest_cluster['len'])))/(new_size-len(groups[sid]['gene_id']))
        old_clusters[nearest_cluster_name]['max_length']=max(int(old_clusters[nearest_cluster_name]['max_length']),int(neartest_cluster['len']))
        old_clusters[nearest_cluster_name]['min_length']=min(int(old_clusters[nearest_cluster_name]['min_length']),int(neartest_cluster['len']))
        del groups[sid]
    #json.dump(clusters, open(os.path.join(out_dir, 'ref_by_clusters.json'), 'w'), indent=4, sort_keys=True)
    with open(un_match_combined_faa_file, 'w') as fh, open(seqs_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            if not newseq.id in top_match.keys():
                SeqIO.write(newseq,fh,'fasta')
    elapsed = datetime.now() - starttime
    logging.info(f'Matched {len(top_match.keys())} new groups by ref to form {len(old_clusters.keys())}  clusters -- time taken {str(elapsed)}')
    return un_match_combined_faa_file,old_clusters,groups
    
            

        

