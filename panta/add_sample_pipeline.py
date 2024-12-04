import os
import logging
import copy
from datetime import datetime
import multiprocessing

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from panta.utils import run_command, parse_cluster_file,getIdentAlignFromCell

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

def match_new_sequence_to_oldcluster(new_seqs_file,old_clusters,out_dir,groups,consensusdb,method='diamond',evalue=1E-6,threads=1):
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
            pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
                
            if pident <= 0.7 or len_diff <= 0.7 or align_short <= 0.7 or align_long <= 0.7:
                continue
            if not cells[0] in top_match.keys():
                top_match[cells[0]]={'len':qlen,'cluster':clustername,'ident':pident} 
            if top_match[cells[0]]['ident']<pident:
                top_match[cells[0]]['cluster']=clustername
                top_match[cells[0]]['ident']=pident
            # if not sid in top_match.keys():
            #     top_match[sid]=[]
            # top_match[sid].append({'len':qlen,'cluster':clustername,'ident':pident} )
        # temseqdir=os.path.join(out_dir, 'temp_seqs')
        # temmatchingdir=os.path.join(out_dir, 'temp_matching')
        # if not os.path.exists(temseqdir):
        #     os.mkdir(temseqdir)
        # if not os.path.exists(temmatchingdir):
        #     os.mkdir(temmatchingdir)
        # with open(new_seqs_file, 'r') as fh:
        #     for seq in SeqIO.parse(fh,'fasta'):
        #         temp_seq= os.path.join(temseqdir, seq.id+".faa")
        #         with open(temp_seq, 'w') as fo:
        #             SeqIO.write(seq,fo,'fasta')
        # #blast with groups in matched clusters
        
        # pool = multiprocessing.Pool(processes=threads)
        # results = []
        # matching_file=open(os.path.join(out_dir,"matching_old_cluster_log.txt"),'w')
        # for sid in top_match.keys():
        #     if not os.path.exists(temmatchingdir+"/"+sid):
            
        #         os.mkdir(temmatchingdir+"/"+sid)
        #     #if len(top_match[sid])<=1:
        #     #    continue
            
        #     for c in top_match[sid]:
        #         diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
        #         cmd = f'./diamond blastp -q {os.path.join(temseqdir,sid+".faa")} -d {os.path.join(cluster_dir,c["cluster"]+"/"+c["cluster"]+".db.dmnd")} -p 1 --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamondout}'
        #         results.append(pool.apply_async(run_command,(cmd, None),error_callback=custom_error_callback))
        # pool.close()
        # pool.join()    
        # for result in results:
        #     if result.get() != 0:
        #         #print(result)
        #         raise Exception('Error running diamond with ref clusters')
            
        added_group=set()
        for sid in top_match.keys():
            neartest_cluster=top_match[sid]
            # best_ident=0
            # #groups[sid]['confident']=1
            # list_matched_groups=[]
            # hash_clustername={}
            # hash_cluster_max_ident={}
            # for c in top_match[sid]:
            #     #print({'c':c})
            #     diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
            #     #print(diamondout)
            #     #read and note max identity group
            #     hash_cluster_max_ident[c['cluster']]={'pident':0,'len_diff':0,'align_short':0,'align_long':0}
            #     for line in open(diamondout, 'r'):
            #         cells = line.rstrip().split('\t')            

            #         sid=cells[0]
            #         groupname=cells[1]
            #         #pident = float(cells[2]) / 100
            #         pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
            #         matching_file.write(sid+'\t'+c['cluster']+'\t'+groupname+'\t'+str(pident)+'\t'+str(len_diff)+'\n')
            #         if pident < 0.98 or len_diff < 0.98 or align_short < 0.98 or align_long < 0.98:
            #             continue
            #         list_matched_groups.append({'c':c,'g':groupname,'pid':pident})
            #         # print({'c':c,'g':groupname,'pid':pident})
            #         if pident>best_ident:
            #             best_ident=pident
            #             neartest_cluster=c
            #         if pident>hash_cluster_max_ident[c['cluster']]['pident']:
            #             hash_cluster_max_ident[c['cluster']]={'match_group':groupname,'pident':pident,'len_diff':len_diff,'align_short':align_short,'align_long':align_long}
            #     hash_clustername[c['cluster']]=c
            # #print(sid +" has "+ str(len(top_match[sid]))+" matched cluster with best identity is "+str(best_ident))
            # if best_ident<0.98:
            #     #give up, not enough envident to add to any clusters
            #     #del top_match[sid]
            #     continue
            # else:
            #     print("consider the undetermined :")
            #     list_matched_groups.sort(key=lambda x: x['pid'], reverse=True)
            #     print(list_matched_groups)
            #     count_matched_cluster={}
            #     #k=5
            #     k=5
            #     if len(list_matched_groups)<5:
            #         k=len(list_matched_groups)
            #     for i in range(k):
            #         if not list_matched_groups[i]['c']['cluster'] in count_matched_cluster:
            #             count_matched_cluster[list_matched_groups[i]['c']['cluster']]=0
            #         count_matched_cluster[list_matched_groups[i]['c']['cluster']]=count_matched_cluster[list_matched_groups[i]['c']['cluster']]+1
            #     max_in_k_nearest=0
            #     cluster_max=list_matched_groups[0]['c']['cluster']
            #     for kc in count_matched_cluster.keys():
            #         if count_matched_cluster[kc]>max_in_k_nearest:
            #             max_in_k_nearest=count_matched_cluster[kc]
            #             cluster_max=kc
            #     neartest_cluster=hash_clustername[cluster_max]
            #     groups[sid]['confident']=hash_cluster_max_ident[cluster_max]['pident']
            #     groups[sid]['match']=hash_cluster_max_ident[cluster_max]
            # else:
            #     groups[sid]['confident']=1
            #     groups[sid]['match']=hash_cluster_max_ident[neartest_cluster['cluster']]
                #print('picked cluster is '+neartest_cluster['cluster']+ ' with '+str(count_matched_cluster[cluster_max])+ ' nearest groups' )
            nearest_cluster_name=neartest_cluster['cluster']
            old_clusters[nearest_cluster_name]['unique_seq'].append(sid)
            old_clusters[nearest_cluster_name]['gene_id'].append(sid)
            old_clusters[nearest_cluster_name]['gene_id'].extend(groups[sid])
            new_size=old_clusters[nearest_cluster_name]['size']+1+len(groups[sid])
            old_clusters[nearest_cluster_name]['size']=new_size
            old_clusters[nearest_cluster_name]['mean_length']=float((int(old_clusters[nearest_cluster_name]['mean_length'])*(new_size-1-len(groups[sid]))+int(neartest_cluster['len'])))/(new_size-len(groups[sid]))
            old_clusters[nearest_cluster_name]['max_length']=max(int(old_clusters[nearest_cluster_name]['max_length']),int(neartest_cluster['len']))
            old_clusters[nearest_cluster_name]['min_length']=min(int(old_clusters[nearest_cluster_name]['min_length']),int(neartest_cluster['len']))
            del groups[sid]
            added_group.add(sid)
        #matching_file.close()
        un_match_combined_faa_file = os.path.join(out_dir,  'unmatch_combined.faa')
        count_unmatched=0
        with open(un_match_combined_faa_file, 'w') as fh, open(new_seqs_file, 'rt') as fi:
            for newseq in SeqIO.parse(fi,'fasta'):
                if not newseq.id in added_group:
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
    
    logging.info(f'Before extend ref clusters, there are {len(old_clusters.keys())} old clusters')
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
        pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
        
        if pident <= 0.7 or len_diff <= 0.7 or align_short <= 0.7 or align_long <= 0.7:
            continue
        
        if not cells[0] in top_match.keys():
            top_match[cells[0]]={'len':qlen,'cluster':clustername,'ident':pident} 
        if top_match[cells[0]]['ident']<pident:
            top_match[cells[0]]['cluster']=clustername
            top_match[cells[0]]['ident']=pident
    un_match_combined_faa_file = os.path.join(out_dir, 'unmatch_combined2.faa')
    matched_cluster=set()
    for sid in top_match.keys():
        c=top_match[sid]['cluster']
        if not c in old_clusters.keys():
            
            old_clusters[c]={}
            old_clusters[c]['gene_id']=[]
            old_clusters[c]['max_length']=0
            old_clusters[c]['mean_length']=0
            old_clusters[c]['min_length']=1E6
            old_clusters[c]['gene_name']=ref_clusters[c]['gene_name']
            old_clusters[c]['product']=ref_clusters[c]['description']
            old_clusters[c]['representative']=''
            old_clusters[c]['source']='reference'
            old_clusters[c]['size']=0
            matched_cluster.add(c)
        old_clusters[c]['gene_id'].append(sid)
        old_clusters[c]['gene_id'].extend(groups[sid])
        new_size=old_clusters[c]['size']+1+len(groups[sid])
        old_clusters[c]['size']=new_size
        old_clusters[c]['mean_length']=float((int(old_clusters[c]['mean_length'])*(new_size-1-len(groups[sid]))+int(top_match[sid]['len'])))/(new_size-len(groups[sid]))
        old_clusters[c]['max_length']=max(int(old_clusters[c]['max_length']),int(top_match[sid]['len']))
        old_clusters[c]['min_length']=min(int(old_clusters[c]['min_length']),int(top_match[sid]['len']))
        del groups[sid]
    #json.dump(clusters, open(os.path.join(out_dir, 'ref_by_clusters.json'), 'w'), indent=4, sort_keys=True)
    count_unmatched=0
    count_num_input=0
    with open(un_match_combined_faa_file, 'w') as fh, open(seqs_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            count_num_input=count_num_input+1
            if not newseq.id in top_match.keys():
                SeqIO.write(newseq,fh,'fasta')
                count_unmatched=count_unmatched+1
    elapsed = datetime.now() - starttime
    logging.info(f'Matching {count_num_input} groups to ref clusters, {len(top_match.keys())} matched, form new {len(matched_cluster)} ref clusters,  {count_unmatched} groups not matched  -- time taken {str(elapsed)}')
    return un_match_combined_faa_file,old_clusters,groups
def realign_group_to_clusters(annotated_clusters):
    groups_need_check=[]
    for cluster_name in annotated_clusters.keys():
        for group in annotated_clusters[cluster_name]['groups']:
            if group['confident']<1:
                groups_need_check.append(group)
    #create seq combined all group to check
def group_new_unique_seq_to_old_cluster(old_unique_seqs_fasta,new_unique_seqs_fasta, old_clusters,map_unique_seq_cluster,new_unique_groups,out_dir ):
    #concat
    new_merged_seq_fasta=os.path.join(out_dir,"unique_concat_seq.fasta")
    cmd=f'cat {old_unique_seqs_fasta} {new_unique_seqs_fasta} > {new_merged_seq_fasta} '
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error concat unique sequences')
    mmseq_cluster_file=os.path.join(out_dir, 'mmseq_unique_seq')
    cmd = f'mmseqs easy-linclust {new_merged_seq_fasta} {mmseq_cluster_file} {out_dir}/tmp --min-seq-id 1 -c 1  --threads {threads} > /dev/null'    
    ret = run_command(cmd)
    c_cursor=0
    cluster_count=0
    unique_groups={}
    matched_seq=set()
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=rep_name                
                unique_groups[c_cursor] = set()
                unique_groups[c_cursor].add(rep_name)
                
                
            else:
                unique_groups[c_cursor].add(member)
    for uid in map_unique_seq_cluster:
        isFound=False
        for uid2 in unique_groups:
            if uid in unique_groups[uid2]:
                
                
                point_cluster=map_unique_seq_cluster[uid]
                set1=set(old_clusters[point_cluster]['unique_seq'])
                set2=unique_groups[uid2]
                set3=set1.union(set2)
                matched_seq.update(set2)
                old_clusters[point_cluster]['unique_seq']=list(set3)
                
                set_g1=set(old_clusters[point_cluster]['gene_id'])
                set_g2=set(new_clusters[max_matched_cluster]['gene_id'])
                old_clusters[point_cluster]['gene_id']=list(set_g1.union(set_g2))
                
                set3=set1.union(set2)
                break
    # remain sequence:
    un_match_unique_seqs_file=os.path.join(out_dir,"unmatched_unique_seqs.fasta")
    with open(un_match_unique_seqs_file, 'w') as fh, open(new_unique_seqs_fasta, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            if not newseq.id in matched_seq:
                SeqIO.write(newseq,fh,'fasta')    
    new_merged_unique_seq_fasta=os.path.join(out_dir,"new_unique_merged_seq.fasta")
    cmd=f'cat {old_unique_seqs_fasta} {un_match_unique_seqs_file} > {new_merged_unique_seq_fasta} '
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error concat unique sequences')
    return new_merged_unique_seq_fasta,old_clusters
        
def merge_clusters_based_on_unique_seq(old_clusters,new_clusters):
    for c1 in old_clusters:
        set1=set(old_clusters[c1]['unique_seq'])
        max_matched_cluster=None
        max_matched_number_unique_seq=1
        for c2 in new_clusters:           
            set2=set(new_clusters[c2]['unique_seq'])
            n=set2.intersection(set1)
            if len(n)>max_matched_number_unique_seq:
                max_matched_cluster=c2
                max_matched_number_unique_seq=len(n)
        if not max_matched_cluster is None and float(max_matched_number_unique_seq)/len(set1)>0.5:
            
            set2=set(new_clusters[max_matched_cluster]['unique_seq'])
            old_clusters[c1]['unique_seq']=list(set1.union(set2))
            set_g1=set(old_clusters[c1]['gene_id'])
            set_g2=set(new_clusters[max_matched_cluster]['gene_id'])
            old_clusters[c1]['gene_id']=list(set_g1.union(set_g2))
            del new_clusters[max_matched_cluster]
    return old_clusters,new_clusters
        
def reinflate_clusters_with_hirachical_group(old_unique_groups,new_unique_groups, combined_unique_group,similar_groups, mcl_file):
    """
    Return
    ------
        - inflated_clusters: list of list of genes
        -clusters: dict(cluster_id->[gene_ids])
    """    
    starttime = datetime.now()
    clusters = {}
    clusters.update(similar_groups)

    inflated_clusters = []
    unique_seq_by_similar_group={}
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = {}
            line = line.rstrip('\n')
            genes = line.split('\t')
            
            for s_gene in genes:
                #inflated_genes.append(gene)
                if s_gene in similar_groups:
                    #inflated_genes.extend(groups[gene]['gene_id'])

                    inflated_genes[s_gene]=[]
                    gene_id_set=set()
                    
                    gene_id_set.add(s_gene)
                    #inflated_genes[gene].extend(unique_groups[gene])
                    unique_seq_by_similar_group[s_gene]=[]
                    if s_gene in old_unique_groups:
                        gene_id_set.update(old_unique_groups[s_gene])
                    if s_gene in new_unique_groups:
                        gene_id_set.update(new_unique_groups[s_gene])
                    for u_seq in combined_unique_group[s_gene]:
                        gene_id_set.add(u_seq)
                        if u_seq in old_unique_groups:
                            gene_id_set.update(old_unique_groups[u_seq])
                        if u_seq in new_unique_groups:
                            gene_id_set.update(new_unique_groups[u_seq])
                    for cu_gene in similar_groups[s_gene]:
                        gene_id_set.add(cu_gene)
                        unique_seq_by_similar_group[s_gene].append(cu_gene)
                        if cu_gene in old_unique_groups:
                            gene_id_set.update(old_unique_groups[cu_gene])
                        if cu_gene in new_unique_groups:
                            gene_id_set.update(new_unique_groups[cu_gene])
                        for u_seq in combined_unique_group[cu_gene]:
                            gene_id_set.add(u_seq)
                            if u_seq in old_unique_groups:
                                gene_id_set.update(old_unique_groups[u_seq])
                            if u_seq in new_unique_groups:
                                gene_id_set.update(new_unique_groups[u_seq])
                        #inflated_genes[gene].extend(unique_groups[g])
                        
                    inflated_genes[s_gene]=list(gene_id_set)
                    del similar_groups[s_gene]
            inflated_clusters.append(inflated_genes)
    
    # Inflate any clusters that were in the clusters file but not sent to mcl
    count_not_mcl=0
    
    for s_gene in similar_groups:
        count_not_mcl=count_not_mcl+1
        #logging.info(f'Not in mcl {s_gene}')
        inflated_genes={}
        inflated_genes[s_gene]=[]
        gene_id_set=set()
        gene_id_set.add(s_gene)
        unique_seq_by_similar_group[s_gene]=[]
        if s_gene in old_unique_groups:
            gene_id_set.update(old_unique_groups[s_gene])
        if s_gene in new_unique_groups:
            gene_id_set.update(new_unique_groups[s_gene])
        for u_seq in combined_unique_group[s_gene]:
            gene_id_set.add(u_seq)
            if u_seq in old_unique_groups:
                gene_id_set.update(old_unique_groups[u_seq])
            if u_seq in new_unique_groups:
                gene_id_set.update(new_unique_groups[u_seq])
        for cu_gene in similar_groups[s_gene]:
            gene_id_set.add(cu_gene)
            unique_seq_by_similar_group[s_gene].append(cu_gene)
            if cu_gene in old_unique_groups:
                gene_id_set.update(old_unique_groups[cu_gene])
            if cu_gene in new_unique_groups:
                gene_id_set.update(new_unique_groups[cu_gene])
            for u_seq in combined_unique_group[cu_gene]:
                gene_id_set.add(u_seq)
                if u_seq in old_unique_groups:
                    gene_id_set.update(old_unique_groups[u_seq])
                if u_seq in new_unique_groups:
                    gene_id_set.update(new_unique_groups[u_seq])
                        #inflated_genes[gene].extend(unique_groups[g])                       
        inflated_genes[s_gene]=list(gene_id_set)
        #logging.info(inflated_genes)         
        inflated_clusters.append(inflated_genes)
    
    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate new {len(inflated_clusters)} clusters with refer to {len(clusters.keys())} groups, {count_not_mcl} groups not found in MCL clustering -- time taken {str(elapsed)}')
    return inflated_clusters, clusters,unique_seq_by_similar_group
def flat_combined_unique_seq(new_unique_groups, old_unique_groups, combined_unique_group):
    flat_unique_seq={}
    for cu_gene in combined_unique_group:
        flat_unique_seq[cu_gene]=[]
        if cu_gene in old_unique_groups:
            flat_unique_seq[cu_gene].extend(old_unique_groups[cu_gene])
        if cu_gene in new_unique_groups:
            flat_unique_seq[cu_gene].extend(new_unique_groups[cu_gene])
        for u_seq in combined_unique_group[cu_gene]:
            flat_unique_seq[cu_gene].append(u_seq)
            if u_seq in old_unique_groups:
                flat_unique_seq[cu_gene].extend(old_unique_groups[u_seq])
            if u_seq in new_unique_groups:
                flat_unique_seq[cu_gene].extend(new_unique_groups[u_seq])
    return flat_unique_seq