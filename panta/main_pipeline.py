import os
import logging
import multiprocessing
from datetime import datetime
from panta.utils import run_command, parse_cluster_file, chunk_fasta_file,getIdentAlignFromCell
from panta.post_analysis import annotate_cluster
import json
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import gc
from panta.utils import *
import faiss
import numpy as np
import torch
import esm
import pandas as pd
from collections import Counter
logger = logging.getLogger(__name__)


def run_cd_hit(faa_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    cd_hit_represent_fasta = os.path.join(out_dir, 'cd-hit.fasta')
    cd_hit_cluster_file = cd_hit_represent_fasta + '.clstr'
    cmd = f'cd-hit -i {faa_file} -o {cd_hit_represent_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running cd-hit')
    cd_hit_clusters = parse_cluster_file(cd_hit_cluster_file)

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity -- time taken {str(elapsed)}')
    return cd_hit_represent_fasta, cd_hit_clusters


def run_cd_hit_with_map(faa_file, map_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    cd_hit_represent_fasta = os.path.join(out_dir, 'cd-hit_tmp.fasta')
    cd_hit_cluster_file = cd_hit_represent_fasta + '.clstr'
    cmd = f'cd-hit -i {faa_file} -o {cd_hit_represent_fasta} -s 0.98 -c 0.98 -T {threads} -M 0 -g 1 -d 256 > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running cd-hit')        

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity part 1 -- time taken {elapsed}')

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    cd_hit_represent_corrected_fasta = os.path.join(out_dir, 'cd-hit.fasta')
    with open(cd_hit_represent_corrected_fasta,'w') as ofh, open(cd_hit_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                ofh.write(f'>{gene_map[line[1:].strip()]}\n')
            else:
                ofh.write(line)      

    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity part 1 -- time taken {elapsed}')
    with open(cd_hit_cluster_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line[0].startswith('>'):
                cluster_name = line[1:]
                clusters[cluster_name] = {'gene_names':[]}                
            else:
                _,_, line = line.partition(', >')            
                gene_name,_,identity = line.partition('... ')
                gene_name, identity                    
                if identity == '*':
                    clusters[cluster_name]['representative'] = gene_map[gene_name]
                elif identity: # make sure it is a valid string
                    clusters[cluster_name]['gene_names'].append(gene_map[gene_name])                
    
    del gene_map    

    # convert to a simple dictionary
    clusters_new = {}
    for cluster_name in clusters:
        clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
    elapsed = datetime.now() - starttime
    logging.info(f'Run CD-HIT with 98% identity -- time taken {str(elapsed)}')
    return cd_hit_represent_corrected_fasta, clusters_new

def run_mmseq_with_map(faa_file, map_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(out_dir, 'mmseq_cluster.tsv')
    cmd = f'mmseqs easy-linclust {faa_file} {out_dir}/mmseq {out_dir}/tmp --min-seq-id 0.98  -c 0.98 --threads {threads} > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseq')        

    elapsed = datetime.now() - starttime
    #logging.info(f'Run mmseq2 with 98% identity part 1 -- time taken {elapsed}')

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    represent_corrected_fasta = os.path.join(out_dir, 'mmseq.fasta')
    with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                ofh.write(f'>{gene_map[line[1:].strip()]}\n')
            else:
                ofh.write(line)    
    #with open(represent_corrected_fasta) as handle:
    #    for record in SeqIO.parse(handle, "fasta"):
    #        SeqIO.write(record, os.path.join(group_dir,record.id+".faa"), "fasta")            

    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 98% identity part 1 -- time taken {elapsed}')
    c_cursor=None
    cluster_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=gene_map[member]
                clusters[gene_map[member]]=[]
                #clusters[gene_map[member]] = {'gene_id':[gene_map[member]]} 
                #clusters[gene_map[member]]['representative'] = gene_map[member]
                cluster_count=cluster_count+1
                
            else:
                clusters[c_cursor].append(gene_map[member])                

           
                
    
    del gene_map    

    # convert to a simple dictionary
    # clusters_new = {}
    # for cluster_name in clusters:
    #     clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 98% identity, {cluster_count} groups -- time taken {str(elapsed)}')
    return represent_corrected_fasta, clusters
def run_mmseq_with_map_unique_seqs(faa_file, map_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(out_dir, 'mmseq_cluster.tsv')
    cmd = f'mmseqs easy-cluster {faa_file} {out_dir}/mmseq {out_dir}/tmp --min-seq-id 1 -c 1 --cov-mode 0 --threads {threads} > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseq')        

    elapsed = datetime.now() - starttime
    #logging.info(f'Run mmseq2 with 98% identity part 1 -- time taken {elapsed}')

    groups = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    represent_corrected_fasta = os.path.join(out_dir, 'unique_seqs.fasta')
    with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                ofh.write(f'>{gene_map[line[1:].strip()]}\n')
            else:
                ofh.write(line)    
    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 100% identity toget uniqued seqs -- time taken {elapsed}')
    c_cursor=None
    cluster_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=gene_map[member]
                groups[gene_map[member]]=[]
                #clusters[gene_map[member]] = {'gene_id':[gene_map[member]]} 
                #clusters[gene_map[member]]['representative'] = gene_map[member]
                cluster_count=cluster_count+1
                
            else:
                groups[c_cursor].append(gene_map[member])                

           
                
    
    del gene_map    

    # convert to a simple dictionary
    # clusters_new = {}
    # for cluster_name in clusters:
    #     clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 100% identity, {cluster_count} groups -- time taken {str(elapsed)}')
    return represent_corrected_fasta, groups
def run_mmseq_unique_seqs(faa_file, out_dir, threads=4, timing_log=None):        
    starttime = datetime.now()
    mem_usage = mem_report(0, "begin run_mmseq_unique_seqs ")
    mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(out_dir, 'mmseq_cluster.tsv')
    cmd = f'mmseqs easy-cluster {faa_file} {out_dir}/mmseq {out_dir}/tmp --min-seq-id 1 -c 1 --cov-mode 0 --threads {threads} > /dev/null'    
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mmseq')        
    mem_usage = mem_report(mem_usage, cmd)
    elapsed = datetime.now() - starttime
    #logging.info(f'Run mmseq2 with 98% identity part 1 -- time taken {elapsed}')

    groups = {}
    
    represent_corrected_fasta = os.path.join(out_dir, 'unique_seqs.fasta')
    with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
        for line in ifh:
            ofh.write(line)    
    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 100% identity toget uniqued seqs -- time taken {elapsed}')
    mem_usage = mem_report(mem_usage, "create represent_corrected_fasta")
    c_cursor=None
    cluster_count=0
    seq_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=member
                groups[member]=[]
                cluster_count=cluster_count+1
                
            else:
                groups[c_cursor].append(member)                

            seq_count=seq_count+1   
                
    
       

    mem_usage = mem_report(mem_usage, "create groups")
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 100% identity, {cluster_count} groups from {seq_count} -- time taken {str(elapsed)}')
    return represent_corrected_fasta, groups
import mmh3
def run_hash_unique_seqs(faa_file, out_dir, threads=4, timing_log=None):        
    starttime = datetime.now()
    #mem_usage = mem_report(0, "begin run_hash_unique_seqs ")
    hash_represent_fasta= os.path.join(out_dir, 'hash_rep_seq.fasta')
    gene_hash={}
    represent_corrected_fasta = os.path.join(out_dir, 'unique_seqs.fasta')
    seq_count=0
    map_rep_hash={}
    with open(faa_file) as faa,open(represent_corrected_fasta,"w") as newfaa:
        for record  in SeqIO.parse(faa,"fasta"):
            hash=str(mmh3.hash128(str(record.seq)))
            seq_count=seq_count+1
            if hash in gene_hash:
                gene_hash[hash].append(record.id)
            else:
                gene_hash[hash]=[record.id]
                map_rep_hash[record.id]=hash
                SeqIO.write(record,newfaa,'fasta')
    
    elapsed = datetime.now() - starttime
    logging.info(f'Run hash , {len(gene_hash)} hash from {seq_count} -- time taken {str(elapsed)}')
    return represent_corrected_fasta, gene_hash,map_rep_hash
def add_hash_unique_seqs(old_gene_hash,faa_file, out_dir, threads=4, timing_log=None):        
    starttime = datetime.now()
    logging.info(f'Start add hash')
    
    
    #gene_hash={}
    remain_new_faa = os.path.join(out_dir, 'remain_seqs.faa')
    seq_count=0
    unmatch_seq_count=0
    match_seq_count=0
    with open(faa_file) as faa,open(remain_new_faa,"w") as newfaa:
        for record  in SeqIO.parse(faa,"fasta"):
            hash=str(mmh3.hash128(str(record.seq)))
            seq_count=seq_count+1
            if hash in old_gene_hash:
                match_seq_count=match_seq_count+1
                old_gene_hash[hash].append(record.id)
            else:
                unmatch_seq_count=unmatch_seq_count+1
                SeqIO.write(record,newfaa,'fasta')
  
    elapsed = datetime.now() - starttime
    logging.info(f'Finishs add hash , add {match_seq_count},  size hash is {len(old_gene_hash)} hash, remain {unmatch_seq_count} seqs -- time taken {str(elapsed)}')
    return  old_gene_hash,remain_new_faa

def convertSeq2KmerCount(k,index_kaa,seq):
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    kmer_counts = Counter(kmers)
    kc_profile=[0]*len(index_kaa)
    for kmer,count in kmer_counts.items():
        kc_profile[index_kaa[kmer]]=count
    return kc_profile
def protein_sequences_to_2mer(fasta_file):
    aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','S','P','F','T','W','Y','V','X']
    k2aa=[]
    for a in aa:
        for b in aa:
            k2aa.append(a+b)
        index_kaa={}
    for i in range(len(k2aa)):
        index_kaa[k2aa[i]]=i
    index_seq_id=[]
    list_vectors=[]
   
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):          
            index_seq_id.append(record.id)          
            list_vectors.append(convertSeq2KmerCount(2,index_kaa,str(record.seq)))
    return index_seq_id,list_vectors
def protein_sequences_to_vector_ems(fasta_file):
    starttime = datetime.now()
    index_seq_id=[]
    list_vectors=[]
    data=[]
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):          
            index_seq_id.append(record.id)          
            data.append((record.id,str(record.seq)))
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    elapsed = datetime.now() - starttime
    
    for i in range(len(batch_tokens)):
        list_vectors.append(batch_tokens[i].numpy().tolist())
    logging.info(f'protein to {len(list_vectors)} vector by ems -- time taken {elapsed}')
    return batch_labels,list_vectors
def transform_to_l1_space(data):
    
    return np.hstack([data, -data])

def run_faiss_unique_seqs(faa_file, out_dir, threads=4, timing_log=None):        
    starttime = datetime.now()
    mem_usage = mem_report(0, "begin run_faiss_unique_seqs ")
    #load sequence to faiss       
    index_seq,vectors=protein_sequences_to_2mer(faa_file)
    array_input = np.array(vectors)
    transformed_data = transform_to_l1_space(array_input)
    index = faiss.IndexFlatL2(transformed_data.shape[1])  # Sử dụng FAISS với L2 Distance
    index.add(transformed_data)  # Thêm dữ liệu vào chỉ số
    mem_usage = mem_report(mem_usage, "faiss input")
    elapsed = datetime.now() - starttime
    logging.info(f'load faiss index -- time taken {elapsed}')
    #search
    distances, indices = index.search(transformed_data, k=transformed_data.shape[0])
    elapsed = datetime.now() - starttime
    logging.info(f'search faiss -- time taken {elapsed}')
    identical_groups = {}
    #print(distances)
    visited = set()
    for i in range(len(transformed_data)):
        if i in visited:
            continue
        identical_groups[i]=[]
        # Find vectors with zero distance (or within the tolerance)
        identical = np.where(distances[i] <= 0)[0]
        #print(identical)
        # Add to the group if there are duplicates
        if len(identical) > 0:
            for j in identical:
                identical_groups[i].append(indices[i][j])
                visited.add(indices[i][j])
    #print(identical_groups)
    groups = {}
    for i in identical_groups.keys():
        groups[index_seq[i]]=[]
        for j in identical_groups[i]:
             groups[index_seq[i]].append(index_seq[j])
    elapsed = datetime.now() - starttime
    logging.info(f'handle faiss result -- time taken {elapsed}')
    represent_corrected_fasta = os.path.join(out_dir, 'unique_seqs.fasta')
    
    
    mem_usage = mem_report(mem_usage, "handle faiss result")
      
                
    with open(faa_file) as handle, open(represent_corrected_fasta,'w') as fo:
        for record in SeqIO.parse(handle, "fasta"): 
            if record.id in groups.keys():
                SeqIO.write(record,fo,'fasta')
       

    
    elapsed = datetime.now() - starttime
    logging.info(f'Run faiss with 100% identity, {len(groups.keys())} groups from {len(index_seq)} sequences-- time taken {str(elapsed)}')
    return represent_corrected_fasta, groups
def run_mmseq_with_map_similar_seqs(faa_file, out_dir, threads=4,timing_log=None, identity=0.98):        
    starttime = datetime.now()
    
    mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(out_dir, 'mmseq_cluster.tsv')
    cmd = f'mmseqs easy-linclust {faa_file} {out_dir}/mmseq {out_dir}/tmp --min-seq-id {identity} -c {identity} --threads {threads} > /dev/null'    
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mmseq')        

    elapsed = datetime.now() - starttime
    #logging.info(f'Run mmseq2 with 98% identity part 1 -- time taken {elapsed}')
    represent_corrected_fasta = os.path.join(out_dir, 'similar_seqs.fasta')
    with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                ofh.write(f'>{line.split(" ",1)[0][1:].strip()}\n')
            else:
                ofh.write(line)    
    
    groups = {}
        
    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 98% identity toget similar seqs -- time taken {elapsed}')
    c_cursor=None
    cluster_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=rep_name
                groups[member]=[]
                #clusters[gene_map[member]] = {'gene_id':[gene_map[member]]} 
                #clusters[gene_map[member]]['representative'] = gene_map[member]
                cluster_count=cluster_count+1
                
            else:
                groups[c_cursor].append(member)                

           
                
    
    #del gene_map    

    # convert to a simple dictionary
    # clusters_new = {}
    # for cluster_name in clusters:
    #     clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 98% identity, {cluster_count} groups -- time taken {str(elapsed)}')
    return represent_corrected_fasta, groups
def run_mmseq_unique_seqs(faa_file, out_dir, threads=4,timing_log=None):        
    starttime = datetime.now()
    
    mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(out_dir, 'mmseq_cluster.tsv')
    cmd = f'mmseqs easy-cluster {faa_file} {out_dir}/mmseq {out_dir}/tmp --min-seq-id 1 -c 1 --threads {threads} > /dev/null'    
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mmseq')        

    elapsed = datetime.now() - starttime
    #logging.info(f'Run mmseq2 with 98% identity part 1 -- time taken {elapsed}')
    represent_corrected_fasta = os.path.join(out_dir, 'unique_seqs.fasta')
    with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
        for line in ifh:
            ofh.write(line)    
    
    groups = {}
        
    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 100% identity to get unique seqs -- time taken {elapsed}')
    c_cursor=None
    cluster_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=rep_name
                groups[member]=[]
                #clusters[gene_map[member]] = {'gene_id':[gene_map[member]]} 
                #clusters[gene_map[member]]['representative'] = gene_map[member]
                cluster_count=cluster_count+1
                
            else:
                groups[c_cursor].append(member)                

           
                
    
    #del gene_map    

    # convert to a simple dictionary
    # clusters_new = {}
    # for cluster_name in clusters:
    #     clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 100% identity, {cluster_count} groups -- time taken {str(elapsed)}')
    return represent_corrected_fasta, groups
def run_diamond_with_map(faa_file, map_file, out_dir, cover=90, threads=4):        
    starttime = datetime.now()
    
    diamond_represent_fasta= os.path.join(out_dir, 'diamond_rep_seq.fasta')
    diamond_cluster_file=os.path.join(out_dir, 'diamond_cluster.tsv')
    cmd = f'./diamond deepclust -d {faa_file} -o {diamond_cluster_file}  --approx-id 98 --round-approx-id 98 --round-coverage {cover} --member-cover {cover} --mutual-cover {cover} > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond')        
    cmd=f'CMD="seqtk subseq {faa_file} <(cut -f1 {diamond_cluster_file} | uniq) > {diamond_represent_fasta}" ; /bin/bash -c "$CMD"'
    #ret = run_command(cmd)
    ret = os.system(cmd)
    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond with 98% identity part 1 -- time taken {elapsed}')

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    represent_corrected_fasta = os.path.join(out_dir, 'diamond.fasta')
    with open(represent_corrected_fasta,'w') as ofh, open(diamond_represent_fasta) as ifh:
        for line in ifh:
            if line[0] == '>':
                ofh.write(f'>{gene_map[line[1:].strip()]}\n')
            else:
                ofh.write(line)      

    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond with 98% identity part 1 -- time taken {elapsed}')
    c_cursor=0
    cluster_count=0
    with open(diamond_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=cluster_count
                
                clusters[c_cursor] = {'gene_names':[]} 
                clusters[c_cursor]['representative'] = gene_map[member]
                cluster_count=cluster_count+1
                
            else:
                clusters[c_cursor]['gene_names'].append(gene_map[member])                

           
                
    
    del gene_map    

    # convert to a simple dictionary
    clusters_new = {}
    for cluster_name in clusters:
        clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond with 98% identity -- time taken {str(elapsed)}')
    return represent_corrected_fasta, clusters_new

def run_diamond_clustering(faa_file, map_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    #mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    diamond_cluster_file=os.path.join(out_dir, 'diamond_clusters')
    cmd = f'diamond linclust -d {faa_file} -o {diamond_cluster_file} --approx-id 98 --member-cover 98> /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond')        

    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond with 98% identity part 1 -- time taken {elapsed}')
   

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    #represent_corrected_fasta = os.path.join(out_dir, 'diamond.fasta')
    # with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
    #     for line in ifh:
    #         if line[0] == '>':
    #             ofh.write(f'>{gene_map[line[1:].strip()]}\n')
    #         else:
    #             ofh.write(line)      

    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond with 70% identity part 1 -- time taken {elapsed}')
    c_cursor=0
    cluster_count=0
    with open(diamond_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=cluster_count
                
                clusters[c_cursor] = {'gene_names':[]} 
                clusters[c_cursor]['representative'] = gene_map[member]
                cluster_count=cluster_count+1
                
            else:
                clusters[c_cursor]['gene_names'].append(gene_map[member])                

           
                
    
    del gene_map    

    # convert to a simple dictionary
   # clusters_new = {}
    inflated_clusters = []
    for cluster_name in clusters:
        #clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
        inflated_genes=[clusters[cluster_name]['representative']]
        inflated_genes.extend(clusters[cluster_name]['gene_names'])
        inflated_clusters.append(inflated_genes)
        

    elapsed = datetime.now() - starttime
    logging.info(f'Run  diamond 70% identity -- time taken {str(elapsed)}')
    return inflated_clusters
def run_diamond_clustering_pipeline(faa_file, map_file, out_dir, threads=4):
    starttime = datetime.now()
    
    #mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    diamond_cluster_file=os.path.join(out_dir, 'diamond_clusters')
    cmd = f'./diamond linclust -d {faa_file} -o {diamond_cluster_file} --approx-id 98 --member-cover 98> /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond')        

    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond with 98% identity part 1 -- time taken {elapsed}')
    diamond_represent_fasta= os.path.join(out_dir, 'diamond_rep_seq.faa')
    
    cmd=f'CMD="seqtk subseq {faa_file} <(cut -f1 {diamond_cluster_file} | uniq) > {diamond_represent_fasta}" ; /bin/bash -c "$CMD"'
    
    ret = os.system(cmd)
    elapsed = datetime.now() - starttime
    logging.info(f'Run seqtk subseq -- time taken {elapsed}')

    diamond_represent_fasta_db= os.path.join(out_dir, 'diamond_rep_seq_db')
    cmd=f'./diamond makedb --in {diamond_represent_fasta} -d {diamond_represent_fasta_db}'
    ret = run_command(cmd)
    #ret = os.system(cmd)
    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond makedb  -- time taken {elapsed}')
    out_blast=os.path.join(out_dir, 'blast_out')
    cmd=f'./diamond blastp -q {diamond_represent_fasta} -d {diamond_represent_fasta_db} -o {out_blast} --fast -f 6 qseqid sseqid qcovhsp scovhsp corrected_bitscore --approx-id 70 -p {threads} --query-cover 90 -k1000 '
    ret = run_command(cmd)
    #ret = os.system(cmd)
    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond blastp  -- time taken {elapsed}')
    
    edge_file=os.path.join(out_dir, 'edge.tsv')
    cmd=f'cat {out_blast} > {edge_file}'
    ret = run_command(cmd)
    #ret = os.system(cmd)
    elapsed = datetime.now() - starttime
    logging.info(f'Run cat blastp output  -- time taken {elapsed}')

    cmd=f'samtools faidx {diamond_represent_fasta}'
    ret = run_command(cmd)
    #ret = os.system(cmd)
    elapsed = datetime.now() - starttime
    logging.info(f'Run samtools faidx  -- time taken {elapsed}')
    diamond_represent_vertex= os.path.join(out_dir, 'diamond_rep_vertex')
    diamond_cluster_file_r2=os.path.join(out_dir, 'diamond_clusters_r2.tsv')
    cmd=f'./diamond greedy-vertex-cover --edges {edge_file} -d {diamond_represent_fasta}.fai --centroid-out {diamond_represent_vertex}  --edge-format triplet -o {diamond_cluster_file_r2}'
    ret = run_command(cmd)
    #ret = os.system(cmd)
    elapsed = datetime.now() - starttime
    logging.info(f'diamond greedy-vertex-cover  -- time taken {elapsed}')

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    #represent_corrected_fasta = os.path.join(out_dir, 'diamond.fasta')
    # with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
    #     for line in ifh:
    #         if line[0] == '>':
    #             ofh.write(f'>{gene_map[line[1:].strip()]}\n')
    #         else:
    #             ofh.write(line)      

    elapsed = datetime.now() - starttime
    logging.info(f'Run diamond clustering pipeline -- time taken {elapsed}')
    c_cursor=0
    cluster_count=0
    p_cluster='-1'
    with open(diamond_cluster_file_r2, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()

            if not rep_name == p_cluster:
                c_cursor=cluster_count
                
                clusters[c_cursor] = {'gene_names':[]} 
                clusters[c_cursor]['representative'] = gene_map[rep_name]
                cluster_count=cluster_count+1
                p_cluster=rep_name
                
            else:
                clusters[c_cursor]['gene_names'].append(gene_map[member])                

           
                
    
    del gene_map    

    # convert to a simple dictionary
   # clusters_new = {}
    inflated_clusters = []
    for cluster_name in clusters:
        #clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
        inflated_genes=[clusters[cluster_name]['representative']]
        inflated_genes.extend(clusters[cluster_name]['gene_names'])
        inflated_clusters.append(inflated_genes)
        

    elapsed = datetime.now() - starttime
    logging.info(f'Run  diamond 70% identity -- time taken {str(elapsed)}')
    return inflated_clusters
def run_mmseq_clustering(faa_file, map_file, out_dir, threads=4):        
    starttime = datetime.now()
    
    #mmseq_represent_fasta= os.path.join(out_dir, 'mmseq_rep_seq.fasta')
    mmseq_cluster_file=os.path.join(out_dir, 'mmseq_clusters')
    cmd = f'mmseqs easy-linclust {faa_file} {mmseq_cluster_file} {out_dir}/tmp --min-seq-id 0.70 -c 0.7  --threads {threads} > /dev/null'    
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseq')        
    mmseq_cluster_file=mmseq_cluster_file+"_cluster.tsv"
    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 70% identity  -- time taken {elapsed}')

    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    #represent_corrected_fasta = os.path.join(out_dir, 'diamond.fasta')
    # with open(represent_corrected_fasta,'w') as ofh, open(mmseq_represent_fasta) as ifh:
    #     for line in ifh:
    #         if line[0] == '>':
    #             ofh.write(f'>{gene_map[line[1:].strip()]}\n')
    #         else:
    #             ofh.write(line)      

    elapsed = datetime.now() - starttime
    logging.info(f'Run mmseq with 70% identity -- time taken {elapsed}')
    c_cursor=0
    cluster_count=0
    with open(mmseq_cluster_file, 'r') as fh:
        for line in fh:
            rep_name,member = line.strip().split()
            rep_name=rep_name.strip()
            member=member.strip()
           
            if rep_name == member:
                c_cursor=cluster_count
                
                clusters[c_cursor] = {'gene_names':[]} 
                clusters[c_cursor]['representative'] = gene_map[member]
                cluster_count=cluster_count+1
                
            else:
                clusters[c_cursor]['gene_names'].append(gene_map[member])                

           
                
    
    del gene_map    

    # convert to a simple dictionary
   # clusters_new = {}
    inflated_clusters = []
    for cluster_name in clusters:
        #clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
        inflated_genes=[clusters[cluster_name]['representative']]
        inflated_genes.extend(clusters[cluster_name]['gene_names'])
        inflated_clusters.append(inflated_genes)
        

    elapsed = datetime.now() - starttime
    logging.info(f'Run  mmseq 70% identity -- time taken {str(elapsed)}')
    return inflated_clusters

def run_blast(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4):
    starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    blast_db = os.path.join(out_dir, 'output_contigs')
    cmd = f"makeblastdb -in {database_fasta} -dbtype prot -out {blast_db} -logfile /dev/null"
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running makeblastdb')
    
    # chunk fasta file
    chunk_dir = os.path.join(out_dir, 'chunk_files')
    chunked_file_list = chunk_fasta_file(query_fasta, chunk_dir)

    # run parallel all-against-all blast
    #blast_cmds_file = os.path.join(out_dir,"blast_cmds.txt")    
    blast_output_file_list = []
    pool = multiprocessing.Pool(processes=threads)
    results = []      

    #with open(blast_cmds_file,'w') as fh:
    for chunked_file in chunked_file_list:
        blast_output_file = os.path.splitext(chunked_file)[0] + '.out'
        blast_output_file_list.append(blast_output_file)
        cmd = f'blastp -query {chunked_file} -db {blast_db} -evalue {evalue} -num_threads 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -max_target_seqs 2000 2> /dev/null 1> {blast_output_file}'
        results.append(pool.apply_async(run_command,(cmd, None)))
    pool.close()
    pool.join()
    
    for result in results:
        if result.get() != 0:
            raise Exception('Error running all-against-all blast')        
    
    # combining blast results
    blast_result = os.path.join(out_dir, 'blast_results')
    if os.path.isfile(blast_result):
        os.remove(blast_result)
    for blast_output_file in blast_output_file_list:
        os.system(f'cat {blast_output_file} >> {blast_result}')
        os.remove(blast_output_file)

    elapsed = datetime.now() - starttime
    logging.info(f'All-against-all BLASTP -- time taken {str(elapsed)}')
    return blast_result


def pairwise_alignment_diamond(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4,timing_log=None,max_seq=2000):
    starttime = datetime.now()
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # make diamond database
    diamond_db = os.path.join(out_dir, 'diamond_db')
    cmd = f'./diamond makedb --in {database_fasta} -d {diamond_db} -p {threads} --quiet'
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running diamond makedb')
    
    # run diamond blastp
    diamond_result = os.path.join(out_dir, 'diamond.tsv')
    cmd = f'./diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --max-target-seqs {max_seq} 2> /dev/null 1> {diamond_result}'
    #subprocess.call(cmd, shell=True)
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running diamond makedb')


    elapsed = datetime.now() - starttime
    logging.info(f'Protein pairwise alignment with Diamond -- time taken {str(elapsed)}')
    return diamond_result

def split_fasta(input_fasta, output_dir, batch_size):
    """Splits a FASTA file into chunks of batch_size sequences each."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    chunk_count = 0
    chunk_file = None
    seq_count = 0

    with open(input_fasta, "r") as infile:
        for line in infile:
            if line.startswith(">"):  # New sequence
                if seq_count >= batch_size:
                    chunk_file.close()
                    chunk_count += 1
                    seq_count = 0
                
                if seq_count == 0:
                    chunk_file = open(os.path.join(output_dir, f"chunk_{chunk_count}.fasta"), "w")

                seq_count += 1
            
            if chunk_file:
                chunk_file.write(line)
    
    if chunk_file:
        chunk_file.close()
    
    return [os.path.join(output_dir, f"chunk_{i}.fasta") for i in range(chunk_count + 1)]

def pairwise_alignment_diamond_split_db(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4, timing_log=None, max_seq=2000, batch_size=100000):
    starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Step 1: Split database FASTA into chunks
    db_chunks_dir = os.path.join(out_dir, "db_chunks")
    db_chunks = split_fasta(database_fasta, db_chunks_dir, batch_size)

    all_results = []
    
    for i, chunk in enumerate(db_chunks):
        logging.info(f"Processing database chunk {i+1}/{len(db_chunks)}: {chunk}")
        
        # Create a DIAMOND database for the chunk
        diamond_db = os.path.join(out_dir, f'diamond_db_chunk_{i}')
        cmd = f'./diamond makedb --in {chunk} -d {diamond_db} -p {threads} --quiet'
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception(f'Error running diamond makedb for chunk {i}')
        
        # Run DIAMOND BLASTP for this chunk
        chunk_result = os.path.join(out_dir, f'diamond_chunk_{i}.tsv')
        cmd = f'./diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --max-target-seqs {max_seq} 2> /dev/null 1> {chunk_result}'
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception(f'Error running diamond blastp for chunk {i}')
        
        all_results.append(chunk_result)

    # Step 3: Merge all chunk results into a single file
    final_result_file = os.path.join(out_dir, "diamond_final.tsv")
    with open(final_result_file, "w") as outfile:
        for i, result_file in enumerate(all_results):
            with open(result_file, "r") as infile:
                if i == 0:  # Copy header from the first file
                    outfile.write(infile.read())
                else:  # Skip headers from subsequent files
                    next(infile)  # Skip first line
                    outfile.write(infile.read())

    # Clean up temporary chunk files and databases
    shutil.rmtree(db_chunks_dir)
    # for chunk in db_chunks:
    #     os.remove(chunk)
    # for i in range(len(db_chunks)):
    #     os.remove(os.path.join(out_dir, f'diamond_db_chunk_{i}.dmnd'))
    #     os.remove(os.path.join(out_dir, f'diamond_chunk_{i}.tsv'))

    elapsed = datetime.now() - starttime
    logging.info(f'Protein pairwise alignment with Diamond completed -- time taken {str(elapsed)}')
    
    return final_result_file
#import os
#import subprocess
from collections import defaultdict
#from Bio import SeqIO

def diamond_cd_hit_2d(reference_fasta, query_fasta, out_dir,threads,evalue,timing_log=None,identity=98, coverage=90):
    """
    Mimics cd-hit-2d using DIAMOND to cluster query sequences against a reference set.
    
    Parameters:
    - reference_fasta: Path to reference protein sequences (FASTA)
    - query_fasta: Path to query protein sequences (FASTA)
    - output_clusters: Path to save clusters (txt)
    - output_unmatched: Path to save unmatched sequences (FASTA)
    - identity: Minimum percent identity for clustering (default 90%)
    - coverage: Minimum query and subject coverage for clustering (default 80%)
    
    Output:
    - A text file with clusters where the representative comes from the reference.
    - A FASTA file with query sequences that did not match any reference sequence.
    """
    
    # Step 1: Create DIAMOND database from reference
   
    starttime = datetime.now()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    diamond_db = os.path.join(out_dir, 'diamond_ref_db')
    cmd = f'./diamond makedb --in {reference_fasta} -d {diamond_db} -p {threads} --quiet'
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running diamond makedb')
    
    # Step 2: Run DIAMOND alignment
    blast_output = os.path.join(out_dir,"diamond_matches.m8")
    cmd = f'./diamond blastp -q {query_fasta} -d {diamond_db} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --max-target-seqs 1 --id {str(identity)} --query-cover {str(coverage)} --subject-cover {str(coverage)} 2> /dev/null 1> {blast_output}'
    #subprocess.call(cmd, shell=True)
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running diamond makedb')

    #subprocess.run([
    #    "diamond", "blastp", "-d", db_path, "-q", query_fasta, "-o", blast_output,
    #    "--id", str(identity), "--query-cover", str(coverage), "--subject-cover", str(coverage),
    #    "--max-target-seqs", "1", "--outfmt", "6"
    #], check=True)

    # Step 3: Process results to form clusters
    clusters = defaultdict(list)
    matched_queries = set()

    with open(blast_output, "r") as f:
        for line in f:
            query_id, ref_id = line.split("\t")[:2]
            clusters[ref_id].append(query_id)
            matched_queries.add(query_id)

    
    # Step 5: Extract unmatched sequences
    unmatched_sequences = []
    for record in SeqIO.parse(query_fasta, "fasta"):
        if record.id not in matched_queries:
            unmatched_sequences.append(record)
    output_unmatched=os.path.join(out_dir,"unmatched_new_unique_seq.fasta")
    SeqIO.write(unmatched_sequences, output_unmatched, "fasta")

    # Cleanup
    #os.remove(blast_output)
    #os.remove(db_path)

    elapsed = datetime.now() - starttime
    logging.info(f'Clustering complete, merge {len(matched_queries)} similar unique seqs to existed groups, remain {len(unmatched_sequences)} seqs, unmatched sequences saved to {output_unmatched} -- time taken {str(elapsed)}')
    
    return clusters,output_unmatched
# Example usage:
# diamond_cd_hit_2d("reference.fasta", "query.fasta", "clusters.txt", "unmatched.fasta")

import sourmash
from sourmash import MinHash, SourmashSignature
def pairwise_alignment_sourmash( query_fasta, out_dir,ksize=3, similarity=0.7,diff_len=0.7,num=100, evalue=1E-6, threads=4,timing_log=None):
    starttime = datetime.now()
   
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    signatures = []
    sequence_names = []
    seq_len=[]
    for record in SeqIO.parse(query_fasta, "fasta"):
        sequence_names.append(record.id)
        # Create a MinHash sketch for the current sequence
        mh = MinHash(ksize=ksize, n=num,  is_protein=True)
        mh.add_sequence(str(record.seq), force=True)
        # Create a SourmashSignature for this MinHash
        sig = SourmashSignature(mh, name=record.id)
        signatures.append(sig)
        seq_len.append(len(record.seq))
    
    
    num_sequences = len(signatures)
    
    output_tsv = os.path.join(out_dir, 'distances_sourmash.tsv')
    with open(output_tsv, "w", newline="") as f:
        for i in range(num_sequences-1):
            for j in range(i,num_sequences):
                sim=signatures[i].similarity(signatures[j])
                #diff=abs(seq_len[i]-seq_len[j])
                p_diff=min(seq_len[i],seq_len[j])/max(seq_len[i],seq_len[j])
                if sim > similarity and p_diff>diff_len: 
                    str_line=sequence_names[i]+"\t"+sequence_names[j]+"\t"+str(sim)
                    f.write(str_line+"\n")
    
   
    
        
  
    elapsed = datetime.now() - starttime
    logging.info(f'Protein pairwise alignment with sourmash -- time taken {str(elapsed)}')
    return output_tsv
def pairwise_alignment_mash(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4,timing_log=None):
    starttime = datetime.now()
    mem_usage=0
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    cmd = f'mash sketch -o database_sketch  {database_fasta}'
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mash')
    cmd = f'mash sketch -o query_sketch  {query_fasta}'
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mash sketch')
    mash_result = os.path.join(out_dir, 'distances_mash.tsv')
    cmd = f'mash dist database_sketch.msh query_sketch.msh >  {mash_result}'
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mash sketch')
  
    elapsed = datetime.now() - starttime
    logging.info(f'Protein pairwise alignment with faiss -- time taken {str(elapsed)}')
    return mash_result
def pairwise_alignment_faiss(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4,timing_log=None):
    starttime = datetime.now()
    mem_usage=0
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    index_seq,vectors=protein_sequences_to_2mer(database_fasta)
    array_input = np.array(vectors)
    transformed_data = transform_to_l1_space(array_input)
    index = faiss.IndexFlatL2(transformed_data.shape[1])  # Sử dụng FAISS với L2 Distance
    index.add(transformed_data) 
    mem_usage = mem_report(mem_usage, "faiss input")
    elapsed = datetime.now() - starttime
    logging.info(f'load faiss index -- time taken {elapsed}')
    distances, indices = index.search(transformed_data, k=2000)
    row_sums = array_input.sum(axis=1, keepdims=True)
    distances = distances / row_sums
    faiss_result = os.path.join(out_dir, 'distances_faiss.tsv')
    visited={}
    for i in range(len(transformed_data)):       
        near_items=np.where(distances[i] <= 4)[0]
        if len(near_items)>1:
            for j in near_items:
                if i!=indices[i][j]:
                    ind=str(i)+","+str(indices[i][j])
                    rev_ind=str(indices[i][j])+","+str(i)
                    if ind not in visited and rev_ind not in visited:
                        visited[ind]=distances[i][j]
    with open(faiss_result,'w') as f:
        for k in visited.keys():
            nodes=k.split(",")
            str_line=index_seq[int(nodes[0])]+" "+index_seq[int(nodes[1])]+" "+str(visited[k])
            f.write(str_line+"\n")
    lapsed = datetime.now() - starttime
    logging.info(f'Protein pairwise alignment with faiss -- time taken {str(elapsed)}')
    return faiss_result
def pairwise_alignment_ems(database_fasta, query_fasta, out_dir, evalue=1E-6, threads=4,timing_log=None):
    starttime = datetime.now()
    mem_usage=0
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    index_seq,vectors=protein_sequences_to_vector_ems(database_fasta)
    mem_usage = mem_report(mem_usage, "ems vector")
    ems_result = os.path.join(out_dir, 'distances_ems.tsv')
    
    # with open(ems_result,'w') as f:
    #     for i in range(len(vectors)):
    #         for j in range(i+1,len(vectors)):
    #             d=l2_distance(vectors[i],vectors[j])
    #             if d<150:
    #                 str_line=index_seq[i]+" "+index_seq[j]+" "+str(d)
    #                 f.write(str_line+"\n")
    #faiss distance
    # array_input = np.array(vectors)
    
    # index = faiss.IndexFlatL2(array_input.shape[1])  # Sử dụng FAISS với L2 Distance
    # index.add(array_input) 
    # mem_usage = mem_report(mem_usage, "faiss input")
    # elapsed = datetime.now() - starttime
    # logging.info(f'load faiss index -- time taken {elapsed}')
    # distances, indices = index.search(array_input, k=2000)
    # row_sums = array_input.sum(axis=1, keepdims=True)
    # #distances = distances / row_sums
    # #faiss_result = os.path.join(out_dir, 'distances_faiss.tsv')
    # visited={}
    # for i in range(len(array_input)):       
    #     near_items=np.where(distances[i] <= 200)[0]
    #     if len(near_items)>1:
    #         for j in near_items:
    #             if i!=indices[i][j]:
    #                 ind=str(i)+","+str(indices[i][j])
    #                 rev_ind=str(indices[i][j])+","+str(i)
    #                 if ind not in visited and rev_ind not in visited:
    #                     visited[ind]=distances[i][j]
    # with open(ems_result,'w') as f:
    #     for k in visited.keys():
    #         nodes=k.split(",")
    #         str_line=index_seq[int(nodes[0])]+" "+index_seq[int(nodes[1])]+" "+str(visited[k])
    #         f.write(str_line+"\n")
    # cosin distance
    visited={}
    distances=calculate_cosine_distances_np(vectors)
    for i in range(len(vectors)):       
        near_items=np.where(distances[i] <= 0.2)[0]
        if len(near_items)>1:
            for j in near_items:
                if i!=j:
                    ind=str(i)+","+str(j)
                    rev_ind=str(j)+","+str(i)
                    if ind not in visited and rev_ind not in visited:
                        visited[ind]=distances[i][j]
    with open(ems_result,'w') as f:
        for k in visited.keys():
            nodes=k.split(",")
            str_line=index_seq[int(nodes[0])]+" "+index_seq[int(nodes[1])]+" "+str(visited[k])
            f.write(str_line+"\n")
    mem_usage = mem_report(mem_usage, "cal distance")
    elapsed = datetime.now() - starttime
    logging.info(f'Protein pairwise alignment with ems -- time taken {str(elapsed)}')
    return ems_result
def pairwise_alignment_partion_diamond(input_fasta, out_dir, threshold=0.7,evalue=1E-6, threads=4):
    starttime = datetime.now()
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    list_len=[]
    with open(input_fasta) as inp:
        for record in SeqIO.parse(inp,'fasta'):
            list_len.append({'id':record.id,'len':len(record.seq)})
           
    list_len.sort(key=lambda x: x['len'], reverse=False)
     
    next_mark=0
    nm=[]
    map_name_order={}
    for i in range(0,len(list_len)):
        map_name_order[list_len[i]['id']]=i
        k=list_len[i]['len']
        if next_mark==0:
            next_mark=k+k*(1-threshold)
        if k>next_mark:
            next_mark=k+k*(1-threshold)
            nm.append(next_mark)
    logging.info(nm)    
    #partion first round:
    index=0
    map_to_partition={}
    c_nm=[]
    c_nm.extend(nm)
    number_mark=len(nm)
    for i in range(0, len(list_len)):
        if len(nm)<=0:
            break
        if list_len[i]['len']<=nm[0]:
            map_to_partition[list_len[i]['id']]=index
        else:
            del nm[0]
            index=index+1
            map_to_partition[list_len[i]['id']]=index
    
    #seq_files_fasta=[]
    with open(input_fasta) as inp:
        for record in SeqIO.parse(inp,'fasta'):
            fasta_file=os.path.join(out_dir,f'p_seq{map_to_partition[record.id]}.fasta')
            #seq_files_fasta.append(fasta_file)
            with open(fasta_file, 'a') as out_fh:
                out_fh.write(SeqIO.FastaIO.as_fasta(record))
    #round 2:
    #seq_files_fasta_overlap=[]
    nm2=[]
    for i in range(1,len(c_nm)-1):
        nm2.append((c_nm[i]+c_nm[i+1])/2)
    logging.info(nm2)    
    index=0
    map_to_partition={}
    number_mark2=len(nm2)
    for i in range(0, len(list_len)):
        if len(nm2)<=0:
            break
        if list_len[i]['len']<=nm2[0]:
            map_to_partition[list_len[i]['id']]=index
        else:
            del nm2[0]
            index=index+1
            map_to_partition[list_len[i]['id']]=index
    with open(input_fasta) as inp:
        for record in SeqIO.parse(inp,'fasta'):
            if not record.id in map_to_partition:
                break
            fasta_file=os.path.join(out_dir,f'op_seq{map_to_partition[record.id]}.fasta')
            #seq_files_fasta_overlap.append(fasta_file)
            with open(fasta_file, 'a') as out_fh:
                out_fh.write(SeqIO.FastaIO.as_fasta(record))
    json.dump(map_to_partition, open(os.path.join(out_dir, 'map_to_partion2.json'), 'w'), indent=4, sort_keys=True) 
    
    for f_index in range(0,number_mark):
         # make diamond database
        diamond_db = os.path.join(out_dir, f'diamond_db{f_index}')
        fasta_file=os.path.join(out_dir,f'p_seq{f_index}.fasta')
        count_seq_in_partion=0
        with open(fasta_file) as inp:
            for record in SeqIO.parse(inp,'fasta'):
                count_seq_in_partion=count_seq_in_partion+1
        logging.info(f'count seq in {fasta_file} : {count_seq_in_partion}')   
        cmd = f'./diamond makedb --in {fasta_file} -d {diamond_db} -p {threads} --quiet'
        #ret = os.system(cmd)
        ret = run_command(cmd)
        if ret != 0:
            raise Exception('Error running diamond makedb')
        
        # run diamond blastp
        diamond_result = os.path.join(out_dir, f'diamond{f_index}.tsv')
        cmd = f'./diamond blastp -q {fasta_file} -d {diamond_db} -p {threads} -b1 --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --max-target-seqs 2000 > {diamond_result}'
        #subprocess.call(cmd, shell=True)
        #ret = os.system(cmd)
        gc.collect()
        ret = run_command(cmd)
        if ret != 0:
            raise Exception('Error running diamond blastp'+str(ret))
    for f_index in range(0,number_mark2):
         # make diamond database
        fasta_file=os.path.join(out_dir,f'op_seq{f_index}.fasta')
        count_seq_in_partion=0
        if not os.path.exists(fasta_file):
            continue
        with open(fasta_file) as inp:
            for record in SeqIO.parse(inp,'fasta'):
                count_seq_in_partion=count_seq_in_partion+1
        logging.info(f'count seq in {fasta_file} : {count_seq_in_partion}')  
        diamond_db = os.path.join(out_dir, f'diamond_db_o{f_index}')
        cmd = f'./diamond makedb --in {fasta_file} -d {diamond_db} -p {threads} --quiet'
        #ret = os.system(cmd)
        ret = run_command(cmd)
        gc.collect()
        if ret != 0:
            raise Exception('Error running diamond makedb')
    
        # run diamond blastp
        diamond_result = os.path.join(out_dir, f'diamond_o{f_index}.tsv')
        cmd = f'./diamond blastp -q {fasta_file} -d {diamond_db} -b1 -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --max-target-seqs 2000 > {diamond_result}'
        #subprocess.call(cmd, shell=True)
        #ret = os.system(cmd)
        ret = run_command(cmd)
        if ret != 0:
            raise Exception('Error running diamond blastp' + str(ret))    
    set_pairwise_ident=set()
    combine_diamond_tsv=os.path.join(out_dir, f'diamond.tsv')
    cdt=open(combine_diamond_tsv,'w')
    for f_index in range(0,number_mark):
        set_pairwise_ident=set()
        with open(os.path.join(out_dir, f'diamond{f_index}.tsv')) as file_tsv:
            for line in file_tsv:
                t=line.split('\t')
                couple=str(map_name_order[t[0]])+'_'+str(map_name_order[t[1]])
                if not couple in set_pairwise_ident:
                    cdt.write(line)
                    set_pairwise_ident.add(couple)
        if f_index < number_mark-1 and os.path.exists(os.path.join(out_dir, f'diamond_o{f_index}.tsv')):
            
                
            with open(os.path.join(out_dir, f'diamond_o{f_index}.tsv')) as file_tsv:
                for line in file_tsv:
                    t=line.split('\t')
                    couple=str(map_name_order[t[0]])+'_'+str(map_name_order[t[1]])
                    if not couple in set_pairwise_ident:
                        cdt.write(line)
                        set_pairwise_ident.add(couple)
        if f_index > 0 and os.path.exists(os.path.join(out_dir, f'diamond_o{f_index-1}.tsv')):
            with open(os.path.join(out_dir, f'diamond_o{f_index-1}.tsv')) as file_tsv:
                for line in file_tsv:
                    t=line.split('\t')
                    couple=str(map_name_order[t[0]])+'_'+str(map_name_order[t[1]])
                    if not couple in set_pairwise_ident:
                        cdt.write(line)
                        set_pairwise_ident.add(couple)

    cdt.close()

    elapsed = datetime.now() - starttime
    logging.info(f'Protein partition pairwise alignment with Diamond -- time taken {str(elapsed)}')
    return  combine_diamond_tsv

def pairwise_alignment_mmseq(database_fasta, query_fasta, out_dir, evalue=1E-6 ,threads=4):
    starttime = datetime.now()
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # make mmseq database
    mmseq_db = os.path.join(out_dir, 'mmseq_db')
    cmd = f'mmseqs createdb {database_fasta} {mmseq_db}'
    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseqs createdb target_db')
    query_db = os.path.join(out_dir, 'query_db')
    cmd = f'mmseqs createdb {query_fasta} {query_db}'
    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseqs createdb query_db')
    cmd = f'mmseqs createindex {mmseq_db} tmp'
    #ret = os.system(cmd)
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running mmseqs createindex ')
    
    # run mmseq blastp
    mmseq_result = os.path.join(out_dir, 'mmseqs.tsv')
    cmd = f'mmseqs search {query_db} {mmseq_db} resultDB tmp'
    #subprocess.call(cmd, shell=True)
    #ret = os.system(cmd)
    ret = run_command(cmd)

    if ret != 0:
        raise Exception('Error running mmseqs search')
    cmd = f'mmseqs convertalis {query_db} {mmseq_db} resultDB -outfmt 6 {mmseq_result}'
    #subprocess.call(cmd, shell=True)
    #ret = os.system(cmd)
    ret = run_command(cmd)

    if ret != 0:
        raise Exception('Error running mmseqs convertalis')
    
    elapsed = datetime.now() - starttime
    logging.info(f'Protein alignment with mmseq -- time taken {str(elapsed)}')
    return mmseq_result



def filter_blast_result(blast_result, 
                        # gene_annotation, 
                        out_dir, identity, length_difference, alignment_coverage_short, alignment_coverage_long):
    filtered_blast_result_file = os.path.join(out_dir, 'filtered_blast_results')

    with open(filtered_blast_result_file, 'w') as fh:
        for line in open(blast_result, 'r'):
            cells = line.rstrip().split('\t')            

            qlen = int(cells[12])# * 3 + 3
            slen = int(cells[13])# * 3 + 3

            pident = float(cells[2]) / 100
            alignment_length = int(cells[3]) # * 3

            short_seq = min(qlen, slen)
            long_seq = max(qlen, slen)
            len_diff = short_seq / long_seq
            align_short = alignment_length / short_seq
            align_long = alignment_length / long_seq
            
            if pident <= identity or len_diff <= length_difference or align_short <= alignment_coverage_short or align_long <= alignment_coverage_long:
                continue

            fh.write(line)

    return filtered_blast_result_file

            
def cluster_with_mcl(blast_result, out_dir, threads=4,inflation=1.5,timing_log=None):
    starttime = datetime.now()
    if threads > 1:
        threads = threads - 1
    
    mcl_file = os.path.join(out_dir, 'mcl_clusters')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_result} 2> /dev/null | mcl - --abc -I {inflation} -te {threads} -o {mcl_file} > /dev/null 2>&1"
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    elapsed = datetime.now() - starttime
    logging.info(f'Cluster with MCL -- time taken {str(elapsed)}')
    return mcl_file
def cluster_with_mcl_from_faiss(file_faiss, out_dir, threads=4,inflation=1.5,timing_log=None):
    starttime = datetime.now()
    if threads > 1:
        threads = threads - 1
    
    mcl_file = os.path.join(out_dir, 'mcl_clusters')
    cmd = f"mcl {file_faiss} --abc -I {inflation} -te {threads} -o {mcl_file} > /dev/null 2>&1"
    #ret = os.system(cmd)
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    elapsed = datetime.now() - starttime
    logging.info(f'Cluster with MCL from faiss -- time taken {str(elapsed)}')
    return mcl_file

def reinflate_clusters(groups, mcl_file):
    """
    Return
    ------
        - inflated_clusters: list of list of genes
        -clusters: dict(cluster_id->[gene_ids])
    """    
    starttime = datetime.now()
    clusters = {}
    clusters.update(groups)

    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    count_unique_seq_in_mcl=0
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = {}
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                #inflated_genes.append(gene)
                if gene in groups:
                    #inflated_genes.extend(groups[gene]['gene_id'])
                    inflated_genes[gene]=[]
                    inflated_genes[gene].append(gene)
                    inflated_genes[gene].extend(groups[gene])
                    count_unique_seq_in_mcl=count_unique_seq_in_mcl+len(groups[gene])+1
                    del groups[gene]
            inflated_clusters.append(inflated_genes)
    
    # Inflate any clusters that were in the clusters file but not sent to mcl
    count_not_mcl=0
    for gene in groups:
        count_not_mcl=count_not_mcl+1+len(groups[gene])

        inflated_genes = {}
        inflated_genes[gene]=[]
        inflated_genes[gene].append(gene)
        #inflated_genes.extend(groups[gene]['gene_id'])
        inflated_genes[gene].extend(groups[gene])
        inflated_clusters.append(inflated_genes)
    #count seq in inflate
    set_unique_seq=set()
    for c in inflated_clusters:
        for k in c:
            set_unique_seq.update(c[k])
    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate new {len(inflated_clusters)} clusters with {count_unique_seq_in_mcl} unique_seq (group), {count_not_mcl} seqs (groups) not found in MCL clustering , {len(set_unique_seq)} seq in inflate clusters-- time taken {str(elapsed)}')
    return inflated_clusters, clusters
def reinflate_clusters_with_unique_seqs(unique_groups,similar_groups, mcl_file):
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
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = {}
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                #inflated_genes.append(gene)
                if gene in similar_groups:
                    #inflated_genes.extend(groups[gene]['gene_id'])

                    inflated_genes[gene]=[]
                    inflated_genes[gene].append(gene)
                    inflated_genes[gene].extend(unique_groups[gene])
                    for g in similar_groups[gene]:
                        inflated_genes[gene].append(g)
                        inflated_genes[gene].extend(unique_groups[g])
                        

                    del similar_groups[gene]
            inflated_clusters.append(inflated_genes)
    
    # Inflate any clusters that were in the clusters file but not sent to mcl
    count_not_mcl=0
    for gene in similar_groups:
        count_not_mcl=count_not_mcl+1

        inflated_genes = {}
        inflated_genes[gene]=[]
        inflated_genes[gene].append(gene)
        inflated_genes[gene].extend(unique_groups[gene])
        for g in similar_groups[gene]:
            inflated_genes[gene].append(g)
            inflated_genes[gene].extend(unique_groups[g])
                    
        inflated_clusters.append(inflated_genes)
    
    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate new {len(inflated_clusters)} clusters with refer to {len(clusters.keys())} groups, {count_not_mcl} groups not found in MCL clustering -- time taken {str(elapsed)}')
    return inflated_clusters, clusters
def reinflate_clusters_by_groups(groups, mcl_file):
    """
    Return
    ------
        - inflated_clusters: list of list of genes
        -clusters: dict(cluster_id->[gene_ids])
    """    
    starttime = datetime.now()
    clusters = {}
    clusters.update(groups)

    inflated_clusters = []
    # Inflate genes from cdhit which were sent to mcl
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                #inflated_genes.append(gene)
                if gene in groups:
                    groups[gene]['confident']=1
                    inflated_genes.append(groups[gene])
                    del groups[gene]
            inflated_clusters.append(inflated_genes)
    
    # Inflate any clusters that were in the clusters file but not sent to mcl
    count_not_mcl=0
    for gene in groups:
        count_not_mcl=count_not_mcl+1

        inflated_genes = []
        #inflated_genes.append(gene)
        groups[gene]['confident']=1
        inflated_genes.append(groups[gene])
        inflated_clusters.append(inflated_genes)
    
    elapsed = datetime.now() - starttime
    logging.info(f'Reinflate new {len(inflated_clusters)} clusters with refer to {len(clusters.keys())} groups, {count_not_mcl} groups not found in MCL clustering -- time taken {str(elapsed)}')
    return inflated_clusters, clusters
def make_clusters_from_mcl(mcl_file, map_file):
    clusters = {}
    gene_map = {}
    count = 0
    with open(map_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            gene_map[f'{count}'] = line
            count += 1
    c_cursor=0
    cluster_count=0
    with open(mcl_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            genes = line.split('\t')
            c_cursor=cluster_count
            clusters[c_cursor] = {'gene_names':[]} 
            clusters[c_cursor]['representative'] = gene_map[genes[0]]
            cluster_count=cluster_count+1
            for i in range(1,len(genes)):
                clusters[c_cursor]['gene_names'].append(gene_map[genes[i]])                
                   
                            

           
                
    
    del gene_map    

    # convert to a simple dictionary
   # clusters_new = {}
    inflated_clusters = []
    for cluster_name in clusters:
        #clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']    
        inflated_genes=[clusters[cluster_name]['representative']]
        inflated_genes.extend(clusters[cluster_name]['gene_names'])
        inflated_clusters.append(inflated_genes)
    return inflated_clusters, clusters
def match_seqs_to_ref(out_dir,seqs_file, groups, refdb,ref_clusters, threads=1, evalue=1E-6):
    starttime = datetime.now()
    diamond_result = os.path.join(out_dir, 'ref_diamond.tsv')
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
        if not cells[0] in top_match.keys():
            top_match[cells[0]]={'len':qlen,'cluster':clustername,'ident':pident} 
        if top_match[cells[0]]['ident']<pident:
            top_match[cells[0]]['cluster']=clustername
            top_match[cells[0]]['ident']=pident
    un_match_combined_faa_file = os.path.join(out_dir, 'unmatch_combined.faa')
    clusters={}
    # gene_map = {}
    # count = 0
    # with open(map_file, 'r') as fh:
    #     for line in fh:
    #         line = line.strip()
    #         gene_map[f'{line}'] = count
    #         count += 1
    matched_cluster=set()
    for sid in top_match.keys():
        c=top_match[sid]['cluster']
        if not c in clusters.keys():
            
            clusters[c]={}
            clusters[c]['gene_id']=[]
            clusters[c]['max_length']=0
            clusters[c]['mean_length']=0
            clusters[c]['min_length']=1E6
            clusters[c]['gene_name']=ref_clusters[c]['gene_name']
            clusters[c]['product']=ref_clusters[c]['description']
            clusters[c]['representative']=''
            clusters[c]['unique_seq']=[]
            clusters[c]['source']='reference'
            clusters[c]['size']=0
            matched_cluster.add(c)
        clusters[c]['unique_seq'].append(sid)
        clusters[c]['gene_id'].append(sid)
        clusters[c]['gene_id'].extend(groups[sid])
        new_size=clusters[c]['size']+1+len(groups[sid])
        clusters[c]['size']=new_size
        clusters[c]['mean_length']=float((int(clusters[c]['mean_length'])*(new_size-1-len(groups[sid]))+int(top_match[sid]['len'])))/(new_size-len(groups[sid]))
        clusters[c]['max_length']=max(int(clusters[c]['max_length']),int(top_match[sid]['len']))
        clusters[c]['min_length']=min(int(clusters[c]['min_length']),int(top_match[sid]['len']))
        del groups[sid]
    """ for c in ref_clusters.keys():
        if not c in clusters.keys():
            clusters[c]={}
            clusters[c]['gene_id']=[]
            clusters[c]['max_length']=0
            clusters[c]['mean_length']=0
            clusters[c]['min_length']=1E6
            clusters[c]['product']=ref_clusters[c]['description']
            clusters[c]['representative']=ref_clusters[c]['representative']
            clusters[c]['size']=0 """
    
    json.dump(clusters, open(os.path.join(out_dir, 'clusters_by_ref.json'), 'w'), indent=4, sort_keys=True)
    count_unmatched=0
    count_num_input=0
    with open(un_match_combined_faa_file, 'w') as fh, open(seqs_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            count_num_input=count_num_input+1
            if not newseq.id in top_match.keys():
                #newseq.id=str(gene_map[newseq.id])
                newseq.description=''
                newseq.name=''
                SeqIO.write(newseq,fh,'fasta')
                count_unmatched=count_unmatched+1
    #del gene_map 
    elapsed = datetime.now() - starttime
    logging.info(f'Remain {len(groups.keys())} not matched')
    logging.info(f'Matching {count_num_input} groups to ref clusters, {len(top_match.keys())} matched, form {len(clusters.keys())} ref clusters,  {count_unmatched} groups not matched  -- time taken {str(elapsed)}')
    return un_match_combined_faa_file,clusters,groups
def match_seqs_to_ref_by_nearest_group(out_dir,seqs_file, groups, refdb,refcdb,ref_clusters, threads=1, evalue=1E-6):
    starttime = datetime.now()
    diamond_result = os.path.join(out_dir, 'ref_diamond.tsv')
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
        if not sid in top_match.keys():
            top_match[sid]=[]
        top_match[sid].append({'len':qlen,'cluster':clustername,'ident':pident} )
       

        
    un_match_combined_faa_file = os.path.join(out_dir, 'unmatch_combined.faa')
    clusters={}
    # gene_map = {}
    # count = 0
    # with open(map_file, 'r') as fh:
    #     for line in fh:
    #         line = line.strip()
    #         gene_map[f'{line}'] = count
    #         count += 1
    matched_cluster=set()
    dict_seqs={}
    temseqdir=os.path.join(out_dir, 'temp_seqs')
    temmatchingdir=os.path.join(out_dir, 'temp_matching')
    os.mkdir(temseqdir)
    os.mkdir(temmatchingdir)
    with open(seqs_file, 'r') as fh:
        for seq in SeqIO.parse(fh,'fasta'):
            temp_seq= os.path.join(temseqdir, seq.id+".faa")
            with open(temp_seq, 'w') as fo:
                SeqIO.write(seq,fo,'fasta')
        #blast with groups in matched clusters
    pool = multiprocessing.Pool(processes=threads)
    results = []
    for sid in top_match.keys():
        #print(sid +' check candidated clusters')
        os.mkdir(temmatchingdir+"/"+sid)
        
        for c in top_match[sid]:
            diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
            cmd = f'./diamond blastp --quiet -q {os.path.join(temseqdir,sid+".faa")} -d {refcdb+"/"+c["cluster"]+".db.dmnd"} -p 1 --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamondout}'
            results.append(pool.apply_async(run_command,(cmd, None)))
        
    pool.close()
    pool.join()
    for result in results:
        if result.get() != 0:
            raise Exception('Error running diamond with ref clusters')
    
    added_group=set()
    for sid in top_match.keys():
        neartest_cluster=top_match[sid][0]
        best_ident=0
        #groups[sid]['confident']=1
        list_matched_groups=[]
        hash_clustername={}
        hash_cluster_max_ident={}
        
        for c in top_match[sid]:
            diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
            #read and note max identity group
            hash_cluster_max_ident[c['cluster']]={'pident':0,'len_diff':0,'align_short':0,'align_long':0}
            for line in open(diamondout, 'r'):
                cells = line.rstrip().split('\t')            

                sid=cells[0]
                groupname=cells[1]
                pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
                if pident < 0.98 or len_diff < 0.98 or align_short < 0.98 or align_long < 0.98:
                    continue
                list_matched_groups.append({'c':c,'g':groupname,'pid':pident})
                        
                if pident>best_ident:
                    best_ident=pident
                    neartest_cluster=c
                    
                if pident>hash_cluster_max_ident[c['cluster']]['pident']:
                    hash_cluster_max_ident[c['cluster']]={'match_group':groupname,'pident':pident,'len_diff':len_diff,'align_short':align_short,'align_long':align_long}
            hash_clustername[c['cluster']]=c
            
        if best_ident<0.98:
            #give up, not enough envident to add to any ref clusters
            
            continue
        # else:
        #     #if group not trully belong to a ref cluster, try to add it to best cluster by k-nearest neighbor and cal confident 
        #     #print("consider the undetermined :")
        #     list_matched_groups.sort(key=lambda x: x['pid'], reverse=True)
        #     print(list_matched_groups)
        #     count_matched_cluster={}         
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

        else:

            groups[sid]['confident']=1
            groups[sid]['match']=hash_cluster_max_ident[neartest_cluster['cluster']]
        nearest_cluster_name=neartest_cluster['cluster']
        if not nearest_cluster_name in clusters.keys():
            
            clusters[nearest_cluster_name]={}
            clusters[nearest_cluster_name]['groups']=[]
            clusters[nearest_cluster_name]['max_length']=0
            clusters[nearest_cluster_name]['mean_length']=0
            clusters[nearest_cluster_name]['min_length']=1E6
            clusters[nearest_cluster_name]['gene_name']=ref_clusters[nearest_cluster_name]['gene_name']
            clusters[nearest_cluster_name]['product']=ref_clusters[nearest_cluster_name]['description']
            clusters[nearest_cluster_name]['representative']='gene_families/sequences/'+nearest_cluster_name+'.fasta'
            clusters[nearest_cluster_name]['size']=0
            clusters[nearest_cluster_name]['source']='reference'
            matched_cluster.add(nearest_cluster_name)
        #groups[sid]['confident']=1
        clusters[nearest_cluster_name]['groups'].append(groups[sid])
        #clusters[neartest_cluster]['gene_id'].extend(groups[sid])
        #TODO: need to recalculate 
        new_size=clusters[nearest_cluster_name]['size']+1+len(groups[sid]['gene_id'])
        clusters[nearest_cluster_name]['size']=new_size
        clusters[nearest_cluster_name]['mean_length']=float((int(clusters[nearest_cluster_name]['mean_length'])*(new_size-1-len(groups[sid]['gene_id']))+int(neartest_cluster['len'])))/(new_size-len(groups[sid]['gene_id']))
        clusters[nearest_cluster_name]['max_length']=max(int(clusters[nearest_cluster_name]['max_length']),int(neartest_cluster['len']))
        clusters[nearest_cluster_name]['min_length']=min(int(clusters[nearest_cluster_name]['min_length']),int(neartest_cluster['len']))
        del groups[sid]
        added_group.add(sid)
    json.dump(clusters, open(os.path.join(out_dir, 'clusters_by_ref.json'), 'w'), indent=4, sort_keys=True)
    count_unmatched=0
    count_num_input=0
    with open(un_match_combined_faa_file, 'w') as fh, open(seqs_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            count_num_input=count_num_input+1
            if not newseq.id in added_group:
                #newseq.id=str(gene_map[newseq.id])
                newseq.description=''
                newseq.name=''
                SeqIO.write(newseq,fh,'fasta')
                count_unmatched=count_unmatched+1
    #del gene_map 
    elapsed = datetime.now() - starttime
    logging.info(f'Remain {len(groups.keys())} not matched')
    logging.info(f'Matching {count_num_input} groups to ref clusters, {len(top_match.keys())} matched, form {len(clusters.keys())} ref clusters,  {count_unmatched} groups not matched  -- time taken {str(elapsed)}')
    return un_match_combined_faa_file,clusters,groups   
# def blastq_with_ref_clusters(sample,our_dir,refdb,threads=1,evalue=1E-6, thredshold=0.7):
#     sample_id = sample['id']
#     sample_dir = os.path.join(out_dir, 'samples', sample_id)
#     seqs_file = os.path.join(sample_dir, sample_id +'.faa')
#     diamond_result = os.path.join(out_dir, 'ref_diamond.tsv')
#     cmd = f'./diamond blastp -q {seqs_file} -d {refdb} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamond_result}'
#     ret = run_command(cmd)
#     if ret != 0:
#         raise Exception('Error running diamond blastp')
#     top_match={}
#     for line in open(diamond_result, 'r'):
#         cells = line.rstrip().split('\t')            

#         sid=cells[0]
#         clustername=cells[1]
#         pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
#         if pident <= thredshold or len_diff <= thredshold or align_short <= thredshold or align_long <= thredshold:
#             continue
#         if not sid in top_match.keys():
#             top_match[sid]=[]
#         top_match[sid].append({'len':qlen,'cluster':clustername,'ident':pident} )
#     return top_match
# def clustering_sample_sequences_by_ref(sample,out_dir,refclusters, refdb,refcdb,threads=1):
#     top_match=blastq_with_ref_clusters(sample,out_dir,refdb,threads)
#     if threads == 0:
#         threads = multiprocessing.cpu_count()
#     clusters={}
#     with multiprocessing.Pool(processes=threads) as pool:
#         results = pool.map(partial(matching_sequence_to_near_clusters,clusters=clusters, refcdb=refcdb,out_dir=out_dir, top_match=top_match), top_match.keys())
# def matching_sequence_to_near_clusters(sid,clusters,refcdb,out_dir,top_match):
#     list_candidate_clusters=top_match[sid]
#     os.mkdir(temmatchingdir+"/"+sid)
        
#     for c in list_candidate_clusters:
#         diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
#         cmd = f'./diamond blastp --quiet -q {os.path.join(out_dir,'samples/'+sid+".faa")} -d {refcdb+"/"+c["cluster"]+".db.dmnd"} -p 1 --evalue 1E-6 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamondout}'
#         run_command(cmd)
#     neartest_cluster=list_candidate_clusters[0]
#     best_ident=0
#     #groups[sid]['confident']=1
#     list_matched_groups=[]
#     hash_clustername={}
#     hash_cluster_max_ident={}
#     for c in list_candidate_clusters:
#         diamondout=temmatchingdir+"/"+sid+"/"+c['cluster']+".tsv"
#         #read and note max identity group
#         hash_cluster_max_ident[c['cluster']]={'pident':0,'len_diff':0,'align_short':0,'align_long':0}
#         for line in open(diamondout, 'r'):
#             cells = line.rstrip().split('\t')            

           
#             groupname=cells[1]
#             pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
#                 #if pident <= 0.7 or len_diff <= 0.7 or align_short <= 0.7 or align_long <= 0.7:
#                 #    continue
#             list_matched_groups.append({'c':c,'g':groupname,'pid':pident})
                        
#             if pident>best_ident:
#                 best_ident=pident
#                 neartest_cluster=c
                    
#             if pident>hash_cluster_max_ident[c['cluster']]['pident']:
#                 hash_cluster_max_ident[c['cluster']]={'pident':pident,'len_diff':len_diff,'align_short':align_short,'align_long':align_long}
#         hash_clustername[c['cluster']]=c
            
#     if best_ident<0.98:
#         #if group not trully belong to a ref cluster, try to add it to best cluster by k-nearest neighbor and cal confident 
#         print("consider the undetermined :")
#         list_matched_groups.sort(key=lambda x: x['pid'], reverse=True)
#         print(list_matched_groups)
#         count_matched_cluster={}         
#         k=5
#         if len(list_matched_groups)<5:
#             k=len(list_matched_groups)
#         for i in range(k):
#             if not list_matched_groups[i]['c']['cluster'] in count_matched_cluster:
#                 count_matched_cluster[list_matched_groups[i]['c']['cluster']]=0
#             count_matched_cluster[list_matched_groups[i]['c']['cluster']]=count_matched_cluster[list_matched_groups[i]['c']['cluster']]+1
#         max_in_k_nearest=0
#         cluster_max=list_matched_groups[0]['c']['cluster']
#         for kc in count_matched_cluster.keys():                   
#             if count_matched_cluster[kc]>max_in_k_nearest:
#                 max_in_k_nearest=count_matched_cluster[kc]
#                 cluster_max=kc
#         neartest_cluster=hash_clustername[cluster_max]
#         groups[sid]['confident']=hash_cluster_max_ident[cluster_max]['pident']
#         groups[sid]['match']=hash_cluster_max_ident[cluster_max]

#      else:
#         groups[sid]['confident']=1
#         groups[sid]['match']=hash_cluster_max_ident[neartest_cluster['cluster']]
#     nearest_cluster_name=neartest_cluster['cluster']
#     if not nearest_cluster_name in clusters.keys():
            
#         clusters[nearest_cluster_name]={}
#         clusters[nearest_cluster_name]['groups']=[]
#         clusters[nearest_cluster_name]['max_length']=0
#         clusters[nearest_cluster_name]['mean_length']=0
#         clusters[nearest_cluster_name]['min_length']=1E6
#         clusters[nearest_cluster_name]['gene_name']=ref_clusters[nearest_cluster_name]['gene_name']
#         clusters[nearest_cluster_name]['product']=ref_clusters[nearest_cluster_name]['description']
#         clusters[nearest_cluster_name]['representative']='gene_families/sequences/'+nearest_cluster_name+'.fasta'
#         clusters[nearest_cluster_name]['size']=0
#         clusters[nearest_cluster_name]['source']='reference'
#         matched_cluster.add(nearest_cluster_name)
#         #groups[sid]['confident']=1
#     clusters[nearest_cluster_name]['groups'].append(groups[sid])
#         #clusters[neartest_cluster]['gene_id'].extend(groups[sid])
#         #TODO: need to recalculate 
#     new_size=clusters[nearest_cluster_name]['size']+1+len(groups[sid]['gene_id'])
#     clusters[nearest_cluster_name]['size']=new_size
#     clusters[nearest_cluster_name]['mean_length']=float((int(clusters[nearest_cluster_name]['mean_length'])*(new_size-1-len(groups[sid]['gene_id']))+int(neartest_cluster['len'])))/(new_size-len(groups[sid]['gene_id']))
#     clusters[nearest_cluster_name]['max_length']=max(int(clusters[nearest_cluster_name]['max_length']),int(neartest_cluster['len']))
#     clusters[nearest_cluster_name]['min_length']=min(int(clusters[nearest_cluster_name]['min_length']),int(neartest_cluster['len']))
def clustering(input_fasta_file,out_dir,threads,evalue=10e6,identity=0.7,LD=0.7,AS=0.7,AL=0.7,timing_log=None):
    starttime = datetime.now()
    groups_representative_fasta, similar_groups = run_mmseq_with_map_similar_seqs(
        faa_file=input_fasta_file,
        
        out_dir=out_dir,      
        threads=threads,
        timing_log=timing_log)
    blast_result = pairwise_alignment_diamond(
      
        #database_fasta = groups_representative_fasta,
        database_fasta = groups_representative_fasta,
        query_fasta = groups_representative_fasta,
        out_dir = out_dir,
        evalue = evalue,
        threads=threads,
        timing_log=timing_log)

    filtered_blast_result = filter_blast_result(
        blast_result=blast_result,
        out_dir = out_dir,
        identity=identity,
        length_difference=LD,
        alignment_coverage_short=AS,
        alignment_coverage_long=AL)

    mcl_file = cluster_with_mcl(
        out_dir = out_dir,
        blast_result = filtered_blast_result,
        threads=threads,
        timing_log=timing_log)
    inflated_clusters, groups = reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)   
    #set_of_representative_id=set()
    count_items=0
    for c in inflated_clusters:
        #print(c)
        
        #for k in c.keys():
        #    #count_items=count_items+len(c[k])
        #    set_of_representative_id.add(k)
        #    break
            #set_of_representative_id.update(c[k])
        for k in c.keys():
            count_items=count_items+len(c[k])
           
    # new_representative_clusters=os.path.join(out_dir,'new_representative_clusters.fasta')
    # with open(groups_representative_fasta,'rt') as fi, open(new_representative_clusters,'w') as fo:
    #     for r in SeqIO.parse(fi,'fasta'):
    #         if r.id in set_of_representative_id:
                
    #             SeqIO.write(r,fo,'fasta')
    elapsed = datetime.now() - starttime
    logging.info(f'Clustering with MCL with {count_items} seqs -- time taken {str(elapsed)}')
    
    return inflated_clusters,groups,groups_representative_fasta
def clustering_faiss(input_fasta_file,out_dir,threads,evalue=10e6,identity=0.7,LD=0.7,AS=0.7,AL=0.7,timing_log=None):
    starttime = datetime.now()
    groups_representative_fasta, similar_groups = run_mmseq_with_map_similar_seqs(
        faa_file=input_fasta_file,
        
        out_dir=out_dir,      
        threads=threads,
        timing_log=timing_log)
    faiss_result_file = pairwise_alignment_faiss(
      
        #database_fasta = groups_representative_fasta,
        database_fasta = groups_representative_fasta,
        query_fasta = groups_representative_fasta,
        out_dir = out_dir,
        evalue = evalue,
        threads=threads,
        timing_log=timing_log)

    

    mcl_file = cluster_with_mcl_from_faiss(
        out_dir = out_dir,
        file_faiss = faiss_result_file,
        threads=threads,
        inflation=4,
        timing_log=timing_log)
    inflated_clusters, groups = reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)   
    set_of_representative_id=set()
    count_items=0
    for c in inflated_clusters:
        #print(c)
        
        for k in c.keys():
            #count_items=count_items+len(c[k])
            set_of_representative_id.add(k)
            break
            #set_of_representative_id.update(c[k])
        for k in c.keys():
            count_items=count_items+len(c[k])
           
    new_representative_clusters=os.path.join(out_dir,'new_representative_clusters.fasta')
    with open(groups_representative_fasta,'rt') as fi, open(new_representative_clusters,'w') as fo:
        for r in SeqIO.parse(fi,'fasta'):
            if r.id in set_of_representative_id:
                
                SeqIO.write(r,fo,'fasta')
    elapsed = datetime.now() - starttime
    logging.info(f'Clustering with MCL with {count_items} seqs -- time taken {str(elapsed)}')
    
    return inflated_clusters,groups,groups_representative_fasta,new_representative_clusters
def clustering_ems(input_fasta_file,out_dir,threads,evalue=10e6,identity=0.7,LD=0.7,AS=0.7,AL=0.7,timing_log=None):
    starttime = datetime.now()
    groups_representative_fasta, similar_groups = run_mmseq_with_map_similar_seqs(
        faa_file=input_fasta_file,
        
        out_dir=out_dir,      
        threads=threads,
        timing_log=timing_log)
    ems_result_file = pairwise_alignment_ems(
      
        #database_fasta = groups_representative_fasta,
        database_fasta = groups_representative_fasta,
        query_fasta = groups_representative_fasta,
        out_dir = out_dir,
        evalue = evalue,
        threads=threads,
        timing_log=timing_log)

    

    mcl_file = cluster_with_mcl_from_faiss(
        out_dir = out_dir,
        file_faiss = ems_result_file,
        threads=threads,
        inflation=4,
        timing_log=timing_log)
    inflated_clusters, groups = reinflate_clusters(
        groups=similar_groups,
        mcl_file=mcl_file)   
    set_of_representative_id=set()
    count_items=0
    for c in inflated_clusters:
        #print(c)
        
        for k in c.keys():
            #count_items=count_items+len(c[k])
            set_of_representative_id.add(k)
            break
            #set_of_representative_id.update(c[k])
        for k in c.keys():
            count_items=count_items+len(c[k])
           
    new_representative_clusters=os.path.join(out_dir,'new_representative_clusters.fasta')
    with open(groups_representative_fasta,'rt') as fi, open(new_representative_clusters,'w') as fo:
        for r in SeqIO.parse(fi,'fasta'):
            if r.id in set_of_representative_id:
                
                SeqIO.write(r,fo,'fasta')
    elapsed = datetime.now() - starttime
    logging.info(f'Clustering with MCL with {count_items} seqs -- time taken {str(elapsed)}')
    
    return inflated_clusters,groups,groups_representative_fasta,new_representative_clusters
def make_representative_clusters_from_similar_groups(annotated_clusters_file,representative_groups,out_dir):
    annotated_clusters=json.load(open(annotated_clusters_file, 'r'))
    represent_clusters_file=os.path.join(out_dir,"representative_clusters.fasta")
    map_represent_seqid_clusterid={}
    for c in annotated_clusters.keys():
        list_geneid=read_array(annotated_clusters[c]['unique_seq'])
        for s in list_geneid:
            map_represent_seqid_clusterid[s]=c
    #json.dump(map_represent_seqid_clusterid, open(os.path.join(out_dir, 'map_represent_seqid_clusterid.json'), 'w'), indent=4, sort_keys=True)

    
    with open(represent_clusters_file,'w') as ofh, open(representative_groups) as ifh:
        for line in ifh:
            if line[0] == '>':
                gene_id=line[1:].strip()
                
                ofh.write(f'>{map_represent_seqid_clusterid[gene_id]}\n')
            else:
                ofh.write(line)  
    return  represent_clusters_file   
def make_representative_clusters(annotated_clusters_file,unique_seqs,out_dir):
    annotated_clusters= json.load(open(annotated_clusters_file, 'r'))
    represent_clusters_file=os.path.join(out_dir,"representative_clusters.fasta")
    map_represent_seqid_clusterid={}
    for c in annotated_clusters.keys():
        map_represent_seqid_clusterid[annotated_clusters[c]['representative']]=c
    #json.dump(map_represent_seqid_clusterid, open(os.path.join(out_dir, 'map_represent_seqid_clusterid.json'), 'w'), indent=4, sort_keys=True)

    
    with open(represent_clusters_file,'w') as ofh, open(unique_seqs) as ifh:
        for r in SeqIO.parse(unique_seqs,'fasta'):
            if r.id in map_represent_seqid_clusterid:
                r.id=map_represent_seqid_clusterid[r.id]
                SeqIO.write(r,ofh,'fasta')
    return  represent_clusters_file   
def identical_matching_old_clusters(old_clusters_file,  old_unique_groups_file,old_unique_sequence,new_unique_groups,new_unique_seqs,out_dir,root_dir,threads,timing_log=None):
    starttime = datetime.now()
    mem_usage = mem_report(0, "begin identical_matching_old_clusters")
    old_clusters= json.load(open(old_clusters_file, 'r'))
    
    map_seqid_to_cluster={}
    for c in old_clusters.keys():
        list_useq=read_array(old_clusters[c]['unique_seq'])
        for s in list_useq:

            map_seqid_to_cluster[s]=c
    mem_usage = mem_report(mem_usage, "identical_matching_old_clusters:map_seqid_to_cluster")
    combined_unique_seq=os.path.join(out_dir,"combined_unique_seq.faa")
    cmd=f'cat {old_unique_sequence} {new_unique_seqs} > {combined_unique_seq}'
    ret = run_command(cmd,timing_log)
    temp2=os.path.join(out_dir,"ident")
    os.mkdir(temp2)
    unique_seqs_fasta, unique_groups = run_mmseq_unique_seqs(
        faa_file=combined_unique_seq,     
        out_dir=temp2,      
        threads=threads,
        timing_log=timing_log)
    json.dump(unique_groups, open(os.path.join(temp2, 'groups_unique_seqss.json'), 'w'), indent=4, sort_keys=True)
    #json.dump(new_unique_groups, open(os.path.join(temp2, 'new_unique_groups1.5.json'), 'w'), indent=4, sort_keys=True)
    mem_usage = mem_report(mem_usage, "identical_matching_old_clusters:unique_groups")
    additon_cluster_geneid={}
    additon_group_geneid={}
    count_ident=0
    count_ident_group=0
    for g  in unique_groups:
        list_ids=[]
        list_ids.append(g)
        for s in unique_groups[g]:
            list_ids.append(s)
        matched=False
        cluster_matched=None
        id_matched=None
        for id in list_ids:
            if id in map_seqid_to_cluster.keys():
                cluster_matched=map_seqid_to_cluster[id]
                id_matched=id
                matched=True
                break
        #print(list_ids)
        if matched:
            for id in list_ids:
                if not id in map_seqid_to_cluster.keys() and id in new_unique_groups.keys() :
                    if cluster_matched not in additon_cluster_geneid:
                        additon_cluster_geneid[cluster_matched]=[]
                    additon_cluster_geneid[cluster_matched].append(id)
                    additon_cluster_geneid[cluster_matched].extend(new_unique_groups[id])
                    
                    #old_clusters[cluster_matched]['gene_id'].append(id)
                    #old_clusters[cluster_matched]['gene_id'].extend(new_unique_groups[id])
                    #old_unique_groups[id_matched].append(id)
                    #old_unique_groups[id_matched].extend(new_unique_groups[id])
                    if id_matched not in additon_group_geneid:
                        additon_group_geneid[id_matched]=[]
                    additon_group_geneid[id_matched].append(id)
                    additon_group_geneid[id_matched].extend(new_unique_groups[id])
                    count_ident=count_ident+len(new_unique_groups[id])+1
                    count_ident_group=count_ident_group+1
                    
                    del new_unique_groups[id]
    del map_seqid_to_cluster
    
    mem_usage = mem_report(mem_usage, "identical_matching_old_clusters:additon_group_geneid")    
    new_reduced_unique_seqs=os.path.join(out_dir,'reduced_new_unique_seqs.fasta')
    with open(new_unique_seqs) as fi, open(new_reduced_unique_seqs,'w')as fo:
        for r in SeqIO.parse(fi,'fasta'):
            if r.id in new_unique_groups.keys():
                SeqIO.write(r,fo,'fasta')
    for c in additon_cluster_geneid.keys():
        list_geneid=read_array(old_clusters[c]['gene_id'])
        list_geneid.extend(additon_cluster_geneid[c])
        write_array(old_clusters[c]['gene_id'],list_geneid)
    #json.dump(old_clusters, open(old_clusters_file, 'w'), indent=4, sort_keys=True)
    del old_clusters
    old_unique_groups=json.load(open(old_unique_groups_file, 'r'))
    for g in additon_group_geneid.keys():
        list_useq=read_array(old_unique_groups[g])
        list_useq.extend(additon_group_geneid[g])
        write_array(old_unique_groups[g],list_useq)
    for k in new_unique_groups.keys():
        old_unique_groups[k]=add_unique_groups(k,new_unique_groups[k],root_dir)
        
    json.dump(old_unique_groups, open(old_unique_groups_file, 'w'), indent=4, sort_keys=True)
    cmd = f'cat {new_reduced_unique_seqs} >> {old_unique_sequence}'
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error concat '+cmd)
    del old_unique_groups
    del additon_cluster_geneid
    del additon_group_geneid
    mem_usage = mem_report(mem_usage, "identical_matching_old_clusters:del old_unique_groups") 
    elapsed = datetime.now() - starttime   
    logging.info(f'identical matching {count_ident} sequences in {count_ident_group} groups, remain {len(new_unique_groups.keys())}   -- time taken {str(elapsed)}')
    return old_clusters_file,old_unique_groups_file,new_unique_groups,new_reduced_unique_seqs
def simple_match_new_seq_to_old_clusters(old_clusters_file,representative_clusters_file, new_sequences_file,old_sequences_file,out_dir, identity=0.7,evalue=0.000001, threads=1,timing_log=None):
    starttime = datetime.now()
    mem_usage = mem_report(0, "begin simple_match_new_seq_to_old_clusters")
    dict_gene_lenght={}
    with open(new_sequences_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            dict_gene_lenght[newseq.id]=len(newseq.seq)
    logging.info(f'number of new sequence before matching: {str(len(dict_gene_lenght.keys()))}')  
    
    diamond_result = os.path.join(out_dir, 'matching_new_sequences_diamond.tsv')
    cmd = f'./diamond blastp -q {new_sequences_file} -d {representative_clusters_file} -p {threads} --evalue {evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 1 2> /dev/null 1> {diamond_result}'
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error running diamond blastp')
    elapsed = datetime.now() - starttime
   
    logging.info(f'diamond blastp  -- time taken {str(elapsed)}')
    top_match={}
    mem_usage = mem_report(mem_usage, "simple_match_new_seq_to_old_clusters:diamond")
    for line in open(diamond_result, 'r'):
        cells = line.rstrip().split('\t')            

        sid=cells[0]
        clustername=cells[1]
        pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
        if pident <= identity or len_diff <= identity or align_short <= identity or align_long <= identity:
            continue
        if not sid in top_match.keys():
            top_match[sid]={'len':qlen,'cluster':clustername,'ident':pident}
        if pident > top_match[sid]['ident']:
            top_match[sid]['ident']=pident
    elapsed = datetime.now() - starttime
   
    logging.info(f'open diamond_result  -- time taken {str(elapsed)}')
    old_clusters= json.load(open(old_clusters_file, 'r'))
    update_clusters={}
    for t in top_match.keys():
        nearest_cluster_name=top_match[t]['cluster']
        if not nearest_cluster_name in update_clusters:
            update_clusters[nearest_cluster_name]=set()
        #list_geneid=read_array(old_clusters[nearest_cluster_name]['gene_id'])
        #list_geneid.append(t)
        update_clusters[nearest_cluster_name].add(t)
        #old_clusters[nearest_cluster_name]['gene_id']=write_array(old_clusters[nearest_cluster_name]['gene_id'],list_geneid)
        #list_useq=read_array(old_clusters[nearest_cluster_name]['unique_seq'])
        #list_useq.append(t)
        #old_clusters[nearest_cluster_name]['unique_seq']=write_array(old_clusters[nearest_cluster_name]['unique_seq'],list_useq)
        #clusters[neartest_cluster]['gene_id'].extend(groups[sid])
        #TODO: need to recalculate 
        # new_size=old_clusters[nearest_cluster_name]['size']+1
        # old_clusters[nearest_cluster_name]['size']=new_size
        # old_clusters[nearest_cluster_name]['mean_length']=float((int(old_clusters[nearest_cluster_name]['mean_length'])*(new_size-1)+int(top_match[t]['len'])))/(new_size-1)
        # old_clusters[nearest_cluster_name]['max_length']=max(int(old_clusters[nearest_cluster_name]['max_length']),int(top_match[t]['len']))
        # old_clusters[nearest_cluster_name]['min_length']=min(int(old_clusters[nearest_cluster_name]['min_length']),int(top_match[t]['len']))
        # old_clusters[nearest_cluster_name]['updated']=1
    #json.dump(old_clusters, open(old_clusters_file, 'w'), indent=4, sort_keys=True)
    count_seq_in_update_cluster=0
    for c in update_clusters.keys():
        if len(update_clusters[c])==0:
            continue
        count_seq_in_update_cluster=count_seq_in_update_cluster+len(update_clusters[c])
        list_geneid=read_array(old_clusters[c]['gene_id'])
        list_useq=read_array(old_clusters[c]['unique_seq'])
        
        #clusters[neartest_cluster]['gene_id'].extend(groups[sid])
        #TODO: need to recalculate 
        #print(update_clusters[c])
        mean_lenght=0
        min_lenght=0
        max_lenght=9999999
        for g in update_clusters[c]:           
            l=dict_gene_lenght[g]
            if l>max_lenght:
                max_lenght=l
            if l<min_lenght:
                min_lenght=l
            mean_lenght=mean_lenght+l
        mean_lenght=mean_lenght/len(update_clusters[c])
        new_size=old_clusters[c]['size']+len(update_clusters[c])
        #print(old_clusters[c])
        old_clusters[c]['size']=new_size
        old_clusters[c]['mean_length']=(old_clusters[c]['mean_length']*len(list_geneid)+mean_lenght*len(update_clusters[c]))/new_size
        #old_clusters[c]['mean_length']=float((int(old_clusters[c]['mean_length'])*(new_size-1)+int(top_match[t]['len'])))/(new_size-1)
        old_clusters[c]['max_length']=max(int(old_clusters[c]['max_length']),max_lenght)
        old_clusters[c]['min_length']=min(int(old_clusters[c]['min_length']),min_lenght)
        old_clusters[c]['updated']=1
        list_geneid.extend(update_clusters[c])
        old_clusters[c]['gene_id']=write_array(old_clusters[c]['gene_id'],list_geneid)
        list_useq.extend(update_clusters[c])
        old_clusters[c]['unique_seq']=write_array(old_clusters[c]['unique_seq'],list_useq)
    elapsed = datetime.now() - starttime   
    logging.info(f'matching  -- time taken {str(elapsed)}')
    mem_usage = mem_report(mem_usage, "simple_match_new_seq_to_old_clusters:old_clusters")
    json.dump(old_clusters, open(old_clusters_file, 'w'), indent=4, sort_keys=True)
   
    del old_clusters
    un_match_combined_faa_file = os.path.join(out_dir, 'unmatch_combined.faa')
    count_num_input=0
    count_unmatched=0
    with open(un_match_combined_faa_file, 'w') as fh, open(new_sequences_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            count_num_input=count_num_input+1
            if not newseq.id in top_match:
                #newseq.id=str(gene_map[newseq.id])
                newseq.description=''
                newseq.name=''
                SeqIO.write(newseq,fh,'fasta')
                count_unmatched=count_unmatched+1
    elapsed = datetime.now() - starttime   
    logging.info(f'un_match_combined_faa_file  -- time taken {str(elapsed)}')
    #del gene_map
    # cmd = f'cat {un_match_combined_faa_file} >> {old_sequences_file}'
    # ret = run_command(cmd,timing_log)
    # if ret != 0:
    #     raise Exception('Error concat '+cmd)
    mem_usage = mem_report(mem_usage, "simple_match_new_seq_to_old_clusters:del old_clusters")
    elapsed = datetime.now() - starttime
   # logging.info(f'Remain {len(groups.keys())} not matched')
    logging.info(f'Matching {count_num_input} groups to ref clusters, {len(top_match.keys())} matched,  {count_unmatched} groups not matched  -- time taken {str(elapsed)}')
    return un_match_combined_faa_file,old_clusters_file   
def hc_match_new_seq_to_old_clusters(old_clusters_file,representative_clusters_file, new_sequences_file,old_sequences_file,out_dir, identity=0.7,evalue=0.000001, threads=1,timing_log=None):
    starttime = datetime.now()
    mem_usage = mem_report(0, "begin hc matching")
    combined_new_and_old_seqs_file=os.path.join(out_dir,'combined_new_and_old_seqs_file.fasta')
    cmd = f'cat {new_sequences_file} {representative_clusters_file} >> {combined_new_and_old_seqs_file}'
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error concat '+cmd)
    old_clusters= json.load(open(old_clusters_file, 'r'))
    dict_gene_lenght={}
    with open(new_sequences_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            dict_gene_lenght[newseq.id]=len(newseq.seq)
    logging.info(f'number of new sequence before matching: {str(len(dict_gene_lenght.keys()))}')  
    update_clusters = {}
    dict_added_genes={}
    #interative matching:
    #for t in ['0.98','0.95','0.90','0.85','0.80','0.75','0.70']:
    #for t in ['0.98','0.95','0.90','0.85']:
    for i in range(98,70,-1):
        t=float(i)/100
        mmseq_cluster_file=os.path.join(out_dir, 'mmseq_groups')
        cmd = f'mmseqs easy-linclust {combined_new_and_old_seqs_file} {mmseq_cluster_file} {out_dir}/tmp --min-seq-id {t} -c {t} --cov-mode 0 --threads {threads} > /dev/null'    
        ret = run_command(cmd)
        #handle output and mapping seq to cluster
        if ret != 0:
            raise Exception('Error mmseqs '+cmd)
        cluster_count=0
        mmseq_cluster_file=mmseq_cluster_file+"_cluster.tsv"
        list_genes=set()
          
        with open(mmseq_cluster_file, 'r') as fh:
            for line in fh:
                rep_name,member = line.strip().split()
                rep_name=rep_name.strip()
                member=member.strip()

                if rep_name == member:
                    if len(list_genes)>0:
                        #handle previous list_genes
                        isFound=False
                        clusters_found=[]

                        for g in list_genes:
                            if g in old_clusters.keys():
                                isFound=True
                                clusters_found.append(g)
                                
                        if isFound and len(clusters_found)>0:
                            
                            if not clusters_found[0] in update_clusters.keys():
                                update_clusters[clusters_found[0]]=set()
                            for c in clusters_found:
                                list_genes.discard(c)
                            for g in list_genes:
                                if g in dict_added_genes:
                                    logging.info(f'gene {g} already added: {dict_added_genes[g]} at {t}')    
                                dict_added_genes[g]=clusters_found[0]
                                update_clusters[clusters_found[0]].add(g)
                        
                        list_genes.clear()            
                list_genes.add(member)                
            if len(list_genes)>0:
                #handle previous list_genes
                isFound=False
                clusters_found=[]

                for g in list_genes:
                    if g in old_clusters.keys():
                        isFound=True
                        clusters_found.append(g)
                                
                if isFound and len(clusters_found)>0:
                            
                    if not clusters_found[0] in update_clusters.keys():
                        update_clusters[clusters_found[0]]=set()
                    for c in clusters_found:
                        list_genes.discard(c)
                    for g in list_genes:
                        dict_added_genes[g]=clusters_found[0]
                        update_clusters[clusters_found[0]].add(g)
                list_genes.clear()
        logging.info(f'gene added: {str(len(dict_added_genes.keys()))}')    
        #create new sequence
        if len(dict_added_genes.keys())>0:
            remain_sequences_file=os.path.join(out_dir, 'remain_seqs.fasta')
            count_remain=0
            with open(new_sequences_file, 'r') as fh, open(remain_sequences_file,'w') as fo:
                for newseq in SeqIO.parse(fh,'fasta'):
                    if not newseq.id in dict_added_genes.keys():
                        count_remain=count_remain+1
                        newseq.description=''
                        newseq.name=''
                        SeqIO.write(newseq,fo,'fasta')
            cmd = f'cat {remain_sequences_file} {representative_clusters_file} > {combined_new_and_old_seqs_file}'
            ret = run_command(cmd,timing_log)
            logging.info(f'remain seq {str(count_remain)}')
    
    #compare update_clusters and dictionary
    list_dup={}
    #update old clusters
    count_seq_in_update_cluster=0
    for c in update_clusters.keys():
        if len(update_clusters[c])==0:
            continue
        count_seq_in_update_cluster=count_seq_in_update_cluster+len(update_clusters[c])
        list_geneid=read_array(old_clusters[c]['gene_id'])
        list_useq=read_array(old_clusters[c]['unique_seq'])
        
        #clusters[neartest_cluster]['gene_id'].extend(groups[sid])
        #TODO: need to recalculate 
        #print(update_clusters[c])
        mean_lenght=0
        min_lenght=0
        max_lenght=9999999
        for g in update_clusters[c]:
            if g in list_dup:
                list_dup[g]=list_dup[g]+1
            else:
                list_dup[g]=1

            l=dict_gene_lenght[g]
            if l>max_lenght:
                max_lenght=l
            if l<min_lenght:
                min_lenght=l
            mean_lenght=mean_lenght+l
        mean_lenght=mean_lenght/len(update_clusters[c])
        new_size=old_clusters[c]['size']+len(update_clusters[c])
        #print(old_clusters[c])
        old_clusters[c]['size']=new_size
        old_clusters[c]['mean_length']=(old_clusters[c]['mean_length']*len(list_geneid)+mean_lenght*len(update_clusters[c]))/new_size
        #old_clusters[c]['mean_length']=float((int(old_clusters[c]['mean_length'])*(new_size-1)+int(top_match[t]['len'])))/(new_size-1)
        old_clusters[c]['max_length']=max(int(old_clusters[c]['max_length']),max_lenght)
        old_clusters[c]['min_length']=min(int(old_clusters[c]['min_length']),min_lenght)
        old_clusters[c]['updated']=1
        list_geneid.extend(update_clusters[c])
        old_clusters[c]['gene_id']=write_array(old_clusters[c]['gene_id'],list_geneid)
        list_useq.extend(update_clusters[c])
        old_clusters[c]['unique_seq']=write_array(old_clusters[c]['unique_seq'],list_useq)
    for s in list_dup:
        if list_dup[s]>1:
            logging.info(f'duplicate in updated_clusters  {s} : {list_dup[s]} times in {dict_added_genes[s]}')
    un_match_combined_faa_file = os.path.join(out_dir, 'unmatch_combined.faa')
    count_num_input=0
    count_unmatched=0
    with open(un_match_combined_faa_file, 'w') as fh, open(new_sequences_file, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            count_num_input=count_num_input+1
            if not newseq.id in dict_added_genes.keys():
                #newseq.id=str(gene_map[newseq.id])
                newseq.description=''
                newseq.name=''
                SeqIO.write(newseq,fh,'fasta')
                count_unmatched=count_unmatched+1
    elapsed = datetime.now() - starttime   
    logging.info(f'un_match_combined_faa_file  -- time taken {str(elapsed)}')
    #del gene_map
    cmd = f'cat {un_match_combined_faa_file} >> {old_sequences_file}'
    ret = run_command(cmd,timing_log)
    if ret != 0:
        raise Exception('Error concat '+cmd)
    mem_usage = mem_report(mem_usage, "hc_match_new_seq_to_old_clusters:del old_clusters")
    json.dump(old_clusters, open(old_clusters_file, 'w'), indent=4, sort_keys=True)
    
    elapsed = datetime.now() - starttime
   # logging.info(f'Remain {len(groups.keys())} not matched')
    logging.info(f'Matching {count_num_input} groups to {len(old_clusters.keys())} old clusters, {len(dict_added_genes.keys())} matched {count_seq_in_update_cluster},  {count_unmatched} groups not matched  -- time taken {str(elapsed)}')
    return un_match_combined_faa_file,old_clusters_file   
def refine_updated_clusters(updated_clusters_file,new_annotated_clusters,old_representatve,new_representative,gene_annotation_fn,unique_sequences,out_dir,root_dir,threads):
    #collect all representative seq of all clusters
    updated_clusters= json.load(open(updated_clusters_file, 'r'))
    starttime = datetime.now()
    mem_usage = mem_report(0, "refine_updated_clusters begin")
    map_seqid_to_cluster={}
    
    for c in new_annotated_clusters.keys():
        new_annotated_clusters[c]['updated']=1
        
        list_useq=read_array(new_annotated_clusters[c]['unique_seq'])
        for s in list_useq:

            map_seqid_to_cluster[s]=c
    for c in updated_clusters.keys():
        rep=updated_clusters[c]['representative']
        if rep in map_seqid_to_cluster.keys():
            map_seqid_to_cluster[rep]=c
    del updated_clusters
    json.dump(new_annotated_clusters, open(os.path.join(out_dir,"new_clusters"), 'w'), indent=4, sort_keys=True)
    #new_representative=os.path.join(out_dir,"combined_representative.fasta")
    logging.info("size of map_seqid_to_cluster: "+str(sys.getsizeof(map_seqid_to_cluster)))
    mem_usage = mem_report(mem_usage, "refine_updated_clusters:map_seqid_to_cluster")
    with open(old_representatve, 'a') as fh, open(new_representative, 'rt') as fi:
        for newseq in SeqIO.parse(fi,'fasta'):
            newseq.description=''
            newseq.name=''
            newseq.id=map_seqid_to_cluster[newseq.id]
            SeqIO.write(newseq,fh,'fasta')    
    #....
    list_groups_clusters_need_refine_file,annotated_clusters_file=split_clusters_by_connected_component(
        annotated_clusters_file=updated_clusters_file,
        representative_clusters_file=old_representatve,
        out_dir=out_dir,
        threads=threads
    )
    mem_usage = mem_report(mem_usage, "refine_updated_clusters:split_clusters_by_connected_component")
    #print(list_groups_clusters_need_refine)
    list_groups_clusters_need_refine= json.load(open(list_groups_clusters_need_refine_file, 'r'))
    new_list_clusters={}
    new_inflate_clusters=[]
    for groups in list_groups_clusters_need_refine:
        new_refined_inflated_clusters=reclustering(groups,unique_sequences,gene_annotation_fn,out_dir=os.path.join(out_dir,"temp_refine"),threads=threads)
        if new_refined_inflated_clusters != None:
            new_inflate_clusters.extend(new_refined_inflated_clusters)
    if len(new_inflate_clusters)>0:
        new_annotated_clusters = annotate_cluster(
            unlabeled_clusters=new_inflate_clusters,
            gene_annotation_fn=gene_annotation_fn)

        annotated_clusters_file=merge_new_cluster_to_old_clusters(new_annotated_clusters,annotated_clusters_file)
    del list_groups_clusters_need_refine
    mem_usage = mem_report(mem_usage, "refine_updated_clusters:merge_new_cluster_to_old_clusters")
    check_clusters(annotated_clusters_file,root_dir)
    elapsed = datetime.now() - starttime
    logging.info(f'Done refine clusters-- time taken {str(elapsed)}')
    return annotated_clusters_file
import networkx as nx
def split_clusters_by_connected_component(annotated_clusters_file,representative_clusters_file,out_dir,threads,identity=0.7):   
    diamond_result = os.path.join(out_dir, 'pairwise_rep_clusters_diamond.tsv')
    cmd = f'./diamond blastp -q {representative_clusters_file} -d {representative_clusters_file} -p {threads} --evalue 1e-6 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --sensitive --max-target-seqs 2000 2> /dev/null 1> {diamond_result}'
    ret = run_command(cmd)
    if ret != 0:
        raise Exception('Error running diamond blastp')
    annotated_clusters= json.load(open(annotated_clusters_file, 'r'))
    G = nx.Graph()
    for c in annotated_clusters.keys():
        G.add_node(c)
    for line in open(diamond_result, 'r'):
        cells = line.rstrip().split('\t')            

        sid1=cells[0]
        sid2=cells[1]
        pident,len_diff,align_short,align_long,qlen=getIdentAlignFromCell(cells)
        if pident <= identity or len_diff <= identity or align_short <= identity or align_long <= identity:
            continue
        G.add_edge(sid1, sid2)
    nx.write_gml(G,os.path.join(out_dir,"clusters_graph.gml"))
    connected_component=nx.connected_components(G)
   #print('connected_component')
    
    list_groups_clusters_need_refine=[]
    for com in connected_component:
        #print('com')
        if len(com)<=1:
            continue
        need_refine=False
        for c in com:
            #print(c)
            #print( annotated_clusters[c])
            if 'updated' in annotated_clusters[c] and  annotated_clusters[c]['updated']==1:
                need_refine=True
                break
        groups_clusters_need_refine=[]
        if need_refine:
            logging.info("size of annotated_clusters before del: "+str(len(annotated_clusters.keys())))
            for c in com:
                groups_clusters_need_refine.append(annotated_clusters[c])
                del annotated_clusters[c]
            logging.info("size of annotated_clusters after del: "+str(len(annotated_clusters.keys())))
        logging.info("need_refiend="+str(need_refine)+",size of groups_clusters_need_refine : "+str(len(groups_clusters_need_refine)))
        if len(groups_clusters_need_refine)>1:
            list_groups_clusters_need_refine.append(groups_clusters_need_refine)
    json.dump(annotated_clusters, open(annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    groups_clusters_need_refine_file=os.path.join(out_dir, 'groups_clusters_need_refine.json')
    json.dump(list_groups_clusters_need_refine, open(groups_clusters_need_refine_file, 'w'), indent=4, sort_keys=True)
    del annotated_clusters
    del list_groups_clusters_need_refine
    return groups_clusters_need_refine_file,annotated_clusters_file
def merge_new_cluster_to_old_clusters(new_annotated_clusters,old_annotated_clusters_file):
    starttime = datetime.now()
    old_annotated_clusters=json.load(open(old_annotated_clusters_file, 'r'))
    num_old_clusters=len(old_annotated_clusters.keys())
    count=0
    for clustername in new_annotated_clusters:
        if clustername in old_annotated_clusters:
            count=count+1
            #search max suffix

            suffix=1
            while clustername+'_'+str(suffix) in old_annotated_clusters:
                suffix=suffix+1
            
            old_annotated_clusters[clustername+'_'+str(suffix)]=new_annotated_clusters[clustername]
        else:
            old_annotated_clusters[clustername]=new_annotated_clusters[clustername]
        
    json.dump(old_annotated_clusters, open(old_annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    elapsed = datetime.now() - starttime
    logger.info(f'Merge {len(new_annotated_clusters.keys())} unmatched cluster to {num_old_clusters} old clusters, {count} duplicate -- time taken {str(elapsed)}')
    return old_annotated_clusters_file    
def reclustering(group_clusters,unique_seqs_file,gene_annotation_fn,out_dir,threads):
    starttime = datetime.now()
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    unique_seq_list=set()
   
    for c in group_clusters:
        list_useq=read_array(c['unique_seq'])
        unique_seq_list.update(list_useq)
    #print('unique_seqs_file='+unique_seqs_file)
    #logger.info(unique_seq_list)
    logger.info(f'len unique_seq_list before reclustering {str(len(unique_seq_list))}')
    if len(unique_seq_list)==0:
        #print(group_clusters)
        return None
    file_input=os.path.join(out_dir,'temp_group_clusters_unique_seqs.fasta')
    with open(unique_seqs_file) as fi, open(file_input,'w') as fo:
        for r in SeqIO.parse(fi,'fasta'):
            if r.id in unique_seq_list:
                SeqIO.write(r,fo,'fasta')
    new_inflated_clusters, similar_groups, representative_similar_groups,new_representative_clusterss=clustering(
        input_fasta_file=file_input,
        out_dir=out_dir,
        threads=threads
        )
    # new_annotated_clusters = annotate_cluster(
    #     unlabeled_clusters=new_inflated_clusters,
    #     gene_annotation_fn=gene_annotation_fn)
    
    logger.info(new_inflated_clusters)
    count_after_reflate=0
    for c in new_inflated_clusters:
        count_after_reflate=count_after_reflate+len(c)
    elapsed = datetime.now() - starttime
    logger.info(f'reclustering {len(group_clusters)} clusters with {len(unique_seq_list)} unique seq, get {count_after_reflate} seqs back -- time taken {str(elapsed)}')
    
    return new_inflated_clusters
def merge_new_cluster_to_old_clusters_by_core_analysis(new_annotated_clusters,old_annotated_clusters_file,old_unique_seqs,new_unique_seqs,gene_present_tab,root_dir,temp_dir,threads, timing_log=None):
    starttime = datetime.now()
    old_annotated_clusters= json.load(open(old_annotated_clusters_file, 'r'))
    num_old_clusters=len(old_annotated_clusters.keys())
    count=0
    cluster_dir=os.path.join(root_dir,'clusters')
    map_seq_new_cluster={}
    map_count_seq_newcluster={}
    for clustername in new_annotated_clusters:
        for us in new_annotated_clusters[clustername]['unique_seq']:
            map_seq_new_cluster[us]=clustername
        map_count_seq_newcluster[clustername]=len(new_annotated_clusters[clustername]['unique_seq'])
    map_seq_old_cluster={}
    sample_size=100
    
    for clustername in old_annotated_clusters:
        list_us=read_array(old_annotated_clusters[clustername]['unique_seq'])
        if len(list_us)>100:
            list_us = random.choices(my_list, k=sample_size)
        for us in list_us:
            map_seq_old_cluster[us]=clustername 
        
    #split core and shell clusters:
    df = pd.read_csv(gene_present_tab, sep="\t")  
    row_sums = df.iloc[:, 1:].sum(axis=1)
    k = df.shape[1] - 1  # Total number of binary columns
    row_means = row_sums / k
    filtered_df = df[row_means >= 0.15]
    core_set = set(filtered_df['Gene']) 
    #make new old unique seq:
    print('core set:'+str(len(core_set)))
    count_cloud_cluster=0
    old_unique_core_seqs=os.path.join(temp_dir,'old_unique_core_seqs.fasta')
    old_unique_cloud_seqs=os.path.join(temp_dir,'old_unique_cloud_seqs.fasta')
    with open(old_unique_seqs,'r') as fi, open(old_unique_core_seqs,'w') as fc, open(old_unique_cloud_seqs,'w') as fh:
        for seq in SeqIO.parse(fi,'fasta'):
            if seq.id in map_seq_old_cluster:
                if map_seq_old_cluster[seq.id] in core_set:
                    SeqIO.write(seq,fc,'fasta') 
                else:
                    count_cloud_cluster=count_cloud_clucster+1
                    SeqIO.write(seq,fh,'fasta')

    combined_unique_core_seqs=os.path.join(temp_dir,'combined_unique_core_seqs.fasta')
    cmd=f'cat {old_unique_core_seqs} {new_unique_seqs} > {combined_unique_core_seqs}'
    ret = run_command(cmd,timing_log) 
    groups_representative_fasta, similar_groups = run_mmseq_with_map_similar_seqs(
        faa_file=combined_unique_core_seqs,    
        out_dir=temp_dir,      
        threads=threads,
        timing_log=timing_log,
        identity=0.95)
    pair_clusters={}
    for g in similar_groups:
        if len(similar_groups[g])>1:
            items=[g]
            items.extend(similar_groups[g])
            items.sort()
            for i in range(len(items)-1):           
                for j in range(i+1, len(items)):
                    oc=None
                    if items[i] in map_seq_old_cluster:
                        oc=map_seq_old_cluster[items[i]]
                    if items[j] in map_seq_old_cluster:
                        oc=map_seq_old_cluster[items[j]]
                    nc=None
                    if items[i] in map_seq_new_cluster:
                        nc=map_seq_new_cluster[items[i]]
                    if items[j] in map_seq_new_cluster:
                        nc=map_seq_new_cluster[items[j]]
                    if oc==None or nc==None:
                        continue
                    else:
                        if oc not in pair_clusters:
                            pair_clusters[oc]={}
                        if nc not in pair_clusters[oc]:
                            pair_clusters[oc][nc]=0
                        pair_clusters[oc][nc]=pair_clusters[oc][nc]+1
    #filter
    for oc in pair_clusters:
        list_del=[]
        for nc in  pair_clusters[oc]:
            if pair_clusters[oc][nc]<map_count_seq_newcluster[nc]*0.5:
                list_del.append(nc)
        for nc in list_del:
            del pair_clusters[oc][nc]
    #merge core cluster:
    """ annotated_clusters[cluster_new_name] = {
            'gene_id':gene_id_list,
            'unique_seq':unique_seq,
            'product':cluster_product,
            'representative': unique_seq[0],
            'min_length': lens[0],
            'max_length': lens[1],
            'mean_length': lens[2],
            'size': lens[3],
            'source':'clustering'
            } """
    set_del_new_clusters=set()
    for oc in pair_clusters:
        for nc in pair_clusters[oc]:
            old_annotated_clusters[oc]=merge_2_annotated_clusters(old_annotated_clusters[oc],new_annotated_clusters[nc])
            set_del_new_clusters.add(nc)
    logger.info(f'Merge {len(set_del_new_clusters)} new clusters to old {len(pair_clusters.keys())} cluster  ')
    
    #merge cloud cluster:
    new_unique_cloud_seqs=os.path.join(temp_dir,'new_unique_cloud_seqs.fasta')
    with open(new_unique_seqs,'r') as fi, open(new_unique_cloud_seqs,'w') as fc:
        for seq in SeqIO.parse(fi,'fasta'):
            if seq.id in map_seq_new_cluster:
                if not map_seq_new_cluster[seq.id] in set_del_new_clusters:
                    SeqIO.write(seq,fc,'fasta') 
    combined_unique_cloud_seqs=os.path.join(temp_dir,'combined_unique_cloud_seqs.fasta')
    cmd=f'cat {old_unique_cloud_seqs} {new_unique_cloud_seqs} > {combined_unique_cloud_seqs}'
    ret = run_command(cmd,timing_log)
    groups_representative_fasta, similar_groups = run_mmseq_with_map_similar_seqs(
        faa_file=combined_unique_cloud_seqs,    
        out_dir=temp_dir,      
        threads=threads,
        timing_log=timing_log,
        identity=0.90)
    pair_clusters={}
    for g in similar_groups:
        if len(similar_groups[g])>1:
            items=[g]
            items.extend(similar_groups[g])
            items.sort()
            for i in range(len(items)-1):           
                for j in range(i+1, len(items)):
                    oc=None
                    if items[i] in map_seq_old_cluster:
                        oc=map_seq_old_cluster[items[i]]
                    if items[j] in map_seq_old_cluster:
                        oc=map_seq_old_cluster[items[j]]
                    nc=None
                    if items[i] in map_seq_new_cluster:
                        nc=map_seq_new_cluster[items[i]]
                    if items[j] in map_seq_new_cluster:
                        nc=map_seq_new_cluster[items[j]]
                    if oc==None or nc==None:
                        continue
                    else:
                        if oc not in pair_clusters:
                            pair_clusters[oc]={}
                        if nc not in pair_clusters[oc]:
                            pair_clusters[oc][nc]=0
                        pair_clusters[oc][nc]=pair_clusters[oc][nc]+1
    for oc in pair_clusters:
        for nc in pair_clusters[oc]:
            old_annotated_clusters[oc]=merge_2_annotated_clusters[old_annotated_clusters[oc],new_annotated_clusters[nc]]
            set_del_new_clusters.add(nc)
    for nc in set_del_new_clusters:
        del new_annotated_clusters[nc]
    
    json.dump(old_annotated_clusters, open(old_annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    old_annotated_clusters_file=merge_new_cluster_to_old_clusters(new_annotated_clusters,old_annotated_clusters_file)
    #json.dump(old_annotated_clusters, open(old_annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    save_update_clusters(old_annotated_clusters_file,root_dir)
    del old_annotated_clusters
    elapsed = datetime.now() - starttime
    logger.info(f'Merge done, remain {len(new_annotated_clusters.keys())} clusters  -- time taken {str(elapsed)}')
    return old_annotated_clusters_file,new_annotated_clusters
def merge_new_cluster_to_old_clusters_by_core_analysis2(new_annotated_clusters,old_annotated_clusters_file,old_unique_seqs,new_unique_seqs,gene_present_tab,root_dir,temp_dir,threads, timing_log=None):
    starttime = datetime.now()
    old_annotated_clusters= json.load(open(old_annotated_clusters_file, 'r'))
    map_rep_seq_old_cluster={}
    
    
    for clustername in old_annotated_clusters.keys():
       
        map_rep_seq_old_cluster[old_annotated_clusters[clustername]['representative']]=clustername
    
    representative_core_old_cluster_seqs=os.path.join(temp_dir,"representative_old_core_cluster.fasta")
    representative_cloud_old_cluster_seqs=os.path.join(temp_dir,"representative_old_cloud_cluster.fasta")
    df = pd.read_csv(gene_present_tab, sep="\t")  
    row_sums = df.iloc[:, 1:].sum(axis=1)
    k = df.shape[1] - 1  # Total number of binary columns
    row_means = row_sums / k
    filtered_df = df[row_means >= 0.15]
    core_set = set(filtered_df['Gene']) 
    logger.info(f'core_set:{len(core_set)}')
    with open(old_unique_seqs,'r') as fi, open(representative_core_old_cluster_seqs,'w') as fc, open(representative_cloud_old_cluster_seqs,'w') as fh:
        for seq in SeqIO.parse(fi,'fasta'):
            if seq.id in map_rep_seq_old_cluster:
                if map_rep_seq_old_cluster[seq.id] in core_set:
                    SeqIO.write(seq,fc,'fasta') 
                else:                   
                    SeqIO.write(seq,fh,'fasta')
    map_rep_seq_new_cluster={}
    for clustername in new_annotated_clusters:
       
        map_rep_seq_new_cluster[new_annotated_clusters[clustername]['representative']]=clustername
    representative_new_cluster_seqs=os.path.join(temp_dir,"representative_new_clusters.fasta")
    with open(new_unique_seqs,'r') as fi, open(representative_new_cluster_seqs,'w') as fc:
        for seq in SeqIO.parse(fi,'fasta'):
            if seq.id in map_rep_seq_new_cluster:
                SeqIO.write(seq,fc,'fasta') 
    #pairwise diamond with rep of old_cluster with rep of new clusters
    blast_result = pairwise_alignment_diamond(
      
        #database_fasta = groups_representative_fasta,
        database_fasta = representative_core_old_cluster_seqs,
        query_fasta = representative_new_cluster_seqs,
        out_dir = temp_dir,
        evalue = 1E-6,
        threads=threads,
        timing_log=timing_log)

    filtered_blast_result_file = filter_blast_result(
        blast_result=blast_result,
        out_dir = temp_dir,
        identity=0.7,
        length_difference=0.7,
        alignment_coverage_short=0.7,
        alignment_coverage_long=0.7)
    max_match_new_cluster_old_cluster={}
    for line in open(blast_result, 'r'):
        
        cells = line.rstrip().split('\t')            
        new_rep_id=cells[0]
        old_rep_id=cells[1]
        old_cluster_name=map_rep_seq_old_cluster[old_rep_id] 
        new_cluster_name=map_rep_seq_new_cluster[new_rep_id]
        pident = float(cells[2]) / 100
        if new_cluster_name in max_match_new_cluster_old_cluster:
            max_match_new_cluster_old_cluster[new_cluster_name][old_cluster_name]=pident
            if pident>max_match_new_cluster_old_cluster[new_cluster_name]['maxi']:
                max_match_new_cluster_old_cluster[new_cluster_name]['maxi']=pident
                max_match_new_cluster_old_cluster[new_cluster_name]['match']=old_cluster_name
        else:
            max_match_new_cluster_old_cluster[new_cluster_name]={}
            max_match_new_cluster_old_cluster[new_cluster_name]['maxi']=pident
            max_match_new_cluster_old_cluster[new_cluster_name]['match']=old_cluster_name
    #for all match cluster, merge:
    num_cluster=len(new_annotated_clusters)
    for new_c in max_match_new_cluster_old_cluster:
        old_annotated_clusters[max_match_new_cluster_old_cluster[new_c]['match']]=merge_2_annotated_clusters(old_annotated_clusters[max_match_new_cluster_old_cluster[new_c]['match']],new_annotated_clusters[new_c])
        del new_annotated_clusters[new_c]
    print("reduce "+str(num_cluster-len(new_annotated_clusters))+" after merge core clusters")
    num_cluster=len(new_annotated_clusters)
    remain_representative_new_cluster_seqs=os.path.join(temp_dir,"remain_representative_new_clusters.fasta")
    with open(new_unique_seqs,'r') as fi, open(remain_representative_new_cluster_seqs,'w') as fc:
        for seq in SeqIO.parse(fi,'fasta'):
            if seq.id in map_rep_seq_new_cluster and map_rep_seq_new_cluster[seq.id] in new_annotated_clusters:
                SeqIO.write(seq,fc,'fasta') 
    #grouping with old cloud cluster:
    #concat :
    concat_new_and_cloud_rep_cluster=os.path.join(temp_dir,"concat_new_and_cloud_rep_cluster.fasta")
    cmd=f'cat {representative_cloud_old_cluster_seqs} {remain_representative_new_cluster_seqs} > {concat_new_and_cloud_rep_cluster}'
    ret = run_command(cmd,timing_log)
    groups_representative_fasta, similar_groups = run_mmseq_with_map_similar_seqs(
        faa_file=concat_new_and_cloud_rep_cluster,    
        out_dir=temp_dir,      
        threads=threads,
        timing_log=timing_log,
        identity=0.50)
    
    pairs = []
    
   
    logger.info(f'similar_groups cloud clusters:{len(similar_groups)}')
    for g in similar_groups:
        if len(similar_groups[g])>1:
            items=[g]
            items.extend(similar_groups[g])
            items_from_A = [item for item in items if item in map_rep_seq_old_cluster.keys()]
            items_from_B = [item for item in items if item in map_rep_seq_new_cluster.keys()]
            for item_a in items_from_A:
                for item_b in items_from_B:
                    pairs.append((map_rep_seq_old_cluster[item_a], map_rep_seq_new_cluster[item_b]))
    logger.info(f'MPair cloud clusters:{len(pairs)}')
    for p in pairs:
        if p[0] in old_annotated_clusters and p[1] in new_annotated_clusters:
            old_annotated_clusters[p[0]]=merge_2_annotated_clusters(map_rep_seq_old_cluster[p[0]],new_annotated_clusters[p[1]])
            del new_annotated_clusters[p[1]]
    print("reduce "+str(num_cluster-len(new_annotated_clusters))+" after merge cloud clusters")
    json.dump(old_annotated_clusters, open(old_annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    old_annotated_clusters_file=merge_new_cluster_to_old_clusters(new_annotated_clusters,old_annotated_clusters_file)
    #json.dump(old_annotated_clusters, open(old_annotated_clusters_file, 'w'), indent=4, sort_keys=True)
    save_update_clusters(old_annotated_clusters_file,root_dir)
    del old_annotated_clusters
    elapsed = datetime.now() - starttime
    logger.info(f'Merge done, remain {len(new_annotated_clusters.keys())} clusters  -- time taken {str(elapsed)}')
    return old_annotated_clusters_file,new_annotated_clusters
def merge_2_annotated_clusters(old_cluster,new_cluster):
    list_unique_seq=[]
    if check_path(old_cluster['unique_seq']):
        list_unique_seq=read_array(old_cluster['unique_seq'])
    else:
        list_unique_seq=old_cluster['unique_seq']
    list_unique_seq.extend(new_cluster['unique_seq'])
    old_cluster['min_length']=min(old_cluster['min_length'],new_cluster['min_length'])
    old_cluster['max_length']=max(old_cluster['max_length'],new_cluster['max_length'])
    old_cluster['mean_length']=(old_cluster['mean_length']*old_cluster['size']+new_cluster['mean_length']*new_cluster['size'])/(old_cluster['size']+new_cluster['size'])

    old_cluster['size']=old_cluster['size']+new_cluster['size']
    old_cluster['unique_seq']=write_array(old_cluster['unique_seq'], list_unique_seq)
    return old_cluster
def reduce_unique_seqs_and_expand_seq_ids(clusters_file, old_groups_file,new_groups,old_seqs_fasta,new_seqs_fasta,root_dir,temp_dir,threads, timing_log):
    starttime = datetime.now()
    old_groups=json.load(open(old_groups_file, 'r'))
    
    combined_unique_seqs=os.path.join(temp_dir,'combined_unique_seqs.fasta')
    cmd=f'cat {old_seqs_fasta} {new_seqs_fasta} > {combined_unique_seqs}'
    ret = run_command(cmd,timing_log)
    unique_seqs_fasta, unique_groups = run_mmseq_unique_seqs(
        faa_file=combined_unique_seqs,     
        out_dir=temp_dir,      
        threads=threads,
        timing_log=timing_log)
    json.dump(unique_groups, open(os.path.join(temp_dir,'combined_group_unique_seqs.json'), 'w'), indent=4, sort_keys=True)
    merged_seq_in_groups={}
    del_unique=set()
    revert_map_seq_group={}
    for g in unique_groups:
        merged_seq_in_groups[g]=[g]
        if g in old_groups:
            merged_seq_in_groups[g].extend(read_array(old_groups[g])) 
        else:
            merged_seq_in_groups[g].extend(new_groups[g]) 
        for s in unique_groups[g]:
            merged_seq_in_groups[g].append(s)
            if s in old_groups:
                
                merged_seq_in_groups[g].extend(read_array(old_groups[s]))    
            else:
                merged_seq_in_groups[g].extend(new_groups[s])
            del_unique.add(s)
        #merged_seq_in_groups=list(set(merged_seq_in_groups))
    #count num of sq
    count_seq=0
    for g in merged_seq_in_groups:
        count_seq=count_seq+len(merged_seq_in_groups[g])
        for id in merged_seq_in_groups[g]:
            revert_map_seq_group[id]=g
    logger.info(f"count_seq in merged groups:{count_seq}")
    logger.info(f'merged_seq_in_groups  {len(merged_seq_in_groups)} clusters')
    
    clusters=json.load(open(clusters_file, 'r'))
    json.dump(clusters, open(clusters_file+".bak", 'w'), indent=4, sort_keys=True)
   
    list_del_cluster=[]
    set_used_seq=set()
    for c in  clusters:
        new_unique=[]
        new_ids=[]
        if check_path(clusters[c]['unique_seq']):
            list_unique_seqs=read_array(clusters[c]['unique_seq'])
        else:
            list_unique_seqs=clusters[c]['unique_seq']
        #list_unique_seqs=read_array(clusters[c]['unique_seq'])
        list_unique_seqs=list(set(list_unique_seqs))
        #print(c)
        for s in list_unique_seqs:
            #print(s)
            ns=revert_map_seq_group[s]
            if ns in merged_seq_in_groups:
                new_unique.append(ns)
                new_ids.extend(merged_seq_in_groups[ns]) 
                del merged_seq_in_groups[ns] 
            #new_unique.append(s)
        
        new_ids=list(set(new_ids))
        os.remove(clusters[c]['gene_id'])
        os.remove(clusters[c]['unique_seq'])
        if len(new_ids)==0 or len(new_unique)==0:
            list_del_cluster.append(c)
        clusters[c]['gene_id']=write_array(clusters[c]['gene_id'],new_ids)
        clusters[c]['unique_seq']=write_array(clusters[c]['unique_seq'],new_unique)
    logger.info(f'remove  {list_del_cluster} clusters')
    logger.info(f'merged_seq_in_groups  {len(merged_seq_in_groups)} clusters')
    for c in list_del_cluster:
        # list_seqs=read_array(clusters[c]['gene_id'])
        del clusters[c]
    save_unique_seqs(merged_seq_in_groups,root_dir)
    json.dump(clusters, open(clusters_file, 'w'), indent=4, sort_keys=True)
    elapsed = datetime.now() - starttime
    logger.info(f'reduce unique seq and expand ids  -- time taken {str(elapsed)}')
    #print(set(merged_seq_in_groups)-set_used_seq)
    return unique_seqs_fasta