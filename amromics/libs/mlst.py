# -*- coding: utf-8 -*-
"""
    Detect genes using blast


Revision history:
----------------
2020-02-13: Amromics created

"""
from __future__ import division, print_function, absolute_import
import subprocess
import os, shutil, glob
import re
def find_mlst(query_file,num_threads=1,blastdb='db/mlst/blast/mlst.fa',mlstdb='db/mlst/pubmlst',identity=95,mincov=10,minscore=50,scheme='',exclude='ecoli_2,abaumannii'):
    """
    Call blastn with params
    :param query_file (in fasta), db (blast indexed db), number of threads and identity
    :return: list BLASTFields objects
    """

    # run blastn
    cmd ='blastn -query {query} -db {db} -num_threads {threads} -ungapped -dust \
    no -word_size 32 -max_target_seqs 10000 -perc_identity {identity} -evalue 1E-20 \
    -outfmt \'6 sseqid slen length nident qseqid qstart qend qseq sstrand\' > temp.tab'.format(
        query=query_file,
        identity=identity,
        db=blastdb,
        threads=num_threads
    )
    # cmd ='any2fasta -q {query} | blastn -db {db} -num_threads {threads} -ungapped -dust \
    # no -word_size 32 -max_target_seqs 10000 -perc_identity {identity} -evalue 1E-20 \
    # -outfmt \'6 sseqid slen length nident qseqid qstart qend qseq sstrand\' > temp.tab'.format(
    #     query=query_file,
    #     identity=identity,
    #     db=db_file,
    #     threads=num_threads
    # )
    #print (cmd)
    os.system(cmd)
    #make exclude set from param
    exclude_set=set()
    tok=exclude.split(',')
    for t in tok:
        exclude_set.add(t.strip())
    #parse result
    f=open('temp.tab')
    line = f.readline()
    res={}
    countline=0
    res_e=0;
    res_c=0;
    res_h=0;
    while line:
        #result.append(line)
        #print (line)
        res_h=res_h+1
        countline=countline+1
        z=re.match(r"(\w+)\.(\w+)[_-](\d+)\t(\d+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)", line.strip())
        if z:
            res_c=res_c+1
            hit={'sch':z[1],'gene':z[2],'num':z[3],'hlen':z[4],'alen':z[5],'nident':z[6]\
            ,'qid':z[7],'qstart':z[8],'qend':z[9],'qseq':z[10],'sstrand':z[11]}
            #print(hit)
            #print (hit['nident'],hit['hlen'])
            if int(hit['nident'])/int(hit['hlen']) < mincov/100:
                line = f.readline()
                #print ('ignor < mincov',countline)
                #print (countline)
                res_e=res_e+1
                continue
            if not scheme=='' and not scheme==hit['sch']:
                line = f.readline()
                res_e=res_e+1
                #print ('ignor not in scheme',countline)
                continue
            if hit['sch'] in exclude_set :
                line = f.readline()
                res_e=res_e+1
                #print (hit['sch'],'ignor in exclue scheme',countline)
                continue

            if not hit['sch'] in res.keys():
                res[hit['sch']]={}
            if int(hit['hlen'])==int(hit['alen']) and int(hit['nident'])==int(hit['hlen']):

                if  hit['gene'] in res[hit['sch']].keys():
                    print('WARNING: found additional exact allele match',hit['sch'],'.',hit['gene'],'-',hit['num'])
                    res[hit['sch']][hit['gene']]=res[hit['sch']][hit['gene']]+','+hit['num']
                else:
                    print('Found exact allele match',hit['sch'],'.',hit['gene'],'-',hit['num'])
                    res[hit['sch']][hit['gene']]=hit['num']
            else:
                if int(hit['alen'])==int(hit['hlen']):
                    if not hit['gene'] in res[hit['sch']]:
                        res[hit['sch']][hit['gene']]='~'+str(hit['num'])
                else:
                    if not hit['gene'] in res[hit['sch']]:
                        res[hit['sch']][hit['gene']]=str(hit['num'])+'?'

            # find the signature with the fewest missing/approximate alleles



            #res.append(hit)

        line = f.readline()
    f.close()
    #if os.path.exists('temp.tab'):
    #    os.remove('temp.tab')

    #print ('res_h:',res_h,'res_e:', res_e,'res_c:',res_c)
    scheme_file_hash=load_scheme(mlstdb)
    list_sig=[]
    map_gene_name={}
    for k in res:
        STs,genes=load_genes_from_tabfile(scheme_file_hash[k])
        #print(STs)
        map_gene_name[k]=genes
        sig=get_signature(genes,res[k])
        print (sig)
        ST=get_sequence_type(STs,sig)
        print (ST)
        if not ST=='-' and '-' in sig:
                sig=sig.replace('-','0')
        nlocii=len(genes)
        score=nlocii
        # # novelish
        score=score-0.3*sig.count('~')
        # approx
        score=score-0.5*sig.count('?')
        # absent
        score=score-1.0*sig.count('-')
        score=int(score*90/nlocii)
        if not ST=='-':
            score=score+10

        if score>=minscore:
            if sig.count('~')>0:
                ST=get_sequence_type(STs,sig.replace('~',''))+'*'
            list_sig.append({'name':k,'ST':ST,'sig':sig,'score':score})
    #sort the list sig by score
    list_sorted=sorted(list_sig,key = lambda i: (i['score'], i['ST']),reverse=True)
    #return the best score
    #make output like mlst
    #sprint(list_sorted)
    ret={}
    if len(list_sorted)>0:
        ret['file']=os.path.basename(query_file)
        ret['scheme']=list_sorted[0]['name']
        ret['st']=str(list_sorted[0]['ST'])
        g=[]
        tok=list_sorted[0]['sig'].split('/')
        genes=map_gene_name[list_sorted[0]['name']]
        for i in range(0,len(genes)):
            g.append(genes[i]+'('+str(tok[i])+')')
        ret['profile']=g
        #return list_sorted[0]
    else:
        #return {'name':'-','ST':'-','sig':'-/-/-/-/-/-/-','score':0}
        ret={'file':os.path.basename(query_file),'scheme':'','st':'-','profile':[]}
    print(ret)
    return ret

def load_scheme(db):
    """
    load all name in mlst db in dict
    """
    scheme={}
    #read all name in mslt folder, each name save text file path with name in dict
    listOfDir = os.listdir(db)
    print(len(listOfDir))
    for entry in listOfDir:
        # Create full path
        fullPath = os.path.join(db, entry,entry+".txt")
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(os.path.join(db, entry)):
            scheme[entry.strip()]=fullPath
    return scheme
def load_genes_from_tabfile(tabfile):
    STs={}
    genes=[]
    f=open(tabfile)
    lines=f.readlines()
    #header
    headers=lines[0].split('\t')
    for h in headers:
        if re.search('^(ST|mlst_clade|clonal_complex|species|CC|Lineage)',h)==None:
            genes.append(h.strip())
    for i in range(1,len(lines)):
        cols=lines[i].split('\t')
        if cols[0].isdigit():
            sig=''
            for j in range(0,len(genes)):
                sig=sig+'/'+cols[j+1]
            sig=sig[1:].strip()
            STs[sig]=cols[0]
    return STs,genes
def get_signature(genes,query):
    sig=''
    for i in range(0,len(genes)):
        if genes[i] in query.keys():
            sig=sig+'/'+query[genes[i]]
        else:
            sig=sig+'/'+'-'
    sig=sig[1:]#remove first /
    return sig
def get_sequence_type(STs,sig):
    st='-'
    if sig in STs.keys():
        return STs[sig]
    return st

def setupdb(dbfile,dbtype=None):
    """
    make blast database from fasta file
    :param : fasta file (with folder is the name of db)
    :return:
    """
    #get name of db from file path:
    name=os.path.basename(os.path.dirname(os.path.dirname(dbfile)))
    #check type:
    if dbtype==None:
        dbtype=sequence_type(dbfile)
    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
        path=dbfile,
        name=name,
        type=dbtype

    )
    os.system(cmd)
