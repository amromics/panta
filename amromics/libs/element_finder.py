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
import sys
import argparse
def search_amr(sample,output,threads=1):
    db='db/ncbi/sequences'
    ret=blast(sample=sample,db=db,output=output,  threads=threads)
    export_file(sample,'ncbi',ret,output)
    return ret
def search_plasmid(sample,output,threads=1):
    db='db/plasmidfinder/sequences'
    ret=blast(sample=sample,db=db,output=output, identity=95,mincov=60, threads=threads)
    export_file(sample,'plasmidfinder',ret,output)
    return ret
def search_integron(sample,output,threads=1):
    db='db/integron/sequences'
    ret=blast(sample=sample,db=db,output=output,mincov=60,  threads=threads)
    export_file(sample,'integron',ret,output)
    return ret
def search_integrall(sample,output,threads=1):
    db='db/integrall/sequences'
    ret=blast(sample=sample,db=db,output=output,mincov=60,  threads=threads)
    export_file(sample,'integrall',ret,output)
    return ret
def search_prophage(sample,output,threads=1):
    db='db/prophage/sequences'
    ret=blast(sample=sample,db=db,output=output,mincov=60,  threads=threads,dbtype='prot')
    export_file(sample,'prophage',ret,output)
    return ret
def search_virulome(sample,output,threads=1):
    db='db/vfdb/sequences'
    ret=blast(sample=sample,db=db,output=output,mincov=60,  threads=threads)
    export_file(sample,'vfdb',ret,output)
    return ret
def export_file(sample,db,result,output):
    f=open(output,'w')
    f.write('FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n')
    for s in result:

        r={}
        r['FILE']=os.path.basename(sample)
        r['SEQUENCE']=s['qseqid']
        r['START']=s['qstart']
        r['END']=s['qend']
        #swap start and end when strain is minus:
        if s['sstrand']=='minus':
            t=s['sstart']
            s['sstart']=s['send']
            s['send']=t
            r['STRAND']='-'
        else:
            r['STRAND']='+'

        pccov = 100 * (int(s['length'])-int(s['gaps'])) / int(s['slen'])
        spl_ssqeid=s['sseqid'].split('~~~')
        if len(spl_ssqeid)>=4:
            r['GENE']=spl_ssqeid[1]
            r['DATABASE']=spl_ssqeid[0]
            r['ACCESSION']=spl_ssqeid[2]
            r['RESISTANCE']=spl_ssqeid[3]
        elif len(spl_ssqeid)>=3 :
            r['GENE']=spl_ssqeid[1]
            r['DATABASE']=spl_ssqeid[0]
            r['ACCESSION']=spl_ssqeid[2]
            r['RESISTANCE']=''
        else:
            r['GENE']=s['sseqid']
            r['DATABASE']=db
            r['ACCESSION']=''
            r['RESISTANCE']=''

        r['PRODUCT']=s['stitle']
        z=re.match('(^\S+\s+)',s['stitle'])
        if z:
            r['PRODUCT']=r['PRODUCT'].replace(z[1],'')
            #print(r['PRODUCT'])
        #$product =~ s/[,\t]//g;  # remove output separators
        #https://github.com/tseemann/abricate/issues/95
        #$product =~ s/^\S+\s+// if $product =~ m/$IDSEP/;
        r['COVERAGE']=str(s['sstart'])+'-'+s['send']+'/'+s['slen']
        #r['COVERAGE_MAP']=minimap(int(s['sstart']),int(s['send']),int(s['slen']),int(s['gapopen']))
        r['%COVERAGE']='{0:.3g}'.format(pccov)
        r['%IDENTITY']=s['pident']
        r['GAPS']=s['gapopen']+'/'+s['gaps']
        f.write(r['FILE']+'\t'+r['SEQUENCE']+'\t'+str(r['START'])+'\t'+str(r['END'])+'\t'+\
            r['STRAND']+'\t'+r['GENE']+'\t'+str(r['COVERAGE'])+'\t'+\
            str(r['GAPS'])+'\t'+str(r['%COVERAGE'])+'\t'+str(r['%IDENTITY'])+'\t'+r['DATABASE']+'\t'+r['ACCESSION']+\
            '\t'+r['PRODUCT']+'\t'+r['RESISTANCE']+'\n')
    f.close()
def blast(sample,db,output, identity=90, threads=1, mincov=0,dbtype='nucl'):
    """
    Call blastn with params
    :param query_file (in fasta), db (blast indexed db), number of threads and identity
    :return: list BLASTFields objects
    """
    #check db is indexed
    #dbfile=os.path.join(db_folder, 'sequences')




    # run blastn
    cmd ='blastn -query {query} -task blastn -dust no -perc_identity {identity} -db {db} -outfmt \'6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle\' -num_threads {threads} -evalue 1E-20 -culling_limit 1 > temp.tab'.format(

        query=sample,
        identity=identity,
        db=db,
        threads=threads

    )

    if dbtype=='prot':
        cmd ='blastp -query {query} -task blastp  -db {db} -outfmt \'6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle\' -num_threads {threads} -evalue 1E-20 > temp.tab'.format(

            query=sample,

            db=db,
            threads=threads

        )
    print(cmd)
    os.system(cmd)
    #parse result
    f=open('temp.tab')
    line = f.readline()
    result=[]
    while line:
        #result.append(line)
        t=line.strip().split('\t')
        blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
         'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
          'evalue':t[9], 'length':t[10], 'pident':t[11], 'gaps':t[12], 'gapopen':t[13],\
           'stitle':t[14]}
        result.append(blast_fields)
        line = f.readline()
    f.close()
    if os.path.exists('temp.tab'):
        os.remove('temp.tab')
    ret=[]
    for s in result:

        pccov = 100 * (int(s['length'])-int(s['gaps'])) / int(s['slen'])
        if pccov<=mincov:
            continue
        ret.append(s)

    return ret





def minimap(x,y,L,gap):
    width=15-(1 if gap else 0)
    broken=gap
    scale=L/width
    on='='
    off='.'
    x=int(x/scale)
    y=int(y/scale)
    L=int(L/scale)
    map=''
    #print (x,' ',y,' ',L,' ',scale,' ',width)
    for i in range(0,width):

        if i>=x and i<=y:
            map=map+on
        else:
            map=map+off
        if broken and i==int(width/2):
            map=map+'/'

    return map

def setupdb():
    """
    make blast database from fasta file in db folder,
    :param : fasta file (with folder'folder is the name of db and filename is 'sequences')
    :return:
    """
    #get name of db from file path:
    for root, dirs, files in os.walk('db'):
        for _file in files:
            if _file.endswith(('sequences')):
                name=os.path.basename(str(root))
                #print (name)
                seqfile=str(root)+'/'+_file
                dbtype=sequence_type(seqfile)
                cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                     path=seqfile,
                     name=name,
                     type=dbtype

                )
                print (cmd)
                os.system(cmd)

def sequence_type(dbfile):
    f=open(dbfile)
    line = f.readline()
    type='nucl'

    count=0
    while line:
        #result.append(line)
        if line.startswith('>') or line.startswith(';'):
            line = f.readline()
            continue

        else:
            if re.search(r"[^ATGCNatcgn-]", line.strip()):
                type='prot'
                break
        line = f.readline()
        count=count+1
        #only check top 100 lines
        if count>=100:
            break
    f.close()
    return type
